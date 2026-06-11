#!/bin/bash
## Build kinship and pca from bed files, it will find all matrix_*.bed files in the
## input directory, sample variants, and merge them to create a representative subset 
## for PCA and kinship calculation. 
## Then it will perform LD pruning and calculate PCA and kinship using plink and gemma.

set -euo pipefail

bedinpath="./bed"
bedoutpath="./QK"
sampling_rate=0.01
thread=8

# Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -bedpath) bedinpath="$2"; shift ;;
        -o) bedoutpath="$2"; shift ;;
        -rate|-bfb) sampling_rate="$2"; shift ;;
        -t) thread="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if ! awk -v rate="$sampling_rate" 'BEGIN{exit !(rate > 0 && rate <= 1)}'; then
    echo "sampling rate must be in the range (0, 1]"
    exit 1
fi

mkdir -p "$bedoutpath"
tempdir="$bedoutpath/QK_temp"
mkdir -p "$tempdir"

# Sample and merge bed files to create a representative subset for PCA and kinship calculation.
mapfile -t bed_files < <(find "$bedinpath" -maxdepth 1 -name 'matrix_*.bed' | sort -V)
echo "Starting sampling with $thread threads..."

Numberbeds=${#bed_files[@]}
if [ "$Numberbeds" -lt 1 ]; then
    echo "No bed files found in: $bedinpath"
    exit 1
fi

# Sample each bed file directly with GNU parallel, then merge the sampled outputs.
export sampling_rate tempdir
printf '%s\n' "${bed_files[@]%.bed}" | parallel -j "$thread" '
    plink2 --allow-extra-chr \
        --bfile {} \
        --thin "$sampling_rate" \
        --make-bed \
        --out "$tempdir"/{/}_sampled \
        --silent
'

mergelist="$bedoutpath/mergelist.txt"
echo -n "" > "$mergelist"

sampled_prefixes=()
for bed_prefix in "${bed_files[@]%.bed}"; do
    bed_base="$(basename "$bed_prefix")"
    sampled_prefix="$tempdir/${bed_base}_sampled"
    if [ -s "$sampled_prefix.bed" ]; then
        printf '%s %s %s\n' "$sampled_prefix.bed" "$sampled_prefix.bim" "$sampled_prefix.fam" >> "$mergelist"
        sampled_prefixes+=("$sampled_prefix")
    fi
done

if [ "${#sampled_prefixes[@]}" -lt 1 ]; then
    echo "Sampling rate produced no variants; increase the rate and try again."
    exit 1
fi

# Prepare the sampled merged bed file for PCA and kinship calculation.
if [ "${#sampled_prefixes[@]}" -eq 1 ]; then
    plink --bfile "${sampled_prefixes[0]}" --make-bed --out "$bedoutpath/QK_bed_origin"
else
    plink --merge-list "$mergelist" --make-bed --out "$bedoutpath/QK_bed_origin"
fi

# plink may change the sample order during merging, we need to sort the fam file to make sure the order is correct for gemma
## Extract sample order from 1st bed file
awk '{print $1, $2}' "${sampled_prefixes[0]}.fam" > "$tempdir/sample_order.txt"
## Sort the QK_bed_origin according to the sample_order.txt
plink --bfile "$bedoutpath/QK_bed_origin" \
      --keep "$tempdir/sample_order.txt" \
      --make-bed \
      --out "$bedoutpath/QK_bed" \
      --silent

# No need for pruning 

#PCA (Q)

plink --bfile "$bedoutpath/QK_bed" --pca 5 --out "$bedoutpath/QK_pca"
awk '{print $3, $4, $5}' "$bedoutpath/QK_pca.eigenvec" > "$bedoutpath/QK_pca1-3.txt"

#Kinship (K)
#give false pheno to makesure gemma accept correct format
awk '{$6=1; print}' "$bedoutpath/QK_bed.fam" > "$bedoutpath/QK_bed.fam.tmp" && \
mv "$bedoutpath/QK_bed.fam.tmp" "$bedoutpath/QK_bed.fam"

cd "$bedoutpath/"
gemma -bfile "$bedoutpath/QK_bed" \
      -gk -miss 1.0 -maf 0.0 -r2 1.0 \
      -o QK_kinship

mv "$bedoutpath/output"/* "$bedoutpath/"
rm -rf "$bedoutpath/output/"
rm -rf "$tempdir"

#OUTPUT in $bedoutpath/