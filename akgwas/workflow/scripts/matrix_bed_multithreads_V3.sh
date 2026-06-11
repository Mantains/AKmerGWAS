#!/bin/bash

## Input file can be a single matrix file or a stream from STDIN. 
## Output is a set of bed files in the specified output directory.
## Example usage:
## 1. From a file:
##    ./matrix_bed_multithreads_V3.sh \
##      -fof sample.fof \
##      -file input_matrix.txt \
##      -o output_dir \
##      -m 0.05 \
##      -b 2G \
##      -t 8
## 2. From a stream:
##    ./matrix_bed_multithreads_V3.sh \
##      -fof sample.fof \
##      -file <(pigz -dc Path/to/matrices/Matrix_part_*.gz) \
##      -o output_dir \
##      -m 0.05 \
##      -b 2G \
##      -t 8

set -euo pipefail

thread=8
blocksize=2G
fof_file=""
input_file=""
output_dir=""
maf=0.05

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -fof) fof_file="$2"; shift ;;
        -file) input_file="$2"; shift ;;
        -o) output_dir="$2"; shift ;;
        -m) maf="$2"; shift ;;
        -b) blocksize="$2"; shift ;;
        -t) thread="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Ensure output directory exists
mkdir -p "$output_dir/bed"
log_file="$output_dir/kmatrix_to_bed.log"

{
start_time=$(date +%s)
echo "Start maf filter and BED conversion: $(date)"

# Calculate sample counts and the limit for filtering
# Note: The original logic filters out k-mers present in > (N - MAF*N) samples
Num_samples=$(wc -l < "$fof_file")
lower_limit=$(awk -v n="$Num_samples" -v m="$maf" 'BEGIN{print n * m}')
upper_limit=$(awk -v n="$Num_samples" -v m="$maf" 'BEGIN{print n * (1 - m)}')

echo "Number of samples: $Num_samples"
echo "Filtering criteria: $lower_limit < presence_count < $upper_limit"

# Prepare .tfam (Standard format for both PLINK 1.9 and 2.0)
# PLINK2 reading TPED requires a TFAM file
cut -d ' ' -f 1 "$fof_file" | awk '{print "0", $1, "0 0 0 0"}' > "$output_dir/bed/common.tfam"

echo "Step: Processing stream and converting to BED chunks via pipe..."

# Parallel processing: 
# - Uses --pipe to distribute the zcat stream
# - awk transforms 0/1 matrix to TPED format on-the-fly
# - plink reads from /dev/stdin and writes directly to binary BED format

parallel --pipe -j $thread --block $blocksize --max-lines 0 '
parallel --pipe -j "$thread" --block "$blocksize" --max-lines 0 '
    # Create a unique temp file name for this worker
    TMP_TPED="tmp_worker_{#}.tped"

    # Step A: awk filters and writes to the local temp file
    awk -v low="'"$lower_limit"'" -v up="'"$upper_limit"'" '\''
    {
        count = 0
        for (i=2; i<=NF; i++) { if ($i == 1) count++ }
        
        if (count > low && count < up) {
            out = "1 " $1 " 0 0"
            for (i=2; i<=NF; i++) {
                out = out ($i == 1 ? " 2 2" : " 1 1")
            }
            print out
        }
    }'\'' > "$TMP_TPED"

    # Step B: Run PLINK2 if the temp file is not empty
    if [ -s "$TMP_TPED" ]; then
        plink2 --tped "$TMP_TPED" \
               --tfam "'"$output_dir"'/bed/common.tfam" \
               --make-bed \
               --out "'"$output_dir"'/bed/matrix_{#}" \
               --threads 1 \
               --silent
    fi

    # Step C: Immediate cleanup of the temp file
    rm -f "$TMP_TPED"
' < "$input_file"

# Validation: Count records in the generated .bim files
echo "Step: Summarizing filtered results..."
find "$output_dir/bed/" -name "matrix_*.bim" | sort -V | while read -r bim; do
    count=$(wc -l < "$bim")
    echo "Matrix $(basename "$bim" .bim) variants: $count"
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo -e "BED build done. End time: $(date)\nTotal time used: ${elapsed_time} seconds"

} 2>&1 | tee -a "$log_file"