#!/bin/bash

INPUT_DIR='/xdisk/masel/uliseshmc/EstimatingUd/Gavin_apes114/intergenicAR/mafft_aligned/unaligned'
OUTPUT_DIR='/xdisk/masel/uliseshmc/EstimatingUd/Gavin_apes114/intergenicAR/mafft_aligned/interative_refined'
TREE_DIR='/xdisk/masel/uliseshmc/EstimatingUd/Gavin_apes114'

module load mafft

# Check arguments
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$TREE_DIR"  ]]; then
    echo "Usage: $0 input_folder output_folder"
    exit 1
fi

# Check input folder exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input folder '$INPUT_DIR' not found."
    exit 1
fi

# Ensure output directory exists
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Directory '$DIR' not found."
    exit 1
fi


# Loop through FASTA files
shopt -s nullglob
for f in "$INPUT_DIR"/* "$INPUT_DIR"/*.fasta; do
    echo "Aligning: $f"

    # Build output filename in OUTPUT_DIR
    base=$(basename "$f")
    out="$OUTPUT_DIR/${base%.*}_aligned.fa"

    mafft --auto --treein "$TREE_DIR"/tree.mafft "$f" > "$out"
    #use the next mafft command instead for iterative refined search  https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html#GLE
    #mafft --auto --maxiterate 1000 --treein "$TREE_DIR"/tree.mafft "$f" > "$out"
done




