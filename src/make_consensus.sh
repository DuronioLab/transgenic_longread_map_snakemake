#!/bin/bash

# Assign command-line arguments to variables
input_fasta=$1
genome=$2
output_bam=$3
output_index=$4

echo "Running medaka_consensus..."
echo $input_fasta
echo $genome
echo $output_bam
echo $output_index

medaka_consensus -i $input_fasta -d $genome -o temp_medaka
mv ./temp_medaka/calls_to_draft.bam $output_bam
mv ./temp_medaka/calls_to_draft.bam.bai $output_index
rm ./temp_medaka/consensus_probs.hdf
rm -R -f temp_medaka