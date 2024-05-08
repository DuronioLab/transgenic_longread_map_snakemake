#!/bin/bash

# Assign command-line arguments to variables
input_fastq=$1
reads=$2
fqscreenPath=$3
threads=$4
fqscreenConf=$5
output_html=$6


total_reads=$(($(wc -l < $input_fastq ) / 4))

echo $total_reads

if (( $total_reads > $reads ))
then
  proportion=$(echo "$reads / $total_reads" | bc -l)
  proportion=$(printf "%.2f" $proportion)
else
  proportion=1.0
fi

echo $proportion

seqkit sample -p $proportion $input_fastq > concat_reduced.fastq

awk '{
    if (NR % 4 == 1 || NR % 4 == 3) {
        print $0
    } else {
        seq_length=length($0)
        start=(length / 2) - 50
        print substr($0, start, 100)
    }
}' concat_reduced.fastq > concat_reduced_100bp.fastq

$fqscreenPath --threads $threads --force --aligner bowtie2 --conf $fqscreenConf concat_reduced_100bp.fastq --outdir ./FQScreen

mv ./FQScreen/concat_reduced_100bp_screen.html $6
rm ./FQScreen/concat_reduced_100bp_screen.txt
rm -rf ./FQScreen
