#!/bin/bash

# Assign command-line arguments to variables
input_fastq=$1
fasta=$2
output_sam=$3
output_bam=$4
query=$5
filtered_fasta=$6

seqkit replace --quiet -p .+ -r "seq_{nr}" $input_fastq > renamed_reads.fastq
seqkit seq --quiet --min-len 5000 renamed_reads.fastq -o out.fastq
seqkit fq2fa --quiet out.fastq -o out.fasta

split --lines=40000 -d --additional-suffix=".fasta" out.fasta out.

touch all_results.psl

for fasta_file in out.*.fasta
do
    blat "$fasta_file" $query -oneOff=3 -noHead output.psl
    cat output.psl >> all_results.psl
done

awk 'NR > 6 {print $14}' all_results.psl > temp_names.txt
sort temp_names.txt | uniq > blat_names.txt

readarray -t File < blat_names.txt
touch filtered_reads.fasta
for f in "${File[@]}"
do
  echo $f > fname.txt
  seqkit grep --quiet -n -f fname.txt out.fasta >> filtered_reads.fasta
done

seqtk seq -F '#' filtered_reads.fasta > filtered_reads.fastq

minimap2 --secondary=no --sam-hit-only -ax map-ont $fasta filtered_reads.fastq > $output_sam

samtools sort $output_sam -o mapped.bam
samtools view -bq 1 mapped.bam > $output_bam
samtools index $output_bam

mv filtered_reads.fasta $filtered_fasta

rm out.*
rm mapped.bam
rm filtered_reads.fastq
rm fname.txt
rm temp_names.txt
rm all_results.psl
rm blat_names.txt
rm output.psl
rm renamed_reads.fastq