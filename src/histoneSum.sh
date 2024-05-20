#!/bin/bash

input_fastq=$1
array_fasta=$2
output_sam=$3
output_bam=$4

#Script that takes a 1xDWT file and uses blat to find all matches from the input_fasta. Then it will make new reads from each match from input_fasta.
#The new reads will be written to a new fasta file.

seqkit replace --quiet -p .+ -r "seq_{nr}" $input_fastq > renamed_reads.fastq
seqkit seq --quiet --min-len 5000 renamed_reads.fastq -o out.fastq
seqkit fq2fa --quiet out.fastq -o out.fasta

split --lines=40000 -d --additional-suffix=".fasta" out.fasta out.

touch subreads.fasta

for fasta_file in out.*.fasta
do
    blat "$fasta_file" $array_fasta -oneOff=3 -noHead output.psl
    samtools faidx "$fasta_file"
    counter=1
    while IFS= read -r line
    do
      read_name=$(echo $line | awk '{print $16}')
      start=$(echo $line | awk '{print $17}')
      end=$(echo $line | awk '{print $18}')

      new_read_name="${read_name}_${counter}"

      sequence=$(samtools faidx "$fasta_file" "$read_name:$start-$end")

      echo ">$new_read_name" >> subreads.fasta
      echo "$sequence" >> subreads.fasta
      ((counter++))
    done < output.psl
done

seqtk seq -F '#' subreads.fasta > subreads.fastq

minimap2 --secondary=no --sam-hit-only -ax map-ont $array_fasta subreads.fastq > $output_sam

samtools sort $output_sam -o mapped.bam
samtools view -bq 1 mapped.bam > $output_bam
samtools index $output_bam