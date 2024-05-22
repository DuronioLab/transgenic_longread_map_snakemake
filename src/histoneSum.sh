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

query_length=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next; } { seqlen = seqlen + length($0)}END{print seqlen}' "$array_fasta")

seqkit subseq --region 1:50 $array_fasta > start.fasta

samtools faidx $array_fasta

#touch repeats.txt
touch subreads.fasta

for fasta_file in out.*.fasta
do
  blat "$fasta_file" start.fasta -oneOff=3 -noHead output.psl
  awk -v len=$query_length 'BEGIN {OFS="\t"} {print $14, $16, $16+len}' output.psl > repeats.txt
  while read -r name start end
  do
    echo ${name} ${start} ${end}
    samtools faidx "$fasta_file" ${name}:${start}-${end} >> subreads.fasta

  done < repeats.txt
done

seqtk seq -F '#' subreads.fasta > subreads.fastq

##TODO check the parameters to ensure that minimap2 is aligning correctly. It has a tendency to have a short (150 bp) alignment as the primary alignment and the full 5kb alignment as secondary.
minimap2 -m 500 -ax map-ont $array_fasta subreads.fastq > $output_sam

samtools sort $output_sam -o mapped.bam
samtools view -bq 1 mapped.bam > $output_bam
samtools index $output_bam

rm out.*
rm mapped.bam
rm subreads.fast*
rm repeats.txt
rm output.psl
rm start.fasta
rm renamed_reads.fastq