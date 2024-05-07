import json
import os
import sys
from src import preProcessSampleConfig as pre
from snakemake.io import expand

#import the config.json file as config
with open('./src/config.json') as json_data:
    config = json.load(json_data)

file_info_path = config['sampleInfo']
basename_columns = config['baseNameCols']

REFGENOME = config['refGenome']
DEFAULTGENOME = config['defaultGenome']
chromSize_Path = config['genome'][DEFAULTGENOME]['chromSize']
genomeSize = config['genome'][DEFAULTGENOME]['genomeSize']

modules = config['module']

#########
# Validation

if not os.path.exists(file_info_path):
    sys.exit('Error: {name} does not exist. Be sure to set `sampleInfo` in config.json.'.format(name = file_info_path))

if type(REFGENOME) is not str:
    sys.exit('Error: refGenome must be a string. Currently set to: {}. Double check `refGenome` in config.json.'.format(REFGENOME))

#########
# Generating sampleSheet outputs

if type(REFGENOME) is str:
    REFGENOME = [REFGENOME]
if type(DEFAULTGENOME) is str:
    DEFAULTGENOME = [DEFAULTGENOME]

genomeList = REFGENOME + DEFAULTGENOME
sampleSheet = pre.makeSampleSheets(file_info_path, basename_columns, delim='-')

#For each unique sampleName, add a column for {sampleName}_concat.fastq
sampleSheet['concat'] = expand("{sample}_concat", sample=sampleSheet.sampleName)
sampleSheet['concat_fastq'] = expand("Fastq/{sample}_concat.fastq", sample=sampleSheet.sampleName)

#Add columns to sample sheet for each genome
for genome in genomeList:
    genome_bam = genome + '_bam'
    sampleSheet[genome_bam] = expand("Bam/{concat_sample}_{genome}.{fileType}", concat_sample=sampleSheet.concat, genome=genome, fileType='bam')

#save sampleSheet as a text file in the working directory
sampleSheet.to_csv('sampleSheet.tsv', sep="\t", index=False)

######################
# Begin the pipeline #
######################


rule all:
    input:
        expand("Fastq/{sample}_concat.fastq", sample=sampleSheet.sampleName),
        expand("NanoPlot/NanoPlot-report_{concat_sample}.html", concat_sample=sampleSheet.concat),
        expand("FQScreen/{concat_sample}_screen.html", concat_sample=sampleSheet.concat),
        expand("Bam/{concat_sample}_{genome}.{fileType}", concat_sample=sampleSheet.concat, genome=genomeList, fileType=['bam','sam']),
        expand("Bam/{concat_sample}_{genome}.bam.bai", concat_sample=sampleSheet.concat, genome=genomeList),
        expand("Stats/{concat_sample}_readDepth.pdf", concat_sample=sampleSheet.concat),
        expand("Medaka/{concat_sample}_consensus.bam", concat_sample=sampleSheet.concat)

#Snakemake rule that concatenates the fastq files for each sample found in each sampleDirectory of the sampleSheet
rule concatFastq:
    input:
        expand("{sampleDir}/*.fastq.gz", sampleDir=sampleSheet.sampleDir)
    output:
        fastq = "Fastq/{sample}_concat.fastq"
    shell:
        """
        zcat {input} > {output.fastq}
        """

rule nanoplot:
    input:
        fastq = "Fastq/{concat_sample}_concat.fastq"
    output:
        html = "NanoPlot/NanoPlot-report_{concat_sample}.html"
    envmodules:
        modules['nanopackVer']
    shell:
        """
        NanoPlot --fastq_rich {input.fastq} --N50 -o {output.html}
        """

rule fqscreen:
    input:
        fastq = "Fastq/{concat_sample}_concat.fastq"
    output:
        html = "FQScreen/{concat_sample}_screen.html"
    envmodules:
        modules['seqkitVer']
    params:
        fqscreenPath = modules['fqscreenPath'],
        fqscreenConf = modules['fqscreenConf'],
        reads = 100000
    threads: 4
    shell:
        """
        total_reads=$(($(wc -l < {input}) / 4))

        if (( total_reads > {reads} ))
        then
          proportion=$(echo "{reads} / $total_reads" | bc -l)
          proportion=$(printf "%.2f" $proportion)
        else
          proportion=1.0
        fi
        
        seqkit sample -p $proportion concat.fastq > concat_reduced.fastq
        
        awk '{
            if (NR % 4 == 1 || NR % 4 == 3) {
                print $0
            } else {
                seq_length=length($0)
                start=(length / 2) - 50
                print substr($0, start, 100)
            }
        }' concat_reduced.fastq > concat_reduced_100bp.fastq

        {params.fqscreenPath} --threads {threads} --force --aligner bowtie2 --conf {params.fqscreenConf} concat_reduced_100bp.fastq --outdir ./FQScreen
        """

rule align:
    input:
        fastq = "Fastq/{concat_sample}_concat.fastq"
    output:
        sam = expand("Bam/{concat_sample}_{genome}.sam", concat_sample=sampleSheet.concat, genome=DEFAULTGENOME),
        bam = expand("Bam/{concat_sample}_{genome}.bam", concat_sample=sampleSheet.concat, genome=DEFAULTGENOME),
        bamIndex = expand("Bam/{concat_sample}_{genome}.bam.bai", concat_sample=sampleSheet.concat, genome=DEFAULTGENOME),
        chrGraph = "Stats/{concat_sample}_readDepth.pdf"
    envmodules:
        modules['minimap2Ver'],
        modules['samtoolsVer']
    threads: 4
    params:
        genome = DEFAULTGENOME
    shell:
        """
        minimap2 --secondary=no --sam-hit-only -ax map-ont {genome} {input.fastq} > {output.sam}
        samtools sort {output.sam} -o {output.bam}
        samtools index {output.bam}
        samtools idxstats {output.bam} | cut -f 1,3 | grep -v "\*" > output.txt
        
        echo "set terminal pdf size 8,6" > plot.gp
        echo "set output 'Stats/{concat_sample}_readDepth.pdf'" >> plot.gp
        echo "set style data histograms" >> plot.gp
        echo "set style fill solid 1.0 border -1" >> plot.gp
        echo "set xtics rotate by -45" >> plot.gp
        echo "plot 'output.txt' using 2:xtic(1) title 'Number of reads aligned to each chromosome'" >> plot.gp
        gnuplot plot.gp
        """

rule query_align:
    input:
        fastq = "Fastq/{concat_sample}_concat.fastq",
        fasta = config['genome'][REFGENOME]['fasta']
    output:
        sam = expand("Bam/{concat_sample}_{genome}.sam", concat_sample=sampleSheet.concat, genome=REFGENOME),
        bam = expand("Bam/{concat_sample}_{genome}.bam", concat_sample=sampleSheet.concat, genome=REFGENOME),
        bamIndex = expand("Bam/{concat_sample}_{genome}.bam.bai", concat_sample=sampleSheet.concat, genome=REFGENOME)
    envmodules:
        modules['samtoolsVer'],
        modules['blatVer'],
        modules['seqkitVer'],
        modules['seqtkVer'],
        modules['minimap2Ver']
    params:
        genome = REFGENOME
    shell:
        """
        seqkit replace --quiet -p .+ -r "seq_{nr}" {input.fastq} > renamed_reads.fastq
        seqkit seq --quiet --min-len $length renamed_reads.fastq -o out.fastq
        seqkit fq2fa --quiet out.fastq -o out.fasta
        
        split --lines=40000 -d --additional-suffix=".fasta" out.fasta out.

        touch all_results.psl

        for fasta_file in out.*.fasta
        do
            blat "$fasta_file" query.fasta -oneOff=3 -noHead output.psl
            cat output.psl >> all_results.psl
        done
        
        awk 'NR > 6 print $14' all_results.psl > temp_names.txt
        sort temp_names.txt | uniq > blat_names.txt

        readarray -t File < blat_names.txt
        touch filtered_reads.fasta
        for f in "$File[@]"
        do
          echo $f > fname.txt
          seqkit grep --quiet -n -f fname.txt out.fasta >> filtered_reads.fasta
        done

        seqtk seq -F '#' filtered_reads.fasta > $ref_basename_filtered.fastq

        minimap2 --secondary=no --sam-hit-only -ax map-ont {input.fasta} ./$ref_basename_filtered.fastq > {output.sam}

        samtools sort ./$ref_basename_mapped.sam -o ./$ref_basename_mapped.bam
        samtools view -bq 1 ./$ref_basename_mapped.bam > output.bam
        samtools index output.bam
        
        """