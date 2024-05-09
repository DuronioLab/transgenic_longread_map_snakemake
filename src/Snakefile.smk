import json
import os
import sys
import preProcessSampleConfig as pre

#import the config.json file as config
with open('./src/config.json') as json_data:
    config = json.load(json_data)

file_info_path = config['sampleInfo']
basename_columns = config['baseNameCols']

REFGENOME = config['refGenome']
DEFAULTGENOME = config['defaultGenome']
#chromSize_Path = config['genome'][DEFAULTGENOME]['chromSize']
#genomeSize = config['genome'][DEFAULTGENOME]['genomeSize']

modules = config['module']

#########
# Validation

if not os.path.exists(file_info_path):
    sys.exit('Error: {name} does not exist. Be sure to set `sampleInfo` in config.json.'.format(name=file_info_path))

if type(REFGENOME) is not str:
    sys.exit('Error: refGenome must be a string. Currently set to: {}. Double check `refGenome` in config.json.'.format(REFGENOME))

#########
# Generating sampleSheet outputs

if type(REFGENOME) is str:
    REFGENOME = [REFGENOME]
if type(DEFAULTGENOME) is str:
    DEFAULTGENOME = [DEFAULTGENOME]

genomeList = REFGENOME + DEFAULTGENOME
sampleSheet = pre.makeSampleSheets(file_info_path,basename_columns,delim='-')

#For each unique sampleName, add a column for {sampleName}_concat.fastq
sampleSheet['concat'] = expand("{sample}_concat",sample=sampleSheet.sampleName)
sampleSheet['concat_fastq'] = expand("Fastq/{sample}_concat.fastq",sample=sampleSheet.sampleName)

#Add columns to sample sheet for each genome
for genome in genomeList:
    genome_bam = genome + '_bam'
    sampleSheet[
        genome_bam] = expand("Alignment/{concat_sample}_{genome}.{fileType}",concat_sample=sampleSheet.concat,genome=genome,fileType='bam')



#save sampleSheet as a text file in the working directory
sampleSheet.to_csv('sampleSheet.tsv',sep="\t",index=False)


######################
# Begin the pipeline #
######################


rule all:
    input:
        expand("Fastq/{sample}_concat.fastq",sample=set(sampleSheet.sampleName)),
        expand("Stats/NanoPlot-report_{concat_sample}.html",concat_sample=set(sampleSheet.concat)),
        expand("Stats/{concat_sample}_screen.html",concat_sample=set(sampleSheet.concat)),
        expand("Alignment/{concat_sample}_{genome}.{fileType}",concat_sample=set(sampleSheet.concat),genome=genomeList,fileType=[
            'bam', 'sam']),
        expand("Alignment/{concat_sample}_{genome}.bam.bai",concat_sample=set(sampleSheet.concat),genome=genomeList),
        expand("Stats/{concat_sample}_{genome}readDepth.pdf",concat_sample=set(sampleSheet.concat),genome=DEFAULTGENOME),
        expand("Medaka/{concat_sample}_consensus.bam", concat_sample=set(sampleSheet.concat)),
        expand("Medaka/{concat_sample}_consensus.bam.bai", concat_sample=set(sampleSheet.concat))

#expand("Medaka/{concat_sample}_consensus.bam", concat_sample=sampleSheet.concat)

#Snakemake rule that concatenates the fastq files for each sample found in each sampleDirectory of the sampleSheet

rule concatFastq:
    input:
        directory=expand("{sampleDir}",sampleDir=sampleSheet.sampleDirectory)
    output:
        fastq="Fastq/{sample}_concat.fastq"
    shell:
        """
        for dir in {input.directory}
        do
            zcat $dir/*.fastq.gz >> {output.fastq}
        done
        """

rule nanoplot:
    input:
        fastq=expand("Fastq/{concat_sample}.fastq",concat_sample=set(sampleSheet.concat))
    output:
        html=expand("Stats/NanoPlot-report_{concat_sample}.html", concat_sample=set(sampleSheet.concat))
    envmodules:
        modules['nanopackVer']
    shell:
        """
        NanoPlot --fastq_rich {input.fastq} --N50 -o NanoPlot
        mv NanoPlot/NanoPlot-report.html {output.html}
        rm -R -f NanoPlot
        """

rule fqscreen:
    input:
        fastq=expand("Fastq/{concat_sample}.fastq",concat_sample=set(sampleSheet.concat))
    output:
        html="Stats/{concat_sample}_screen.html"
    envmodules:
        modules['seqkitVer']
    params:
        fqscreenPath=modules['fqscreenPath'],
        fqscreenConf=modules['fqscreenConf'],
        reads=100000
    threads: 4
    shell:
        """
        bash ./src/fqscreen.sh {input.fastq} {params.reads} {params.fqscreenPath} {threads} {params.fqscreenConf} {output.html}

        """

rule align:
    input:
        fastq=expand("Fastq/{concat_sample}.fastq",concat_sample=set(sampleSheet.concat))
    output:
        sam=expand("Alignment/{concat_sample}_{genome}.sam",concat_sample=set(sampleSheet.concat),genome=DEFAULTGENOME),
        bam=expand("Alignment/{concat_sample}_{genome}.bam",concat_sample=set(sampleSheet.concat),genome=DEFAULTGENOME),
        bamIndex=expand("Alignment/{concat_sample}_{genome}.bam.bai",concat_sample=set(sampleSheet.concat),genome=DEFAULTGENOME)
    envmodules:
        modules['minimap2Ver'],
        modules['samtoolsVer']
    threads: 4
    params:
        genome=[config['genome'][genome]['fasta'] for genome in DEFAULTGENOME]
    shell:
        """
        minimap2 --secondary=no --sam-hit-only -ax map-ont {params.genome} {input.fastq} > {output.sam}
        samtools sort {output.sam} -o {output.bam}
        samtools index {output.bam}
        """

rule chrom_graph:
    input:
        expand("Alignment/{concat_sample}_{genome}.bam",concat_sample=set(sampleSheet.concat),genome=DEFAULTGENOME)
    output:
        expand("Stats/{concat_sample}_{genome}readDepth.pdf",concat_sample=set(sampleSheet.concat),genome=DEFAULTGENOME)
    envmodules:
        modules['samtoolsVer']
    shell:
        """
        bash ./src/chrom_graph.sh {input} {output}
        """

rule query_align:
    input:
        fastq=expand("Fastq/{concat_sample}.fastq",concat_sample=set(sampleSheet.concat)),
        fasta=[config['genome'][genome]['fasta'] for genome in REFGENOME]
    output:
        sam=expand("Alignment/{concat_sample}_{genome}.sam",concat_sample=set(sampleSheet.concat),genome=REFGENOME),
        bam=expand("Alignment/{concat_sample}_{genome}.bam",concat_sample=set(sampleSheet.concat),genome=REFGENOME),
        bamIndex=expand("Alignment/{concat_sample}_{genome}.bam.bai",concat_sample=set(sampleSheet.concat),genome=REFGENOME),
        filteredFasta=expand("Alignment/{concat_sample}_{genome}_filtered.fasta",concat_sample=set(sampleSheet.concat),genome=REFGENOME)
    envmodules:
        modules['samtoolsVer'],
        modules['blatVer'],
        modules['seqkitVer'],
        modules['seqtkVer'],
        modules['minimap2Ver']
    params:
        query="query.fasta"
    shell:
        """
        bash ./src/query_align.sh {input.fastq} {input.fasta} {output.sam} {output.bam} {params.query}
        """

rule consensus:
    input:
        fasta=expand("Alignment/{concat_sample}_{genome}_filtered.fasta",concat_sample=set(sampleSheet.concat),genome=REFGENOME)
    output:
        bam=expand("Medaka/{concat_sample}_consensus.bam",concat_sample=set(sampleSheet.concat)),
        index=expand("Medaka/{concat_sample}_consensus.bam.bai",concat_sample=set(sampleSheet.concat))
    params:
        genome=[config['genome'][genome]['fasta'] for genome in REFGENOME]
    envmodules:
        modules['medakaVer']
    shell:
        """
        medaka_consensus -i {input.fasta} -d {params.genome} -o temp_medaka
        mv ./temp_medaka/calls_to_draft.bam {output.bam}
        mv ./temp_medaka/calls_to_draft.bam.bai {output.index}
        rm ./temp_medaka/consensus_probs.hdf
        rm -R -f temp_medaka
        """
