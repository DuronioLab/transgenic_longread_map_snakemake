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

modules = config['module']

##############
# Validation #
##############

# Check whether important files/names are specified in the config.json file. Then check if the files exits.
# if not, exit the program.

if not os.path.exists(file_info_path):
    sys.exit('\nError: {name} does not exist. Be sure to set `sampleInfo` in config.json.\n'.format(name=file_info_path))

if type(REFGENOME) is not str:
    sys.exit('\nError: refGenome must be a string. Currently set to: {}. Double check `refGenome` in config.json.\n'.format(REFGENOME))

if isinstance(REFGENOME, list):
    for genome in REFGENOME:
        if genome not in config['genome']:
            sys.exit('\nError: Your `refGenome` {name} is not found as an entry in the `genome` section of config.json.\n'.format(name=REFGENOME))
else:
    if REFGENOME not in config['genome']:
        sys.exit('\nError: Your `refGenome` {name} is not found as an entry in the `genome` section of config.json.\n'.format(name=REFGENOME))

if isinstance(DEFAULTGENOME, list):
    for genome in DEFAULTGENOME:
        if genome not in config['genome']:
            sys.exit('\nError: Your `defaultGenome` is not found as an entry in the `genome` section of config.json.\n')
else:
    if DEFAULTGENOME not in config['genome']:
        sys.exit('\nError: Your `defaultGenome` is not found as an entry in the `genome` section of config.json.\n')

if not os.path.exists(config['genome'][DEFAULTGENOME]['fasta']):
    sys.exit('\nError: defaultGenome FASTA {name} does not exist. Be sure to set `fasta` in config.json.\n'.format(name=config['genome'][DEFAULTGENOME]['fasta']))

if not os.path.exists(config['genome'][REFGENOME]['fasta']):
    sys.exit('\nError: The refGenome FASTA {name} does not exist. Be sure to set `fasta` in config.json.\n'.format(name=config['genome'][REFGENOME]['fasta']))

if not os.path.exists('query.fa') and not os.path.exists('query.fasta'):
    sys.exit('\nError: query.fa or query.fasta does not exist. Please provide a query fasta file in the working directory.\n')

##################################
# Generating sampleSheet outputs #
##################################

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

# Snakemake rule that runs all rules and specified desired final outputs.
# Anything not listed will either not be generated or will be deleted after the pipeline is run.
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

# Snakemake rule that runs NanoPlot on the concatenated fastq files for each sample in the sampleSheet
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

# Snakemake rule that runs FastQScreen on the concatenated fastq files for each sample in the sampleSheet.
# The configuration for the queried databases is specified in the /src/fastq_screen.conf file.
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

# Rule that aligns the concatenated fastq files to the specified genome(s) using minimap2 and then
# sorts the resulting sam file with samtools into a bam file.
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

# Snakemake rule that makes a read depth graph for each sample in the sampleSheet for each DEFAULTGENOME specifiec.
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
# Snakemake rule that filters the concat FASTQ for each sample for reads that contain sequences in the query.fa file and
# aligns the filtered reads to the specified genome(s) using minimap2. The resulting sam file is then sorted with samtools into a bam file.
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
        bash ./src/query_align.sh {input.fastq} {input.fasta} {output.sam} {output.bam} {params.query} {output.filteredFasta}
        """

# Snakemake rule that runs medaka_consensus on the filtered FASTA files for each sample in the sampleSheet.
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
        bash ./src/make_consensus.sh {input.fasta} {params.genome} {output.bam} {output.index}
        """
