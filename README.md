# Cadillac locus long read sequencing
## Author: Markus Nevil

**Table of Contents**
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Required files](#required-files)
   - [query.fa](#query)
   - [sampleInfo.tsv](#sampleInfo)
   - [config.json](#config)
- [Optional files](#optional-files)
  - [Features](#features)
- [Required directory structure](#required-directory-structure)
- [Expected Output](#expected-output)
- [Getting your data to UNC Longleaf](#getting-your-data-to-unc-longleaf)
   - [Set up Globus on your computer](#set-up-globus-on-your-computer-only-needs-to-be-done-once)
   - [Moving data from the GridIon](#moving-data-from-the-gridion-to-a-safe-location)
   - [Moving data to Longleaf](#moving-data-from-a-safe-location-to-longleaf)
- [Acknowledgements](#acknowledgements)

## Introduction:
This pipeline is designed to take whole-genome long-read sequencing data from an Oxford Nanopore sequencer to perform QC steps and align the reads. 
The alignment is done to one-or-more "default genomes", as well as an "anchored" alignment to one-or-more "reference genomes" which are small regions of interest, likely ones with transgenic insertions.
The "anchoring" step is done to ensure that the reads aligned to the transgenic region maintains the structure of the locus. Specficially this is designed with
transgenic histone gene arrays in mind, where traditional methods of alignment fail due to the repetitive nature of the transgenic histone
gene locus. The pipeline is designed to be run on the UNC Longleaf cluster, but can be run on any system with Snakemake and the 
appropriate programs installed.

## Quick Start:

Clone pipeline
```
git clone https://github.com/markus-nevil/transgenic_longread_map_snakemake.git && mv transgenic_longread_map_snakemake/* . && rm -rf transgenic_longread_map_snakemake
```

> [!IMPORTANT]
> Review this section carefully, especially after cloning the repository. If you had files of the same name, they may have been rewritten.

1. **Required:** Add/edit the `query.fa` file with the sites to be used as anchors for the alignment ([see below](#query)).

2. **Required:** Edit `sampleInfo.tsv` to include the sample information ([see below](#sampleInfo)).

3. **Required:** Edit `config.json` to include the reference genome(s) and default genome ([see below](#config)).

4. **Optional:** Edit `slurmConfig.json` to configure default parameters if necessary.

Start the pipeline submission with:
```
sh slurmSubmission.sh
```
## Required files

### query

The query.fa (or query.fasta) is a FASTA file that contains the sites to be used as anchors for the alignment. This file should contain multiple 40-50 bp sequences, each of unique sequence flanking repeat regions.

```
>site1
ACGTACTCTGCTACTCTGACCCCATACGACAC
>GFP_promoter
CCCATAGGAAATCATATCGGATCTAGCGACAC
>my_gene
AAACATATATAATACCGCGATCGCGATAATGA
```

### sampleInfo

> [!NOTE]
> Paths are relative to the location of the `sampleInfo.tsv` file.

Ensure that the FASTA file with sites to be used as anchors is named properly in the `querySites` column.

```
sampleDirectory	sampleName	rep	refGenome	querySites
barcode11	4xCad	1	dHisC_Cadillac_D4x_P4x	query.fa
barcode12	4xCad	2	dHisC_Cadillac_D4x_P4x	query.fa
barcode10	6xCad	1	dHisC_Cadillac_D6x_P6x	query.fa
```
For more complex experiments, you may want to include multiple refGenomes and querySites to be analyzed
at the same time. For example:

```
sampleDirectory	sampleName	rep	refGenome	querySites
barcode10	4xCad	1	dHisC_Cadillac_D4x_P4x	query.fa
barcode10	4xCad	1	dHisC_Cadillac_D6x_P6x	query.fa
barcode11	6xCad	1	dHisC_Cadillac_D6x	query_oneSide.fa
barcode12	6xCad	1	dHisC_Cadillac_D6x_P6x	query.fa
barcode13	6xCad	1	dHisC_Cadillac_D12x	query_12x.fa
```
Explanation for each column:
1. **sampleDirectory:** The directory containing the FASTQ files.
   1. At the moment, the pipeline assumes all FASTQ files are .fastq.gz files.
2. **sampleName:** A short, easy, and unique (if required) identifier for the sample.
3. **rep:** The replicate number for the sample.
   1. By default replicates are combined.
4. **refGenome:** The reference genome (usually the small transgenic locus) to align the reads to.
   1. This refGenome name must match the refGenome name in the `config.json` file.
   2. To change the `Default Genome` (the entire genome) edit the `config.json` file.
5. **querySites:** The FASTA file containing the sites to be used as anchors for the alignment.
	1. This FASTA file must be multiple 40-50 bp sequences, each of unique sequence flanking repeat regions.

> [!NOTE]
> Each row indicates a different experiment (i.e. "barcode11" is a folder holding FASTQ files from a single run).
> Each row can have different values for each column, or the same values. These determine how the files in the `sampleDirectory` are processed.

### config
The `config.json` file contains the configuration for the pipeline and is found in the `src` directory.
The first half of the config file needs to be edited to match the names of your genomes. The last half should only
be edited if a different program version is required. In the following examples, only the first half is shown.

Edit the `defaultGenome` and `refGenome` fields to match the names of your entire genome(s) and reference of your transgenic locus/loci respectively.

Ensure that `genome` has an entry for each `defaultGenome` and `refGenome`, and each has an entry for the `fasta` file.

**Optional:** Add one-or-more 1x arrays of your reference genome to the `singleArray` field. More information in [the optional files section](#optional-files).
This may otherwise remain empty.

```
{
  "defaultGenome": "dm6_mt_wMel",
  "refGenome": "dHisC_Cadillac_D4x_P4x",
  "singleArray": "",
  "genome" : {
    "dm6_mt_wMel" : {
      "fasta" : "dm6_mt_wMel.fasta",
    },
    "dHisC_Cadillac_D4x_P4x" : {
      "fasta" : "dHisC_Cadillac_D4x_P4x.fa"
    }
  }
```

Sometimes, you may want to align the same data to multiple genomes, or you may wish to process multiple
samples with different genomes. In these cases, you can use a list (by using square brackets [ ] ) for `defaultGenome` and `refGenome`.

An example of a case where you want multiple `defaultGenome` and `refGenome` entries:

```
{
  "defaultGenome": ["dm6_mt_wMel", "dm6"],
  "refGenome": ["dHisC_Cadillac_D4x_P4x", "dHisC_Cadillac_D6x_P6x"],
  "singleArray": "",
  "genome" : {
    "dm6_mt_wMel" : {
      "fasta" : "dm6_mt_wMel.fasta",
    },
    "dm6" : {
      "fasta" : "dm6.fasta",
    },
    "dHisC_Cadillac_D4x_P4x" : {
      "fasta" : "dHisC_Cadillac_D4x_P4x.fa"
    },
    "dHisC_Cadillac_D6x_P6x" : {
      "fasta" : "dHisC_Cadillac_D6x_P6x.fa"
    }
  }
}
```
## Optional files
### Features
The pipeline is able to generate a GFF annotation file for easier viewing of aligned data in genome viewers (like IGV).
The GFF files are generated for each `Reference Genome` only if the a `Features from X.txt` file from SnapGene is provided.
This file should be generated from the same FASTA file used to generate the `Reference Genome` FASTA file and `X` should be
the name of the reference genome.

Follow these steps to generate the `Features from X.txt` file:
1. Open your annotated genome file in Snapgene.
2. Go to `Features` -> `Export Feature Data`
3. Select all 5 feature data types and export the file as a `.txt` file.
   1. Ensure the name of the file is `Features from X.txt` where `X` is the name of the reference genome.
4. Upload this file to your `Project_Folder` directory.

If the pipeline does not find the `Features from X.txt` file or the pipeline does not find a file that matches an 
expected `Reference Genome`, it will not generate the GFF file.

The resulting GFF file will be named *[ReferenceGenome]*.gff and will be moved to the Alignment directory. Open this file in IGV
to view the features of the reference genome.

### Single Array
Determining whether any single array repeat is incorrect is difficult, but can be partially overcome by examining all the
arrays over a single 1x reference sequence. By including one-or-more 1x array references in the directory and listing them
in the `config.json` file, the pipeline will automatically align the reads to this array. This is done by finding 1x arrays
in each read and creating "subreads" of each 1x array in the read. Thus, a read spanning 5 repeats should generate 5, 1x
repeat subreads.

Requirements:
1. A 1x reference genome FASTA file added to the `Project_Folder` directory.
   1. More than one file can be added. Only those listed in `config.json` will be used.
2. The name of the 1x reference genome in the `config.json` file under the `singleArray` field.
3. The FASTA file name must be added in the `config.json` file under the `genome` field.

Example `config.json` file with a single 1x array:
```
{
  "defaultGenome": "dm6_mt_wMel",
  "refGenome": "dHisC_Cadillac_D4x_P4x",
  "singleArray": "wild-type_array",
  "genome" : {
    "dm6_mt_wMel" : {
      "fasta" : "dm6_mt_wMel.fasta",
    },
    "dHisC_Cadillac_D4x_P4x" : {
      "fasta" : "dHisC_Cadillac_D4x_P4x.fa"
    },
    "wild-type_array" : {
      "fasta" : "wild-type_array.fasta"
    }
  }
}
```

## Required directory structure

After cloning the repository, and placing your own `Default Genome` and `Reference Genome` in the `Project_Folder` directory, you should
expect this structure:

```
Project_Folder/
├── src/
│   ├── chrom_graph.sh
│   ├── config.json
│   ├── fastq_screen.conf
│   ├── fqscreen.sh
│   ├── preProcessSampleConfig.py
│   ├── query_align.sh
│   ├── Snakefile.smk
│   └── slurmConfig.json
├── dm6_mt_wMel.fasta (example "Default Genome")
├── dHisC_Cadillac_P6x_D6x.fasta (example "Reference Genome")
├── query.fasta
├── sampleInfo.tsv
├── slurmSubmission.sh
└── README.md
```
You may copy or move your raw data folders into the `Project_Folder` directory.

## Expected Output

1. **Alignment files:** The aligned reads in BAM and SAM format.
   1. Alignment/*[SampleName]*_*[DefaultGenome(s)]*.sam
   2. Alignment/*[SampleName]*_*[DefaultGenome(s)]*.bam
   3. Alignment/*[SampleName]*_*[DefaultGenome(s)]*.bam.bai
   4. Alignment/*[SampleName]*_*[ReferenceGenome(s)]*.sam
   5. Alignment/*[SampleName]*_*[ReferenceGenome(s)]*.bam
   6. Alignment/*[SampleName]*_*[ReferenceGenome(s)]*.bam.bai
   7. Alignment/*[SampleName]*_queries.bed
   
   **Optional Alignment Outputs**
   1. Alignment/*[ReferenceGenome]*.gff
   2. Alignment/*[SampleName]*_*[singleArray]*.sam
   3. Alignment/*[SampleName]*_*[singleArray]*.bam
   4. Alignment/*[SampleName]*_*[singleArray]*.bam.bai
   

2. **Stats files:** Stats about the sequencing and alignments.
   1. Stats/*[SampleName]*\_concat\_*[DefaultGenome(s)]* readDepth.pdf
   2. Stats/NanoPlot-report_*[SampleName]*_concat.html
   3. Stats/*[SampleName]*\_concat\_screen.html


3. **Consensus sequence:** The consensus sequence of the aligned reads to the Reference Genome(s).
   1. Consensus/*[SampleName]*_*[ReferenceGenome(s)]*.fasta

## Getting your data to UNC Longleaf

### Set up Globus on your computer (only needs to be done once!).

Via [this guide](https://docs.globus.org/globus-connect-personal/install/), download [Globus Personal Connect](https://app.globus.org/collections/gcp) onto your computer.

After installation and configuration with your UNC account, follow the appropriate instructions ([Windows](https://docs.globus.org/globus-connect-personal/install/windows/#configuration), [Mac](https://docs.globus.org/globus-connect-personal/install/mac/#configuration))
to set up a shared folder on your computer. This is required to move files to-and-from different "endpoints".

Set up any personal shared folders you may want (your home directory is default), but also connect to Pierre and make your Pierre folder shared as well.

This can also include a removable usb drive, but the drive must be connected for data to transfer.

### Moving data from the GridIon to a safe location.

Connect to Pierre or your removable usb drive on your computer. Ensure Globus Personal Connect is running (it may not show up as it silently works in the background).

Go to the GridIon machine and navigate to the Globus online interface. Log in with your UNC account. Use the web interface to upload the read folders to Pierre or removable drive using the Globus `File Manager`.

You can verify that the transfer has completed on the web interface (don't close until completed!) and on the location you've moved the files.

### Moving data from a safe location to Longleaf.

Connect to Pierre or your removable usb drive on your computer. Ensure Globus Personal Connect is running (it may not show up as it silently works in the background).

Navigate to the Globus online interface and log in with your UNC account.

Using the web interface `File Manager`, connect to the `UNC, Research Computing DataMover` collection endpoint. On the second pane, connect to your Pierre or drive endpoint.

Using the `path` field, navigate to the folders you wish to move data to/from. For Longleaf, usually you will move data to `/work/users/u/s/username/Project_Folder` where `Project_Folder` is a new folder where you will run the pipeline.

Select the folders/files that you want to transfer to Longleaf and click the `Start` button.


## Acknowledgements
The pipeline is based on the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system and the Snakemake pipelines written by Spencer Nystrom and Chris Uyehara.


