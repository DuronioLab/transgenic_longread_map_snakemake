# Cadillac locus long read sequencing
## Author: Markus Nevil

## Quick Start:

Clone pipeline (Current stable version: v1.7.2)
```
git clone https://github.com/markus-nevil/transgenic_longread_map_snakemake.git --depth 1 && cd cutNrun-pipeline/ && rm -rf .git
```

Create `sampleInfo.tsv` ([see below](#sampleInfo)) with descriptive columns of data.
**Note:** Paths are relative to the location of the `sampleInfo.tsv` file.

Ensure that the FASTA file with sites to be used as anchors is named properly in the `querySites` column.

```
sampleDirectory	sampleName	rep	refGenome	querySites
barcode11	4xCad	1	dHisCcadillacP4xD4x	query.fa
barcode12	4xCad	2	dHisCcadillacP4xD4x	query.fa
barcode10	6xCad	1	dHisCcadillacP6xD6x	query.fa
```

Set desired default genome (the entire genome) and reference genome (transgenic locus)  `config.json` ([see below](#config)). Set multiple genomes by passing as an array (e.g. `["dm6", "droYak2"]`).

Include the file names for the genomes in the `genome` section of `config.json`.

```
{
  "sampleInfo" : "sampleInfo.tsv",
  "sampleDirs" : "sampleDirectory",
  "baseNameCols": ["sampleName", "rep"],
  "defaultGenome": "dm6_mt_wMel",
  "refGenome": "dHisCcadillacP4xD4x",
  "querySites": "querySites",
  "genome" : {
	"dm6_mt_wMel" : {
	  "fasta" : "dm6_mt_wMel.fasta",
	},
	"dHisCcadillacP4xD4x" : {
	  "fasta" : "dHisC_Cadillac_D4x_P4x.fa"
	}
  }
```

Edit `slurmConfig.json` to configure default parameters if necessary.

Load python on Longleaf with `module load python`

Start the submission with `sh slurmSubmission.sh`

## Required files and directory structure

After cloning the repository, and placing your own `Default Genome` and `Reference Genome` in the `config.json` file, the directory should look like:

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


