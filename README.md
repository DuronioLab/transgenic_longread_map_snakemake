# Cadillac locus long read sequencing
## Author: Markus Nevil

## Quick Start:

Clone pipeline
```
git clone https://github.com/markus-nevil/transgenic_longread_map_snakemake.git && mv transgenic_longread_map_snakemake/* . && rm -rf transgenic_longread_map_snakemake

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

**Note:** Each row indicates a different experiment (i.e. "barcode11" is a folder holding FASTQ files from a single run).
Each row can have different values for each column, or the same values. These determine how the files in the `sampleDirectory` are processed.



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


