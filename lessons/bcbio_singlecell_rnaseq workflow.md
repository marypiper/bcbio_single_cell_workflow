# bcbio_singlecell_rnaseq workflow

## Setting up for bcbio single cell RNA-Seq analysis

1. Ask client for the following:
	- How many samples were sequenced?
	- What were the sample indices used?
	- How many cells were encapsulated and sequenced per sample?
	- What is the main experimental question - does it require clustering using markers and/or cell trajectory analyses?

2. Acquire data from sequencing core:

	- **Bauer sequencing core:** uses Basespace and to download the sequencing files, follow the code below:

		```
		wget https://da1s119xsxmu0.cloudfront.net/sites/knowledgebase/API/08052014/Script/BaseSpaceRunDownloader_v2.zip
		unzip BaseSpaceRunDownloader_v2.zip
		python BaseSpaceRunDownloader_v2.py -r <Run ID> -a <access token>
		```
		
		The option `-r` is the number in the basespace url and the [access token](https://developer.basespace.illumina.com/docs/content/documentation/authentication/obtaining-access-tokens) is something you have to get for your basespace account. 
		
		The files output will be BCL files that require demultiplexing with the `bcl2fastq` tool (instructions below).

	- **DFCI sequencing center (Zach):** will output the FASTQ files (already demultiplexed).
		
	- **Biopolymers sequencing facility:** will sometimes output BCL and sometimes FASTQ, so necessary to check the files - good idea to ask for the BCL files
		
	- **Broad Institute:** has their own single cell distribution platform

3. If downloaded sequencing files are BCL format, then need to convert to FASTQ. To do this log on to Orchestra or O2 to run `bcl2fastq`. Currently the module is only available on Orchestra.

	- Change directories to the sequencing folder downloaded from the facility. The folder should be arranged according to the image below for NextSeq or MiniSeq:
	

		<img src="../img/sequencing_dir_org.png" width="400">
	
		*Image acquired from [bcl2fastq documentation](../docs/bcl2fastq2_guide_15051736_v2.pdf).* 
	
	- Load `bcl2fastq` module and convert files to FASTQ by using the following command:
		
		```
		bcl2fastq \
		--use-bases-mask y*,y*,y*,y* \
		--mask-short-adapter-reads 0 \
		--minimum-trimmed-read-length 0
		```
		
		More information regarding the `bcl2fastq` command and directory structures for other sequencing machines can be found in the [documentation](../docs/bcl2fastq2_guide_15051736_v2.pdf). 
		
4. The output files should be in the `BaseCalls` directory. For each file of sequenced reads, there should be four associated FASTQ files (R1-R4). These 
	
	- **R1 (61 bp Read 1):** sequence of the read
	- **R2 (8 bp Index Read 1 (i7)):** cellular barcode - which cell read originated from
	- **R3 (8 bp Index Read 2 (i5)):** library index - which sample read originated from
	- **R4 (14 bp Read 2):** read 2 and barcode/UMI - remaining cellular barcode and UMI - which transcript read originated from (to find PCR duplicates)

	The reads for each sequence are depicted in the image below:

	<img src="../img/sc_seq_method.png" width="800">
	
	*Image credit: Sarah Boswell, Harvard Staff Scientist for Sequencing Technologies*

	
5. To quickly view the counts for the barcodes with the top five highest counts based on the first 10,000 reads in a file:

	```
	gzip -cd filename_R3.fq.gz | head -40000 | awk 'NR % 4 == 2' | sort | uniq -c | awk 	'{ print $2 "," $1}' | sort -t"," -n --key=2 | tail -5
	```
	
	*NOTE: `awk 'NR % 4 == 2'` gets every 4th line starting from the 2nd, which is a useful trick when you want to count up FASTQ file entries (Rory's code)*

	The reverse complement sequences of the sample indices given by the client should correspond to the most abundant indices in the file.
	
6. Use the `cat` command to concatenate all of the files for a given sample across lanes:

	```
	cat Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L002_R1_001.fastq.gz Undetermined_S0_L003_R1_001.fastq.gz Undetermined_S0_L004_R1_001.fastq.gz > cat_R1.fastq.gz
	```
	Do the same for the R2, R3, and R4 files.

7. Create metadata file as normal for bcbio run.

	```
	fileName,description
	cat,run1
	```


8. Download the most recent transcriptome FASTA and GTF (patched scaffold) files:

	```
	# Most recent mouse FASTA from Ensembl FTP
	wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
	
	# Most recent mouse GTF from Ensembl FTP
	wget ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf.gz
	
	# Perform the checksums
	sum Mus_musculus.GRCm38.cdna.all.fa.gz
	sum Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf.gz
	
	# Decompress FASTA and GTF to run in bcbio
	gzip -d Mus_musculus.GRCm38.cdna.all.fa.gz
	gzip -d Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf.gz
	```
## Overview of bcbio single cell RNA-Seq workflow on O2

The bcbio single cell RNA-Seq pipeline will perform the following steps:

1. Identify the sample barcodes in the R3 read, which were provided in the `config` file with the `sample_barcodes` parameter. A single mismatch between known sample barcodes and sequences is allowed.

2. Identify the cellular barcodes by parsing the R2 and R4 reads. The cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. Upon treatment of UV and cell lysis, all components mix together inside the droplet and reverse transcription proceeds, followed by droplet breakup and linear amplification for library preparation. **While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. We often see all possible combinations of cellular barcodes at a low level, leading to a higher number of cellular barcodes than cells.**

3. Identify the unique molecular identifiers (UMIs) by parsing R4 read.

4. Filter out the sequence data with cellular barcodes matching less than 1000 reads (indicating poor quality cells due to encapsulation of free floating RNA from dying cells, small cells, or set of cells that failed for some reason). The threshold for the number of matching reads used for filtering can be specified in the `config` file with the `minimum_barcode_depth` parameter.

5. Align reads with Rapmap tool.

6. Take reads that mapped to more than one transcript and divide the count between all of the transcripts to which the reads aligned.

## Running bcbio single cell RNA-Seq workflow on O2

1. Create configuration template for single cell run:

```
details:
  - analysis: scRNA-seq
    algorithm:
      transcriptome_fasta: /n/data1/cores/bcbio/PIs/PI_name/meta/Mus_musculus.GRCm38.cdna.all.fa
      transcriptome_gtf: /n/data1/cores/bcbio/PIs/PI_name/meta/Mus_musculus.GRCm38.90.gtf
      umi_type: harvard-indrop-v3
      minimum_barcode_depth: 1000
      cellular_barcode_correction: 1
      sample_barcodes: /n/data1/cores/bcbio/PIs/PI_name/meta/hbc02055-sample-barcodes-rc.txt
    genome_build: mm10
```

2. Normal bcbio configuration file creation:

	```
	bcbio_nextgen.py -w template ../config/scRNAseq_config_template.yaml ../meta/PI_name.csv ../hbcXXXXX/seq_dir/Data/Intensities/BaseCalls/cat*fastq.gz
	```

3. Create script (below) to run job on O2 and run with `sbatch ../../runJob-PI_name-scRNAseq.slurm`:

```
#!/bin/sh
#SBATCH -p medium
#SBATCH -J win-full
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 4-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=piper@hsph.harvard.edu

export PATH=/n/app/bcbio/tools/bin:$PATH

/n/app/bcbio/dev/anaconda/bin/bcbio_nextgen.py ../config/PI_name.yaml -n 48 -t ipython -s slurm -q medium -r t=4-00:00
```

## Analyzing bcbio single cell RNA-Seq output using the bcbioSingleCell R package

12. Use the information from the client to construct the metadata table to use with bcbioSingleCell R package according to the specifications detailed at [https://github.com/hbc/bcbioSingleCell](https://github.com/hbc/bcbioSingleCell).
	- Note that the `Sequence` column for the inDrop metadata is the **Forward** sequence, not the same as the sequences present in the `sample_barcodes` file, which is the reverse complement.

### Quality control report - 
13. Choose the quality control template.

14. Bring in data from bcbio:
	
	```r
	loadSingleCell("~/bcbio/PIs/path/to/final/", 
               interestingGroups = "sampleName", 
               sampleMetadataFile = "path/to/metadata", 
               gtfFile = "~/bcbio/PIs/path/to/Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf")
	```
