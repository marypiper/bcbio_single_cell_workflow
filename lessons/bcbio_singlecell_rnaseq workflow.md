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
	
	>**NOTE:** `awk 'NR % 4 == 2'` gets every 4th line starting from the 2nd, which is a useful trick when you want to count up FASTQ file entries (Rory's code)

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

2. Identify the cellular barcodes by parsing the R2 and R4 reads. 

3. Identify the unique molecular identifiers (UMIs) by parsing R4 read.

4. Filter out the sequence data with cellular barcodes matching less than 1000 reads (indicating poor quality cells due to encapsulation of free floating RNA from dying cells, small cells, or set of cells that failed for some reason). The threshold for the number of matching reads used for filtering can be specified in the `config` file with the `minimum_barcode_depth` parameter.

5. Align reads with [Rapmap](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985/RapMap-a-rapid-sensitive-and-accurate-tool-for) tool.

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
	- **Important:** the `Sequence` column for the inDrop metadata is the **Forward** sequence, not the same as the sequences present in the `sample_barcodes` file, which is the reverse complement.

### Quality control report

#### Setting up

13. Choose the quality control template.

14. Edit the information in the files `_header.Rmd` and `_footer.Rmd` with experiment-specific information.

15. Install `bcbioSingleCell` and load the library:
	
	```r
	# devtools::install_github("hbc/bcbioSingleCell", ref = "develop", dep = FALSE)
	
	library(bcbioSingleCell)
	```
	
16. Bring in data from bcbio:
	
	```r
	bcbio <- loadSingleCell("~/bcbio/PIs/path/to/final/",
                        interestingGroups = "sampleName",
                        sampleMetadataFile = "~/path/to/metadata", 
                        gtfFile = "~/bcbio/PIs/path/to/Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf")
	
	save(bcbio_output, file="data/bcb.rda")
	```
	
	> **NOTE:** Reading in the GTF file can take a long time.

17. Follow template - run entire `r setup` chunk by clicking on the green triangle at the top of the setup chunk (if you clear your environment, you need to run the chunk this way to make the `params` reappear.

	```r
	# Shared RMarkdown settings
	prepareSingleCellTemplate()
	if (file.exists("setup.R")) {
	    source("setup.R")
	}

	# Directory paths
	dataDir <- file.path(params$outputDir, "data")

	# Load bcbioSingleCell object
	bcbName <- load(params$bcbFile)
	bcb <- get(bcbName, inherits = FALSE)
	```
	
	```r
	eval=file.exists("_header.Rmd")
	```

	```r
	sampleMetadata(bcb)
	```

18. For the count alignment, be sure to update the **linked Ensembl** to be accurate for the organism. This information is present in the file: `_footer.Rmd`. 

19. To explore the raw data stored inside the `bcb` object, the following functions can be helpful:
	
	```r
	# Access metadata for each sample: "sampleID", "sampleName", "description", "fileName", "index", "sequence", "revcomp"
	sampleMetadata(bcb)
	
	# Access metadata for each cell: "nCount", "nUMI", "nGene", "nCoding", "nMito", "log10GenesPerUMI", "mitoRatio" 
	colData(bcb)
	
	# Access raw counts - each column represents a single cell
	counts <- as.data.frame(as.matrix((assay(bcb))))
	
	# Can return cells from a particular sample by using metadata information about which sample corresponds to each barcode
	unsort_counts <- counts[, str_detect(colnames(counts), "run1_ATTAGACG")] # Return only the counts for the `Unsorted` sample
	
	# Extract information associated with each gene including "ensgene", "symbol", "description", "biotype", "broadClass"
	rowData(bcb)
	
	# Return the genes that are used to determine mitochondrial contamination
	subset(rowData(bcb), broadClass == "mito")
	```

#### Quality Control Metrics

##### Reads per cell

20. Evaluate the number of reads per cell:

	```r
	plotReadsPerCell(bcb, filterCells = FALSE)
	```
	
	The three plots give different ways of looking at the number of reads per cell. Generally you would like to see a large peak at around 10,000 reads per cell, and you hope your filtering threshold of 1,000 reads per cell used in bcbio has removed the poor quality cells with few number of reads. The filtering threshold of 1,000 is represented by the vertical dotted line.

	For example, in the figures below, the yellow sample is worrisome because we see a small peak at 10,000 reads per cell, but a much larger peak at 1,000 reads per cell. The larger peak merges into the poor quality cells with few reads per cell.
	
	<img src="../img/sc_qc_reads_ridgeline.png" width="600">
	
	The proportional histogram looks a bit better, as you hope to see all of the samples with peaks in relatively the same location between 10,000 and 100,000 reads per cell. However, the yellow sample still has this shoulder, which is indicative of many poor quality cells. If this were the only issue with the data, we may want to set the threshold to be more strict to ~10,000 reads per cell to get rid of the cells constituting the shoulder in the yellow sample.

	<img src="../img/sc_qc_cellcounts.png" width="500">
	
##### Cell counts

21. Determine the number of cells detected per sample:

	```r
	plotCellCounts(bcb, filterCells = FALSE)
	```

	The cell counts are determined by the number of unique cellular barcodes detected. During the inDrop protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. Upon treatment of UV and cell lysis, all components mix together inside the droplet and reverse transcription proceeds, followed by droplet breakup and linear amplification for library preparation. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. We often see all possible combinations of cellular barcodes at a low level, leading to a higher number of cellular barcodes than cells.

	You expect the number of unique cellular barcodes to be around the number of sequenced cells (determined in step 1) or greater due to some hydrogels having more than one cellular barcode. The yellow sample below seems to have at least double the number of cellular barcodes as the other samples.

	<img src="../img/sc_qc_reads_histogram.png" width="500">

##### UMI counts per cell

22. Determine the number of UMI counts (transcripts) per cell:

	```r
	plotUMIsPerCell(
    		bcb,
    		filterCells = FALSE,
	    	min = params$minUMIs)
	```

	<img src="../img/sc_qc_umisPerCell.png" width="500">
	
##### Genes detected per cell

23. Discover the number of genes detected per cell:

	```r
	plotGenesPerCell(
	    bcb,
	    filterCells = FALSE,
	    min = params$minGenes,
	    max = params$maxGenes)
	```

	<img src="../img/sc_qc_genesDetected.png" width="500">
	
##### UMIs vs. genes detected

24. Identify whether large number of poor quality cells present in any samples with low UMI/genes detected:

	```r
	plotUMIsVsGenes(bcb, filterCells = FALSE)
	```

	<img src="../img/sc_qc_UMIsVsGenesDetected.png" width="500">
	
##### Mitochondrial counts ratio

25. Identify whether there is a large amount of mitochondrial contamination from dead or dying cells:

	```r
	plotMitoRatio(
	    bcb,
	    filterCells = FALSE,
	    max = params$maxMitoRatio)
	```

	<img src="../img/sc_qc_mitoRatio.png" width="500">
	
##### Novelty

26. Explore the novelty for contamination with low complexity cell types:

	```r
	plotNovelty(
	    bcb,
	    filterCells = FALSE,
	    min = params$minNovelty)
	```
	
	<img src="../img/sc_qc_novelty.png" width="500">
	

##### Filtered results

27. One main plot to look at to determine the success of the filtering criteria is the number of cell counts. You should expect roughly the number of sequenced cells per sample. We found out from the client that they had sequenced 2000-3000 cells, so the final numbers were around our expectations. If the number of cells sequenced is vastly different than the number returned after filtering, then you may need to re-visit the threshold criteria used for filtering.
	
	**Cell counts**
	
	<img src="../img/sc_qc_filtered_cellcounts.png" width="500">
	
	In addition, it is a good idea to explore all of the quality plots for the filtered data. All plots should be much improved for the number of reads per cell, genes detected, UMIs per cell, mitochondrial ratio, and novelty. The plots below show the filtered plots from the example data. Since the `Unsorted` sample was a poor quality sample, the filter will remove nearly all of the cells for this sample. 
	
	** Reads per cell**
	
	<img src="../img/sc_qc_filtered_reads.png" width="500">
	
	**Genes detected**
	
	<img src="../img/sc_qc_filtered_genesDetected.png" width="500">
	
	**UMIs per cell**
	
	<img src="../img/sc_qc_filtered_umisPerCell.png" width="500">
	
	**UMIs versus genes detected**
	
	<img src="../img/sc_qc_filtered_UMIsVsGenesDetected.png" width="500">
	
	**Mitochondrial ratio**
	
	<img src="../img/sc_qc_filtered_mitoRatio.png" width="500">
	
	**Novelty**
	
	<img src="../img/sc_qc_filtered_novelty.png" width="500">
	

	
