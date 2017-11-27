Documentation for all functions available from the bcbioSingleCell package is available at [http://bioinformatics.sph.harvard.edu/bcbioSingleCell/reference/index.html](http://bioinformatics.sph.harvard.edu/bcbioSingleCell/reference/index.html)

# bcbioSingleCell QC Report

1. Use the information from the client to construct the metadata table to use with bcbioSingleCell R package according to the specifications detailed at [https://github.com/hbc/bcbioSingleCell](https://github.com/hbc/bcbioSingleCell).
	- **Important:** the `Sequence` column for the inDrop metadata is the **Forward** sequence, not the same as the sequences present in the `sample_barcodes` file, which is the reverse complement.

### Quality control report

#### Setting up

2. Choose the quality control template.

3. Edit the information in the files `_header.Rmd` and `_footer.Rmd` with experiment-specific information.

4. Install `bcbioSingleCell` and load the library:
	
	```r
	# devtools::install_github("hbc/bcbioSingleCell" # Add argument `ref = "develop"` if need development branch
	
	library(bcbioSingleCell)
	```
	
5. Bring in data from bcbio:
	
	```r
	bcbio <- loadSingleCell("~/bcbio/PIs/path/to/final/",
                        interestingGroups = "sampleName",
                        sampleMetadataFile = "~/path/to/metadata", 
                        gtfFile = "~/bcbio/PIs/path/to/Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf")
	
	save(bcbio_output, file="data/bcb.rda")
	```
	
	> **NOTE:** Reading in the GTF file can take a long time.

6. Follow template - run entire `r setup` chunk by clicking on the green triangle at the top of the setup chunk (if you clear your environment, you need to run the chunk this way to make the `params` reappear.

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

7. For the count alignment, be sure to update the **linked Ensembl** to be accurate for the organism. This information is present in the file: `_footer.Rmd`. 

8. To explore the raw data stored inside the `bcb` object, the following functions can be helpful:
	
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

9. Evaluate the number of reads per cell:

	```r
	plotReadsPerCell(bcb, filterCells = FALSE)
	```
	
	The three plots give different ways of looking at the number of reads per cell. Generally you would like to see a large peak at around 10,000 reads per cell, and you hope your filtering threshold of 1,000 reads per cell used in bcbio has removed the poor quality cells with few number of reads. The filtering threshold of 1,000 is represented by the vertical dotted line.

	For example, in the figures below, the yellow sample is worrisome because we see a small peak at 10,000 reads per cell, but a much larger peak at 1,000 reads per cell. The larger peak merges into the poor quality cells with few reads per cell.
	
	<img src="../img/sc_qc_reads_ridgeline.png" width="600">
	
	The proportional histogram looks a bit better, as you hope to see all of the samples with peaks in relatively the same location between 10,000 and 100,000 reads per cell. However, the yellow sample still has this shoulder, which is indicative of many poor quality cells. If this were the only issue with the data, we may want to set the threshold to be more strict to ~10,000 reads per cell to get rid of the cells constituting the shoulder in the yellow sample.

	<img src="../img/sc_qc_reads_histogram.png" width="500">
	
##### Cell counts

10. Determine the number of cells detected per sample:

	```r
	plotCellCounts(bcb, filterCells = FALSE)
	```

	The cell counts are determined by the number of unique cellular barcodes detected. During the inDrop protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. Upon treatment of UV and cell lysis, all components mix together inside the droplet and reverse transcription proceeds, followed by droplet breakup and linear amplification for library preparation. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. We often see all possible combinations of cellular barcodes at a low level, leading to a higher number of cellular barcodes than cells.

	You expect the number of unique cellular barcodes to be around the number of sequenced cells (determined in step 1) or greater due to some hydrogels having more than one cellular barcode. The yellow sample below seems to have at least double the number of cellular barcodes as the other samples.

	<img src="../img/sc_qc_cellcounts.png" width="500">

##### UMI counts per cell

11. Determine the number of UMI counts (transcripts) per cell:

	```r
	plotUMIsPerCell(
    		bcb,
    		filterCells = FALSE,
	    	min = params$minUMIs)
	```

	<img src="../img/sc_qc_umisPerCell.png" width="500">
	
##### Genes detected per cell

12. Discover the number of genes detected per cell:

	```r
	plotGenesPerCell(
	    bcb,
	    filterCells = FALSE,
	    min = params$minGenes,
	    max = params$maxGenes)
	```

	<img src="../img/sc_qc_genesDetected.png" width="500">
	
##### UMIs vs. genes detected

13. Identify whether large number of poor quality cells present in any samples with low UMI/genes detected:

	```r
	plotUMIsVsGenes(bcb, filterCells = FALSE)
	```

	<img src="../img/sc_qc_UMIsVsGenesDetected.png" width="500">
	
##### Mitochondrial counts ratio

14. Identify whether there is a large amount of mitochondrial contamination from dead or dying cells:

	```r
	plotMitoRatio(
	    bcb,
	    filterCells = FALSE,
	    max = params$maxMitoRatio)
	```

	<img src="../img/sc_qc_mitoRatio.png" width="500">
	
##### Novelty

15. Explore the novelty for contamination with low complexity cell types:

	```r
	plotNovelty(
	    bcb,
	    filterCells = FALSE,
	    min = params$minNovelty)
	```
	
	<img src="../img/sc_qc_novelty.png" width="500">
	

##### Filtered results

16. Run the filtering criteria and explore the plots again. The metrics should have improved greatly after removing low gene/UMI cells and high mitochondrial cells.

	```r
	bcbFiltered <- filterCells(bcb,
	minUMIs = params$minUMIs,
	minGenes = params$minGenes,
	maxGenes = params$maxGenes,
	maxMitoRatio = params$maxMitoRatio,
	minNovelty = params$minNovelty,
	minCellsPerGene = params$minCellsPerGene)
	```

	One main plot to look at to determine the success of the filtering criteria is the number of cell counts. You should expect roughly the number of sequenced cells per sample. We found out from the client that they had sequenced 2000-3000 cells, so the final numbers were around our expectations. If the number of cells sequenced is vastly different than the number returned after filtering, then you may need to re-visit the threshold criteria used for filtering.
	
	**Cell counts**
	
	<img src="../img/sc_qc_filtered_cellcounts.png" width="500">
	
	In addition, it is a good idea to explore all of the quality plots for the filtered data. All plots should be much improved for the number of reads per cell, genes detected, UMIs per cell, mitochondrial ratio, and novelty. The plots below show the filtered plots from the example data. Since the `Unsorted` sample was a poor quality sample, the filter will remove a large number of the cells for this sample. 
	
	**Reads per cell**
	
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
	
17. When you are satisfied with the filtered results, save the filtered data.

	```r
	assignAndSaveData(name = "bcbFiltered", object = bcbFiltered, dir = dataDir)
	```
