# Seurat clustering of cells

We will perform clustering using the output from our QC analysis. To use the filtered data, fill in the `params` for `bcbFile` with the path to the filtered output data.

```r
title: "Seurat Clustering"
author: "`r getOption('author')`"
date: "`r Sys.Date()`"
bibliography: bibliography.bib
params:
    bcbFile: "data/bcbFiltered.rda"
    seuratName: "seurat"
    pcCompute: 20
    pcUse: FALSE
    varsToRegress: !r c("nUMI")
    resolution: 0.8
    outputDir: "."
---
```

```{r setup, cache=FALSE, message=FALSE, warning=FALSE}
library(bcbioSingleCell)

# Shared RMarkdown settings
prepareSingleCellTemplate()
if (file.exists("setup.R")) {
    source("setup.R")
}

# Directory paths
dataDir <- file.path(params$outputDir, "data")
markersDir <- file.path(params$outputDir, "results", "markers")
dir.create(markersDir, recursive = TRUE, showWarnings = FALSE)

# Load bcbioSingleCell object
bcbName <- load(params$bcbFile)
bcb <- get(bcbName, inherits = FALSE)

# Seurat object names and file paths
seuratFile <- file.path(dataDir, paste0(params$seuratName, ".rda"))
seuratPreregressName <- paste(params$seuratName, "preregress", sep = "_")
seuratPreregressFile <- file.path(dataDir, paste0(seuratPreregressName, ".rda"))

# Check to see if seurat objects exist
if (exists(c(params$seuratName, seuratPreregressName), inherits = FALSE)) {
    evalSeurat <- FALSE
} else {
    if (all(file.exists(seuratFile, seuratPreregressFile))) {
        loadDataAsName(
            c(seurat = seuratFile,
              seuratPreregress = seuratPreregressFile))
        evalSeurat <- FALSE
    } else {
        evalSeurat <- TRUE
    }
}

# Vector to use for plot looping
if (length(unique(sampleMetadata(bcb)$sampleName)) == 1) {
    groupBy <- "ident"
} else {
    groupBy <- c("ident", "sampleName", "phase", interestingGroups(bcb))
}


# knitr arguments (for `rmarkdown::render()` looping)
opts_chunk$set(
    cache.path = paste(
        params$seuratName,
        "clustering",
        "cache/",
        sep = "_"),
    fig.path = paste(
        params$seuratName,
        "clustering",
        "files/",
        sep = "_")
)
```

```{r header, child="_header.Rmd", eval=file.exists("_header.Rmd")}
```



```{r sample_metadata}
sampleMetadata(bcb)
```

This workflow is adapted from the following sources:

- Satija Lab: [Seurat v2 Guided Clustering Tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html)
- Paul Hoffman: [Cell-Cycle Scoring and Regression](http://satijalab.org/seurat/cell_cycle_vignette.html)



* * *



# Initialize [Seurat][] (`r bcbName`)

Prior to any clustering analysis, the raw counts need to be normalized using global-scaling normalization. Global-scaling normalization (1) normalizes the gene expression measurements for each cell by the total expression, (2) multiplies this by a scale factor (10,000 by default), and (3) log-transforms the result. Following normalization, the average expression and dispersion for each gene is calculated, which places these genes into bins, and then a z-score for dispersion within each bin is calculated. This helps control for the relationship between variability and average expression. Finally, the genes are scaled and centered.
