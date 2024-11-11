
# cellXY      [![DOI](https://zenodo.org/badge/535519225.svg)](https://zenodo.org/doi/10.5281/zenodo.11329594)

<!-- badges: start -->
<!-- badges: end -->

The cellXY package currently contains trained models to classify cells as male or 
female and to predict whether a cell is a male-female doublet or not. 

The classifySex function takes a count matrix as input, computes required features and predict the sex label of each cell with the trained model. We have trained models for human and mouse cells seperately, and you need to specify the genome type of your data. 
Similarly, the findMfDoublet function uses trained machine learning models to identify 
male-female doublet cells in the dataset. 

## Installation

If you would like to view the cellXY vignette, you can install the released 
version of cellXY from [github](https://github.com/phipsonlab/cellXY) using the 
following commands:

``` r
# devtools/remotes won't install Suggested packages from Bioconductor
BiocManager::install(c("CellBench", "BiocStyle", "scater"))

remotes::install_github("phipsonlab/cellXY", build_vignettes = TRUE, 
dependencies = "Suggest")
```

In order to view the vignette for cellXY use the following command:

``` r
browseVignettes("cellXY")
```

If you don't care to view the glorious vignette you can also install cellXY as 
follows:

``` r
library(devtools)
devtools::install_github("phipsonlab/cellXY")
```

## Sex label prediction example 

This is a basic example which shows you how to obtain a sex label prediction for each cell. 

``` r
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(org.Hs.eg.db)
library(scRNAseq)
sc_data <- load_sc_data()
sc_10x <- sc_data$sc_10x

counts <- counts(sc_10x)
ann <- select(org.Hs.eg.db, keys=rownames(sc_10x),
              columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
m <- match(rownames(counts), ann$ENSEMBL)
rownames(counts) <- ann$SYMBOL[m]

sex <- classifySex(counts, genome="Hs")

table(sex$prediction)
boxplot(counts["XIST",]~sex$prediction)


# Mouse data example 
sce <- fetchDataset("zilionis-lung-2019", "2023-12-20", path="mouse")
mouse_cm <- counts(sce)
# make sure the row names are the gene symbols
row.names(mouse_cm) <- row.names(sce)
# make sure the column (cell) names are unique
colnames(mouse_cm) <- paste("cell", 1:ncol(sce), sep="_")

mouse_pred <- classifySex(mouse_cm, genome="Mm")
table(mouse_pred$prediction)
```

