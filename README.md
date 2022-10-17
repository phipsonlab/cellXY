
# cellXY

<!-- badges: start -->
<!-- badges: end -->

The cellXY package currently contains trained models to classify cells as male or 
female and to predict whether a cell is a male-female doublet or not. 

The propeller function performs statistical tests for differences in cell
type composition in single cell data. In order to test for differences in cell
type proportions between multiple experimental conditions at least one of the 
groups must have some form of biological replication (i.e. at least two 
samples). For a two group scenario, the absolute minimum sample size is thus 
three. Since there are many technical aspects which can affect cell type 
proportion estimates, having biological replication is essential for a 
meaningful analysis.

The propeller function takes a SingleCellExperiment or Seurat object as input,
extracts the relevant cell information, and tests whether the cell type 
proportions are statistically significantly different between experimental
conditions/groups. The user can also explicitly pass the cluster, sample and 
experimental group information to propeller. P-values and false discovery rates 
are outputted for each cell type. 

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

This is a basic example which shows you how to obtain a sex label preidction for each cell. 

``` r
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(org.Hs.eg.db)

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
```
Please note that this basic implementation is for when you are only modelling
group information. When you have additional covariates that you would like to 
account for, please use the propeller.ttest() and propeller.anova() functions
directly. Please read the vignette for examples on how to model a continuous 
variable, account for additional covariates and include a random effect in the 
analysis. 


