
```{r}
## Install recount from Bioconductor
install.packages("BiocManager")
BiocManager::install('recount')

## Browse the vignetets for a quick description of how to use the package
library('recount')
browseVignettes('recount')

## Download the RangedSummarizedExperiment object at the gene level for 
## study SRP009615
url <- download_study('SRP009615')

## View the url for the file by printing the object url
url

## Load the data
load(file.path('SRP009615', 'rse_gene.Rdata'))

## Scale counts
rse <- scale_counts(rse_gene)

## Then use your favorite differential expression software

## For more details, check the recount package vignette at
## http://bioconductor.org/packages/recount
```
