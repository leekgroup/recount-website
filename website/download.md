
```{r}
## If needed
install.packages('devtools')

## Install recount from GitHub
devtools::install_github('leekgroup/recount')
## Will be available from Bioconductor in the future
# source('http://bioconductor.org/biocLite.R')
# biocLite('recount')

## Browse the vignetets for a quick description of how to use the package
library('recount')
browseVignettes('recount')

## Download the RangedSummarizedExperiment object at the gene level for 
## study SRP009615
url <- download_study('SRP009615')

## This is the original url for the file
url

## Load the data
load(file.path('SRP009615', 'rse_gene.Rdata'))

## Scale counts
rse <- scale_counts(rse_gene)

## Then use your favorite differential expression software
```
