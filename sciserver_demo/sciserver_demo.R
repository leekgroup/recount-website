## This code is part of the recount vignette section describing how to use
## SciServer compute http://www.sciserver.org/ to access the recount data.

# Install recount
install.packages("BiocManager")
BiocManager::install(c('recount', 'DESeq2'))

# Load package
library('recount')

## Find a project of interest
project_info <- abstract_search('GSE32465')

## Path to the recount files in SciServer
scipath <- '/home/idies/workspace/recount01'

## Load the data locally on SciServer
load(file.path(scipath, project_info$project, 'rse_gene.Rdata'))

## Data in R
rse_gene

## View GEO ids
colData(rse_gene)$geo_accession

## Extract the sample characteristics
geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)

## Note that the information for this study is a little inconsistent, so we
## have to fix it.
geochar <- do.call(rbind, lapply(geochar, function(x) {
    if('cells' %in% colnames(x)) {
        colnames(x)[colnames(x) == 'cells'] <- 'cell.line'
        return(x)
    } else {
        return(x)
    }
}))

## We can now define some sample information to use
sample_info <- data.frame(
    run = colData(rse_gene)$run,
    group = ifelse(grepl('uninduced', colData(rse_gene)$title), 'uninduced', 'induced'),
    gene_target = sapply(colData(rse_gene)$title, function(x) { strsplit(strsplit(x,
        'targeting ')[[1]][2], ' gene')[[1]][1] }),
    cell.line = geochar$cell.line
)

## Scale counts by taking into account the total coverage per sample
rse <- scale_counts(rse_gene)

## Add sample information for DE analysis
colData(rse)$group <- sample_info$group
colData(rse)$gene_target <- sample_info$gene_target

## Perform differential expression analysis with DESeq2
library('DESeq2')

## Specify design and switch to DESeq2 format
dds <- DESeqDataSet(rse, ~ gene_target + group)

## Perform DE analysis
dds <- DESeq(dds, test = 'LRT', reduced = ~ gene_target, fitType = 'local')
res <- results(dds)

## Explore results
plotMA(res, main="DESeq2 results for SRP009615")

## Expressed-regions example
regions <- expressed_regions('SRP009615', 'chrY', cutoff = 5L, 
    maxClusterGap = 3000L, outdir = file.path(scipath, 'SRP009615'))
regions

## Reproducibility info
library('devtools')
options(width = 120)
session_info()
