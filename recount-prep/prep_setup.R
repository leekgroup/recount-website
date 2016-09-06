## Download required files for prep_sample.R
suppressMessages(library('downloader'))
suppressMessages(library('devtools'))

## Check that output files are not present
stopifnot(!dir.exists('rse_temp'))
stopifnot(!any(file.exists('rse_gene.Rdata', 'rse_exon.Rdata', 'rse_jx.Rdata')))
stopifnot(!any(file.exists('counts_exon.tsv.gz', 'counts_gene.tsv.gz')))

## Helper function
down_file <- function(file_name, repo = 'recount-website') {
    if(!file.exists(file_name)) {
        url <- ifelse(repo == 'recount-website',
           'https://github.com/leekgroup/recount-website/raw/master/genes/',
           ifelse(repo == 'recount-introns',
               'https://github.com/leekgroup/recount-website/raw/master/rse/',
               'https://github.com/nellore/runs/raw/master/gtex/')
        )
        download(paste0(url, file_name), dest = file_name, mode = 'wb')
    } else {
        message(paste(Sys.time(), 'using previously downloaded', file_name))
    }
}

## Download files if necessary
bed <- 'ucsc-knowngene-hg38.bed'
gene <- 'ucsc-knowngene-hg38-genes-bp-length.Rdata'
exon <- 'ucsc-knowngene-hg38-exons.Rdata'
count_groups_file <- 'count_groups.Rdata'
hg38 <- 'hg38.sizes'
introns <- 'introns_unique.Rdata'

down_file(bed)
down_file(gene)
down_file(exon)
down_file(count_groups_file)
down_file(hg38, repo = 'runs')
down_file(introns, repo = 'recount-introns')

## Reproducibility info
proc.time()
options(width = 120)
session_info()
