## Download required files for prep_sample.R
library('downloader')
library('devtools')

## Helper function
down_file <- function(file_name) {
    if(!file.exists(file_name)) {
        download(
        paste0('https://github.com/leekgroup/recount-website/raw/master/genes/',
            file_name), dest = file_name, mode = 'wb')
    } else {
        message(paste(Sys.time(), 'using previously downloaded', file_name))
    }
}

## Download files if necessary
bed <- 'ucsc-knowngene-hg38.bed'
gene <- 'ucsc-knowngene-hg38-genes-bp-length.Rdata'
exon <- 'ucsc-knowngene-hg38-exons.Rdata'

down_file(bed)
down_file(gene)
down_file(exon)

## Reproducbility info
options(width = 120)
session_info()
