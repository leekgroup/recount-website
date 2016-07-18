library('SummarizedExperiment')
library('Hmisc')

## Locate rse files
files_locate <- function(type) {
    f <- dir('/dcl01/leek/data/recount-website/rse/rse_sra',
        pattern = paste0('rse_', type, '.Rdata'), recursive = TRUE,
        full.names = TRUE)
    names(f) <- dir('/dcl01/leek/data/recount-website/rse/rse_sra',
        pattern = paste0('rse_', type, '.Rdata'), recursive = TRUE)
    return(f)
}


files_load <- function(f, type) {
    message(paste(Sys.time(), 'loading file', f))
    load(f)
    if(type == 'gene') {
        return(rse_gene)
    } else if (type == 'exon'){
        return(rse_exon)
    } else if (type == 'jx'){
        return(rse_jx)
    }
}

files_main <- function(type) {
    message(paste(Sys.time(), 'locating files'))
    rse_files <- files_locate(type)
    
    ## Split files into groups of ~20
    rse_files_list <- split(rse_files, cut2(seq_len(length(rse_files)), m = 20))
    
    message(paste(Sys.time(), 'loading files'))
    rse_list <- lapply(rse_files_files, function(x) {
        res <- lapply(x, files_load, type = type)
        message(paste(Sys.time(), 'merging RSE objects (group)'))
        do.call(cbind, res)
    })
    
    message(paste(Sys.time(), 'saving rse_list'))
    save(rse_list, file = paste0(
        '/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_list_',
        type, '.Rdata'))
    
    message(paste(Sys.time(), 'merging RSE objects from groups'))
    varname <- paste0('rse_', type)
    assign(varname, do.call(cbind, rse_list))
    
    ## Save results
    dir.create('/dcl01/leek/data/recount-website/rse/rse_sra/all',
        showWarnings = FALSE)
    rse_file <- paste0('/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_',
        type, '.Rdata')
    save(list = varname, file = rse_file)
    return(rse_file)
}

message(paste(Sys.time(), 'processing exon files'))
files_main('exon')

message(paste(Sys.time(), 'processing gene files'))
files_main('gene')

## Junctions won't work asis since not all the rows are the same, plus it'd be huge
#message(paste(Sys.time(), 'processing exon files'))
#files_main('jx')

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
devtools::session_info()
