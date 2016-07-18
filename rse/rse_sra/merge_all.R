library('SummarizedExperiment')

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
    
    message(paste(Sys.time(), 'loading files'))
    rse_list <- lapply(rse_files, files_load, type = type)
    
    message(paste(Sys.time(), 'merging RSE objects'))
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
