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
    
    
    rse_list_file <- paste0(
        '/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_list_',
        type, '.Rdata')
    rse_list_sets_file <- paste0(
        '/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_list_sets_',
        type, '.Rdata')
    dir.create('/dcl01/leek/data/recount-website/rse/rse_sra/all',
        showWarnings = FALSE)
    
    if(!file.exists(rse_list_sets_file)) {
        if(!file.exists(rse_list_file)) {
            message(paste(Sys.time(), 'loading individual RSE files at the', type,
                'level'))
            rse_list <- lapply(rse_files_list, function(x) {
                res <- lapply(x, files_load, type = type)
                message(paste(Sys.time(), 'merging RSE objects (group)'))
                do.call(cbind, res)
            })
    
            message(paste(Sys.time(), 'saving rse_list'))
            save(rse_list, file = rse_list_file)
        } else {
            message(paste(Sys.time(),
                'loading previously computed rse_list at the', type, 'level'))
            load(rse_list_file)
        }
        
        message(paste(Sys.time(),
            'merging RSE objects by sets of groups at the', type, 'level'))
        group.sets <- cut2(seq_len(length(rse_list)), m = 5)
        rse_list_sets <- lapply(levels(group.sets), function(group) {
            do.call(cbind, rse_list[group.sets == group])
        })
        message(paste(Sys.time(), 'saving rse_list_sets'))
        save(rse_list_sets, file = rse_list_sets_file)
        rm(rle_list)
    } else {
        message(paste(Sys.time(),
            'loading previously computed rse_list_sets at the', type, 'level'))
        load(rse_list_sets_file)
    }
    
    message(paste(Sys.time(), 'merging RSE objects from sets at the', type,
        'level'))
    ## Extract data from sets
    message(paste(Sys.time(), 'extracting information from sets'))
    col_data <- do.call(rbind, lapply(rse_list_sets, colData))
    row_ranges <- rowRanges(rse_list_sets[[1]])
    counts_list <- lapply(rse_list_sets, function(x) assays(x)$counts)
    rm(rse_list_sets)
    
    ## Prepare counts matrix
    message(paste(Sys.time(), 'preparing the counts matrix'))
    counts_n <- sapply(counts_list, ncol)
    counts_adj <- c(0, cumsum(counts_n[1:(length(counts_n) - 1)]))
    counts_idx <- mapply(function(n, adj) { seq_len(n) + adj}, counts_n,
        counts_adj, SIMPLIFY = FALSE)
    
    ## Initialize count matrix
    message(paste(Sys.time(), 'initializing counts matrix'))
    counts <- matrix(0, nrow = nrow(counts_list[[1]]),
        ncol = sum(counts_n), dimnames = list(rownames(counts_list[[1]]),
        unlist(lapply(counts_list, colnames))))
    
    ## Fill in matrix
    for(i in seq_len(length(counts_idx))) {
        message(paste(Sys.time(), 'filling in counts matrix with set', i))
        counts[, counts_idx[[i]]] <- counts_list[[i]]
    }
    
    message(paste(Sys.time(), 'creating final rse object'))
    varname <- paste0('rse_', type)
    assign(varname,
        SummarizedExperiment(assays = list('counts' = counts),
            rowRanges = row_ranges, colData = col_data)
    )

    
    ## Save results
    message(paste(Sys.time(), 'saving the final rse object'))
    rse_file <- paste0('/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_',
        type, '.Rdata')
    save(list = varname, file = rse_file)
    return(rse_file)
}

if(!file.exists('rse_gene.Rdata')) {
    message(paste(Sys.time(), 'processing gene files'))
    files_main('gene')
} else {
    message(paste(Sys.time(), 'rse_gene.Rdata already exists'))
}


if(!file.exists('rse_exon.Rdata')) {
    message(paste(Sys.time(), 'processing exon files'))
    files_main('exon')
} else {
    message(paste(Sys.time(), 'rse_exon.Rdata already exists'))
}

## Junctions won't work asis since not all the rows are the same, plus it'd be huge
#message(paste(Sys.time(), 'processing exon files'))
#files_main('jx')

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
devtools::session_info()
