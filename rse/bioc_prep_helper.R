## Define some functions
bioc_load <- function(f) {
    load(f)
    return(res)
}
## Function for one run
bioc_run <- function(run, sample, jx_n) {
    message(paste(Sys.time(), 
        'extracting info from jx_project_tab for run', run))
    
    resdir <- '/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682/resdir'
    dir.create(resdir, showWarnings = FALSE)

    sample_reads <- lapply(jx_project_tab, function(jx_proj_tab) {
        subset(jx_proj_tab, sample_id == sample)
    })
    sample_reads <- do.call(rbind, sample_reads)
    if(nrow(sample_reads) == 0)  {
        message(paste(Sys.time(),
            'found no junction counts for run', run))
        next
    }
    jx_map <- match(jx_project$jx_id, sample_reads$jx_id)
    x <- sample_reads$reads[jx_map[!is.na(jx_map)]]
    i <- which(!is.na(jx_map))
    j <- rep(1, length(i))
    stopifnot(length(i) == length(x))
    res <- sparseMatrix(i = i, j = j, x = x, dims = c(jx_n, 1))
    colnames(res) <- run

    ## Write result
    res_file <- file.path(resdir, paste0('res_', run, '.Rdata'))
    save(res, file = res_file)
    return(res_file)
}
