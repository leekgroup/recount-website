library('SummarizedExperiment')

message(paste(Sys.time(), 'loading exon level RSE object'))
load('/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682/rse_exon.Rdata')

message(paste(Sys.time(), 'loading gene level RSE object'))
load('/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682/rse_gene.Rdata')

fix_tissue <- function(rse) {
    ## Fix 5 samples that are missing smts info but have smtsd info
    colData(rse)$smts[colData(rse)$smts == '' & colData(rse)$smtsd == 'Esophagus - Mucosa'] <-  'Esophagus'
    colData(rse)$smts[colData(rse)$smts == '' & colData(rse)$smtsd == 'Skin - Sun Exposed (Lower leg)'] <-  'Skin'
    colData(rse)$smts[colData(rse)$smts == '' & colData(rse)$smtsd == 'Stomach'] <-  'Stomach'
    
    ## Print number of samples by tissue
    print(sort(table(colData(rse)$smts)))
    
    ## Finish
    return(rse)
}

split_by_tissue <- function(rse, type) {
    ## Split data by tissue
    tissues <- unique(colData(rse)$smts)
    
    message(paste(Sys.time(), 'splitting', type, 'level information by tissue'))
    rse_split <- lapply(tissues, function(tissue) {
        subset(rse, select = colData(rse)$smts == tissue)
    })
    names(rse_split) <- tissues
    
    ## Run some checks
    stopifnot(sum(sapply(rse_split, ncol)) == ncol(rse))
    stopifnot(all(sort(sapply(rse_split, ncol)) - sort(table(colData(rse)$smts)) == 0))
    
    ## For file names
    tissues <- gsub(' ', '_', tolower(tissues))
    
    res <- mapply(function(rse_tissue, tissue, type) {
        
        message(paste(Sys.time(), 'saving', type,
            'level information for tissue', tissue))
        
        rse_file <- paste0('/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682/rse_', type, '_', tissue, '.Rdata')
        
        ## Make sure the variable name is rse_exon or rse_gene
        varname <- paste0('rse_', type)
        assign(varname, rse_tissue)
        
        ## Save
        save(list = varname, file = rse_file)
        return(rse_file)
    }, rse_split, tissues, MoreArgs = list(type = type), SIMPLIFY = FALSE)
}

message(paste(Sys.time(), 'splitting by tissue at the exon level'))
split_by_tissue(fix_tissue(rse_exon), 'exon')

message(paste(Sys.time(), 'splitting by tissue at the gene level'))
split_by_tissue(fix_tissue(rse_gene), 'gene')

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
devtools::session_info()
