## Merge the RSE exon and gene objects and create the jx RSE object

## For specifying parameters
library('getopt')

## Specify parameters
spec <- matrix(c(
    'bigwith_path', 'b', 1, 'character',
    'Path to the directory with the bigwig files'
	'jx_path', 'j', 1, 'character',
    'Path to the directory with the junction files',
    'calculate_mean', 'm', 2, 'logical', 'Whether to calculate the AUC',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Load libraries
suppressMessages(library('SummarizedExperiment'))
suppressMessages(library('Hmisc'))
suppressMessages(library('devtools'))

## For testing
if(FALSE) {
    opt <- list(
        'bigwig_path' = '/dcl01/leek/data/sunghee_analysis/processed/coverage_bigwigs/'
        'jx_path' = '/dcl01/leek/data/sunghee_analysis/processed/junctions_and_indels',
        'calculate_mean' = TRUE
    )
}

## Set some defaults
if(is.null(opt$calculate_mean)) opt$calculate_mean <- TRUE

## Helper function for loading rse file
load_rse <- function(rse_file, type = 'exon') {
    message(paste(Sys.time(), 'loading file', rse_file))
    load(rse_file)
    if(type == 'exon') {
        return(rse_exon)
    } else if (type == 'gene') {
        return(rse_gene)
    }
}

## Locate exon rse objects, load them, merge them and save results
exon_files <- dir('rse_temp', 'rse_exon_', full.names = TRUE)
rse_exon <- do.call(cbind, lapply(exon_files, load_rse))
message(paste(Sys.time(), 'saving rse_exon.Rdata'))
save(rse_exon, file = 'rse_exon.Rdata')

## Same for gene rse objects
gene_files <- dir('rse_temp', 'rse_gene_', full.names = TRUE)
rse_gene <- do.call(cbind, lapply(gene_files, load_rse))
message(paste(Sys.time(), 'saving rse_gene.Rdata'))
save(rse_gene, file = 'rse_gene.Rdata')

## Calculate the mean bigwig if necessary
if(opt$calculate_mean) {
    ## Get metadata information
    metadata <- colData(rse_gene)
    
    ## Name resulting mean.bw file
    outbw <- 'mean.bw'
    outwig <- 'mean.wig'
    

    scaleWig <- function(m) {
        paste(paste('scale', round(1e6*100*40 / m$auc, digits = 17),
            file.path(opt$bigwig_path, m$bigwig_file)), collapse = ' ')
    }
    runCmd <- function(cmd, i = NULL) {
        if(is.null(i)) {
            shell_name <- paste0('.createWig.sh')
        } else {
            shell_name <- paste0('.createWig_part', i,
                '.sh')
        }    
        cat(cmd, file = shell_name)
        system(paste('sh', shell_name))
    }

    ## Calculate mean bigwig
    if(nrow(metadata) < 100) {
        ## Scale commands
        cmd <- scaleWig(metadata)
        ## Calculate mean wig file
        message(paste(Sys.time(), 'creating file', outwig))
        cmd <- paste('wiggletools write', outwig, 'mean', cmd)
        system.time( runCmd(cmd) )
    } else {
        ## Define subsets to work on
        sets <- cut2(seq_len(nrow(metadata)), m = 50)
        meta <- split(metadata, sets)
        names(meta) <- seq_len(length(meta))
    
        ## Calculate sums per subsets
        system.time( tmpfiles <- mapply(function(m, i) {
            cmd <- scaleWig(m)
        
            tmpdir <- tempdir()     
            tmpwig <- file.path(tmpdir, paste0('sum_part', i, '.wig'))
            message(paste(Sys.time(), 'creating file', tmpwig))
            cmd <- paste('wiggletools write', tmpwig, 'sum', cmd)
            runCmd(cmd, i)
            return(tmpwig)
        }, meta, names(meta)) )
    
        ## Calculate final mean
        cmd <- paste('wiggletools write', outwig, 'scale', 1/nrow(metadata),
            'sum', paste(tmpfiles, collapse = ' '))
        system.time( runCmd(cmd) )
        
        ## Clean up
        sapply(tmpfiles, unlink)
    }

    ## Transform to bigwig file
    message(paste(Sys.time(), 'creating file', outbw))
    cmd2 <- paste('wigToBigWig', outwig, 'hg38.sizes', outbw)
    system.time( system(cmd2) )
    
    ## Clean up
    unlink(outwig)
}


## Add code for rse_jx.Rdata

## Clean up
to_clean <- c('ucsc-knowngene-hg38.bed',
    'ucsc-knowngene-hg38-genes-bp-length.Rdata',
    'ucsc-knowngene-hg38-exons.Rdata', 'count_groups.Rdata', 'hg38.sizes')
sapply(to_clean, unlink)

## Reproducibility info
proc.time()
options(width = 120)
session_info()
