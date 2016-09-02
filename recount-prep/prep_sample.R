## Process the data for a given sample and create the data needed for recount

## For specifying parameters
library('getopt')

## Specify parameters
spec <- matrix(c(
	'bigwig_file', 'f', 1, 'character', 'Path to the bigwig file',
    'counts_file', 'c', 1, 'character', 'Path to the counts.tsv.gz file',
    'bwtool', 'b', 2, 'character',
    "Path to bwtool. If not provided, it's assumed that it is on the $PATH",
    'wiggletools', 'w', 2, 'character',
    "Path to wiggletools. If not provided, it's assumed that it is on the $PATH",
    'paired', 'p', 2, 'logical', 'Is the sample paired end?',
    'calculate_auc', 'a', 2, 'logical', 'Whether to calculate the AUC',
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
suppressMessages(library('recount'))
suppressMessages(library('devtools'))

## For testing
if(FALSE) {
    opt <- list('bigwig_file' = '/dcl01/leek/data/sunghee_analysis/processed/coverage_bigwigs/JH-13_GGCTAC_L006.bw',
    'counts_file' = '/dcl01/leek/data/sunghee_analysis/processed/cross_sample_results/counts.tsv.gz',
    'bwtool' = '/dcl01/leek/data/bwtool/bwtool-1.0/bwtool',
    'wiggletools' = 'wiggletools',
    'calculate_auc' = TRUE,
    'paired' = FALSE
    )
}

## Are we on JHPCE? Print some helpful info
jhpce <- grepl('compute-', Sys.info()['nodename'])
if(jhpce & is.null(opt$wiggletools)) {
    message(paste(Sys.time(), 'Note that you can use wiggletools with:
    module load wiggletools/default
    '))
}

## Set some defaults
if(is.null(opt$bwtool)) opt$bwtool <- 'bwtool'
if(is.null(opt$calculate_auc)) opt$calculate_auc <- FALSE
if(is.null(opt$wiggletools)) opt$wiggletools <- 'wiggletools'
if(is.null(opt$paired)) opt$paired <- as.logical(NA)
    
## Print options used
message(paste(Sys.time(), 'options used:'))
print(opt)
    
## Get sample name
bw <- opt$bigwig_file
names(bw) <- gsub('.*/|.bw', '', bw)

## Read sample info
counts <- read.table(opt$counts_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
counts <- counts[counts$X == names(bw), ]

if(!opt$calculate_auc) {
    opt$calculate_auc <- !'auc' %in% colnames(counts)
    if(opt$calculate_auc) {
        message(paste(Sys.time(), 'detected that the AUC information is missing and will calculate this'))
    }
}

if(opt$calculate_auc) {
    ## Choose name for temporary file
    auc_file <- file.path(tempdir(), paste0(names(bw), '.auc'))
    system(paste(opt$wiggletools, 'AUC', auc_file, bw))
    counts$auc <- as.numeric(readLines(auc_file))
    ## Clean up
    unlink(auc_file)
}

## Check that files have been downloaded
bed <- 'ucsc-knowngene-hg38.bed'
gene <- 'ucsc-knowngene-hg38-genes-bp-length.Rdata'
exon <- 'ucsc-knowngene-hg38-exons.Rdata'
count_groups_file <- 'count_groups.Rdata'
if(any(!file.exists(c(bed, gene, exon, count_groups_file)))) {
    message(paste(Sys.time(), 'downloading missing files'))
    source('prep_setup.R')
}

## Define the metadata
metadata <- colData(rse_gene_SRP009615)[1, ]
metadata[1, ] <- NA
rownames(metadata) <- names(bw)
metadata$bigwig_file <- gsub('.*/', '', bw)
metadata$auc <- counts$auc
metadata$mapped_read_count <- as.integer(
    strsplit(counts$total.mapped.reads, ',')[[1]][1])
metadata$reads_downloaded <- as.integer(
    strsplit(counts$total.reads, ',')[[1]][1])
metadata$paired_end <- opt$paired

## Run bwtool to count at the exon level
bw_tsv <- file.path(tempdir(), paste0(names(bw), '.tsv'))
cmd_bwtool <- paste(opt$bwtool, 'summary', bed, bw, 
    "/dev/stdout -fill=0 -with-sum | cut -f1-3,10 | awk -v CONVFMT=%.17g '{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4}' >",
    bw_tsv)
message(paste(Sys.time(), 'running bwtool on the exons'))
message(paste(Sys.time(), 'command used:', cmd_bwtool))
system(cmd_bwtool)

## Read the counts information from bwtool
exon_counts <- mapply(function(tsvFile, sampleName) {
    message(paste(Sys.time(), 'reading file', tsvFile))
    res <- read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL, NULL, 'numeric'))
    colnames(res) <- sampleName
    return(as.matrix(res))
}, bw_tsv, names(bw), SIMPLIFY = FALSE)
exon_counts <- do.call(cbind, exon_counts)

## Save exon counts
dir.create('rse_temp', showWarnings = FALSE)
message(paste(Sys.time(), 'writing file', file.path('rse_temp', 
    paste0('counts_exon_', names(bw), '.tsv'))))
write.table(as.data.frame(exon_counts), file = file.path('rse_temp',
    paste0('counts_exon_', names(bw), '.tsv')),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)


## Load exons info
load(exon)

## Create rse_exon
exons_all <- unlist(exons)
rse_exon <- SummarizedExperiment(assays = list('counts' = exon_counts),
    colData = DataFrame(metadata), rowRanges = exons_all)
message(paste(Sys.time(), 'writing file', file.path('rse_temp', 
    paste0('rse_exon_', names(bw), '.Rdata'))))
save(rse_exon, file = file.path('rse_temp', paste0('rse_exon_', names(bw),
    '.Rdata')))

## Summarize counts at gene level
load(count_groups_file)
load(gene)
counts_gene <- lapply(split(as.data.frame(exon_counts), count_groups), colSums)
counts_gene <- do.call(rbind, counts_gene)
rownames(counts_gene) <- names(genes)

## Save gene counts
message(paste(Sys.time(), 'writing file', file.path('rse_temp', 
    paste0('counts_gene_', names(bw), '.tsv'))))
write.table(as.data.frame(counts_gene), file = file.path('rse_temp',
    paste0('counts_gene_', names(bw), '.tsv')),
    sep = '\t', row.names = FALSE, quote = FALSE,
    col.names = TRUE)

## Create gene level rse
rse_gene <- SummarizedExperiment(assays = list('counts' = counts_gene),
    colData = DataFrame(metadata), rowRanges = genes)
message(paste(Sys.time(), 'writing file', file.path('rse_temp', 
    paste0('rse_gene_', names(bw), '.Rdata'))))
save(rse_gene, file = file.path('rse_temp', paste0('rse_gene_', names(bw),
    '.Rdata')))

## Reproducibility info
proc.time()
options(width = 120)
session_info()



