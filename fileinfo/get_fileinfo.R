## Load libraries
library('getopt')
library('tools')
library('Hmisc')
library('GenomicRanges')

## Specify parameters
spec <- matrix(c(
    'project', 'p', 1, 'character', 'Project ID',
	'metadata', 'm', 1, 'character', 'Metadata file name',
	'projectid', 'i', 1, 'character', 'Project ID',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) {
    opt <- list(project = 'sra', projectid = 'DRP000499',
        'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata'
    )
    ## Largest one, to find memory needed
    opt <- list(project = 'sra', projectid = 'SRP025982',
        'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata'
    )
    ## GTEx
    opt <- list(project = 'gtex', projectid = 'SRP012682',
        'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_gtex.Rdata'
    )
}

## Create output dir
dir.create(paste0('fileinfo_', opt$project), showWarnings = FALSE)

## Load metadata
load(opt$metadata)

## Subset to project of interest
metadata <- subset(metadata, project == opt$projectid)
if(nrow(metadata) == 0) stop(paste('Invalid project id', opt$projectid))

## Subset to only use the samples that have tsv files
metadata <- metadata[!is.na(metadata$tsv_path), ]
if(nrow(metadata) == 0) stop(paste('No samples have bwtool tsv files for project', opt$projectid))
rownames(metadata) <- NULL

## Locate counts and rse files
rse_path <- file.path('/dcl01/leek/data/recount-website/rse',
    paste0('rse_', opt$project), opt$projectid)
rse_files <- dir(rse_path, full.names = TRUE)
names(rse_files) <- dir(rse_path)

if(opt$projectid != 'SRP025982') {
    rse_up <- c('counts_exon.tsv.gz', 'counts_gene.tsv.gz', 'rse_exon.Rdata',
        'rse_gene.Rdata', 'counts_jx.tsv.gz', 'rse_jx.Rdata')
} else {
    rse_up <- c('counts_exon.tsv.gz', 'counts_gene.tsv.gz', 'rse_exon.Rdata',
        'rse_gene.Rdata')#, 'rse_jx.Rdata')
}

if(!all(rse_up %in% names(rse_files))) {
    stop(paste('Missing counts/rse files for project', opt$projectid))
}

## Locate metadata file
meta_path <- file.path('/dcl01/leek/data/recount-website/metadata',
    paste0('project_metadata_', opt$project), paste0(opt$projectid, '.tsv'))
names(meta_path) <- paste0(opt$projectid, '.tsv')
if(!file.exists(meta_path)) stop(paste('Missing metadata file for project', opt$projecid))

## Find mean bigwig file
mean_bw <- file.path('/dcl01/leek/data/recount-website/mean', 
    paste0('means_', opt$project), paste0('mean_', opt$projectid, '.bw'))
if(!file.exists(mean_bw)) stop(paste('Missing mean bigwig file for project', opt$projecid))
names(mean_bw) <- paste0('mean_', opt$projectid, '.bw')

## Create output dir for the project
outdir <- file.path('/dcl01/leek/data/recount-website/fileinfo', 
    paste0('fileinfo_', opt$project), opt$projectid)
dir.create(outdir, showWarnings = FALSE)
stopifnot(dir.exists(outdir))

## Find bigwig files
upload_bigwig <- metadata$bigwig_path
names(upload_bigwig) <- metadata$bigwig_file


## Build paths for junction raw files
jx_raw <- c(
    file.path('/dcl01/leek/data/recount_junctions',
        paste0(opt$projectid, '.junction_id_with_transcripts.bed.gz')),
    file.path('/dcl01/leek/data/recount_junctions',
        paste0(opt$projectid, '.junction_coverage.tsv.gz'))
)
names(jx_raw) <- c(
    paste0(opt$projectid, '.junction_id_with_transcripts.bed.gz'),
    paste0(opt$projectid, '.junction_coverage.tsv.gz')
)

## List of files to upload
upload_files <- c(rse_files[rse_up], meta_path, mean_bw, upload_bigwig, jx_raw)

## Compute md5sums
print('Time spent calculating md5sums')
system.time( files_md5sum <- md5sum(upload_files) )

## Calculate total file size
print('Time spent calculating file sizes')
fileSize <- function(files) {
    as.numeric(system(paste('du -l', paste(files, collapse = ' '), 
        '| cut -f1'), intern = TRUE))
}
up_sets <- cut2(seq_len(length(upload_files)), m = 1000)
up_files <- split(upload_files, up_sets)

system.time( file_size <- unlist(lapply(up_files, fileSize)) )
names(file_size) <- names(upload_files)
project_size <- sum(file_size / 1024)
if(project_size < 1024) {
    print(paste('The total file size for project', opt$projectid, 'is', round(project_size, digits = 3), 'Mb.'))
} else {
    print(paste('The total file size for project', opt$projectid, 'is', round(project_size / 1024, digits = 3), 'Gb.'))
}

## Save md5sums and file size.
## The only thing missing is the actual file URL, but that can only be completed
## after uploading the files.
write.table(data.frame(file = names(upload_files), md5sum = files_md5sum,
    size = file_size, stringsAsFactors = FALSE),
    file = file.path(outdir, 'files_info.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)
    
## Save final list of files to upload
upload_files <- c(upload_files, 'files_info.tsv' = file.path(outdir,
    'files_info.tsv'))
save(upload_files, file = file.path(outdir, 'upload_files.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
