## Prepare metadata

library('getopt')
library('GenomicRanges')
library('recount')
library('BiocParallel')
library('parallel')
library('XML')

## Specify parameters
spec <- matrix(c(
	'project', 'p', 1, 'character', 'Project ID: sra or gtex',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

if(opt$project == 'sra') {
    ## Load sra info
    ## ignoring GSM column since it has a bit of everything on it, not only GEO
    ## accesion ids
    metadata <- read.table('/dcl01/leek/data/gtex_work/runs/sra/v2/recount2_metadata.tsv',
    header = TRUE, sep = '\t', stringsAsFactors = FALSE, quote = "", colClasses = rep(c('character', 'integer', 'numeric', 'integer', 'numeric', 'character', 'character'), c(4, 2, 1, 3, 1, 5, 1)))
    colnames(metadata) <- tolower(colnames(metadata))
    colnames(metadata) <- gsub('\\.', '_', colnames(metadata))

    ## Change logical variables to TRUE/FALSE
    metadata$paired_end <- as.logical(metadata$paired_end)
    metadata$sra_misreported_paired_end <- as.logical(metadata$sra_misreported_paired_end)
    
    ## Load SRA-run info
    sra <- read.csv('/dcl01/leek/data/gtex_work/runs/sra/v2/hg38/SraRunInfo.csv',
        header = TRUE, stringsAsFactors = FALSE)
    colnames(sra) <- tolower(colnames(sra))

    ## Find average read length
    i <- match(metadata$run, sra$run)
    stopifnot(identical(metadata$run, sra$run[i]))
    metadata$avg_read_length <- sra$avglength[i]
    metadata$geo_accession <- NA
    i <- grepl('GSM', metadata$gsm)
    metadata$geo_accession[i] <- metadata$gsm[i]
    
    ## Fix some column names
    colnames(metadata)[colnames(metadata) == 'reads_aligned'] <- 'reads_downloaded'
    colnames(metadata)[colnames(metadata) == 'proportion_of_reads_reported_by_sra_aligned'] <- 'proportion_of_reads_reported_by_sra_downloaded'
    
    ## Drop GSM
    metadata <- metadata[, -which(colnames(metadata) == 'gsm')]
} else if (opt$project == 'gtex') {
    ## Load SRA metadata (metadata object)
    stopifnot(file.exists('/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata'))
    load('/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata')
    
    ## Load GTEx metadata (pheno object)
    load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_complete.Rdata')
    
    ## Fix GTEx data
    pheno$auc <- pheno$SumCoverage
    colnames(pheno) <- tolower(colnames(pheno))
    pheno$avg_read_length <- pheno$avglength
    pheno$bigwig_path <- pheno$bigwigpath
    pheno$bigwig_file <- gsub('.*coverage_bigwigs/', '', pheno$bigwigpath)
    pheno$paired_end <- TRUE ## They are all paired-end
    pheno$project <- as.character(pheno$srastudy)
    pheno$sample <- as.character(pheno$sample)
    pheno$experiment <- as.character(pheno$experiment)
    
    
    ## Get mapped reads by Rail-RNA
    ## Find count files
    counts_files <- file.path(dir('/dcl01/leek/data/gtex', pattern = 'batch',
        full.names = TRUE), 'cross_sample_results', 'counts.tsv.gz')
    names(counts_files) <- dir('/dcl01/leek/data/gtex', pattern = 'batch')

    ## Read in counts info
    counts <- lapply(counts_files, function(i) {
        read.table(i, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    })
    counts <- do.call(rbind, counts)
    counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads,
        ','), '[[', 1))

    ## Match files to counts
    map <- match(gsub('.bw', '', pheno$bigwig_file), counts$X)
    counts <- counts[map, ]
    stopifnot(identical(counts$X, gsub('.bw', '', pheno$bigwig_file)))
    
    ## Save mapped read counts (from Rail-RNA)
    pheno$mapped_read_count <- counts$totalMapped
    
    ## Find reads_downloaded
    reads_info <- read.csv('/dcl01/leek/data/gtex_work/runs/gtex/SraRunInfo.csv', stringsAsFactors = FALSE)
    map_reads <- match(pheno$run, reads_info$Run)
    ## All GTEx samples are paired end (although SRA mis-repors 18 of them as 
    ## single-end.
    pheno$read_count_as_reported_by_sra <- reads_info$spots[map_reads] * 2
    
    ## For the complete samples, reads_downloaded is the same as
    ## read_count_as_reported_by_sra
    pheno$reads_downloaded <- pheno$read_count_as_reported_by_sra
    pheno$proportion_of_reads_reported_by_sra_downloaded <- 1
    
    ## Label mis-reported samples
    pheno$sra_misreported_paired_end <- reads_info$LibraryLayout[map_reads] == 'SINGLE'
    
    ## Fix the incomplete samples
    incomplete <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/incomplete.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    colnames(incomplete) <- tolower(colnames(incomplete))
    colnames(incomplete) <- gsub('\\.', '_', colnames(incomplete))
    map_inc <- match(incomplete$run, pheno$run)
    
    ## Now actually fix them
    pheno$proportion_of_reads_reported_by_sra_downloaded[map_inc] <-  incomplete$proportion_of_reads_reported_by_sra_aligned
    pheno$reads_downloaded[map_inc] <- incomplete$reads_aligned
    pheno$read_count_as_reported_by_sra[map_inc] <- incomplete$read_count_as_reported_by_sra
    
    
        
    ## Store column names
    sra <- colnames(metadata)
    gtex <- colnames(pheno)
    
    ## Create GTEx metadata
    m <- as.data.frame(matrix(NA, nrow = nrow(pheno), ncol = length(sra)))
    colnames(m) <- sra
    for(i in sra) {
        if(i %in% gtex) m[, i] <- pheno[, i]
    }
    metadata <- m
    
    ## Read in public data available via http://www.gtexportal.org/home/datasets
    ## after you register at the GTEx portal
    meta <- read.table('GTEx_Data_V6_Annotations_SampleAttributesDS.txt',
        header = TRUE, sep = '\t', quote = "", stringsAsFactors = FALSE)
    
    map_meta <- match(pheno$sampid, meta$SAMPID)
    meta <- meta[map_meta, ]
    stopifnot(nrow(meta) == nrow(metadata))
    
    ## Lower case the variable names
    colnames(meta) <- tolower(colnames(meta))
    
    ## Check that none of the metadata columns are duplicated
    stopifnot(all(!colnames(metadata) %in% colnames(meta)))
    
    ## Add the information
    metadata <- cbind(metadata, meta)
} else {
    stop("Invalid 'project' choice. Use gtex or sra")
}

## Change sharq_* to sharq_beta_*
colnames(metadata)[colnames(metadata) == 'sharq_tissue'] <- 'sharq_beta_tissue'
colnames(metadata)[colnames(metadata) == 'sharq_cell_type'] <- 'sharq_beta_cell_type'

## Save project ids in a file
write.table(unique(metadata$project), file = paste0('project_ids_',
    opt$project, '.txt'), sep = '\t', row.names = FALSE, quote = FALSE,
    col.names = FALSE)
    
## Find bigwig files
if(opt$project == 'sra') {
    bigwigs <- system(paste0('cut -f 5 -d " " /dcl01/leek/data/recount-website/bwtool/bwtool_cmds_',
        opt$project, '.txt'), intern = TRUE)
    names(bigwigs) <- gsub('.*coverage_bigwigs/|.bw', '', bigwigs)
    j <- match(metadata$run, names(bigwigs))

    ## Matches number of bigwig files
    stopifnot(sum(is.na(j)) == sum(is.na(metadata$auc)))
    metadata$bigwig_path <- bigwigs[j]
    metadata$bigwig_file <- gsub('.*coverage_bigwigs/', '',
        metadata$bigwig_path)
}


## Locate tsv files
tsv_dir <- ifelse(opt$project == 'sra', '/dcl01/leek/data/recount2/coverage', '/dcl01/leek/data/recount2/coverage_gtex')
tsv <- dir(tsv_dir, pattern = 'tsv', full.names = TRUE)
names(tsv) <- gsub('.sum.tsv', '', dir(tsv_dir, pattern = 'tsv'))
k <- match(metadata$run, names(tsv))

## Number of tsv files matches number of bigwig files
stopifnot(sum(is.na(k)) == sum(is.na(metadata$auc)))
metadata$tsv_path <- tsv[k]

## Function for re-running in case something fails
runMyFun <- function(f, ...) {
    res <- 'trying'
    while(res == 'trying') {
        res <- tryCatch(f(...), error = function(e) {
            Sys.sleep(round(runif(1, min =1 , max = 6), 0))
            return('trying')
        })
        if(is.na(res)) break
    }
    return(res)
}

## Find GEO accession id
## Not all cases have GEO id's, like:
# find_geo('DRR000897')
if(!'geo_accession' %in% colnames(metadata)) {
    if(opt$project == 'gtex') {
        metadata$geo_accession <- NA
    } else {
        bp <- MulticoreParam(workers = 3, outfile = Sys.getenv('SGE_STDERR_PATH'))
        metadata$geo_accession <- unlist(bplapply(metadata$run,
            function(runid) {
                runMyFun(find_geo, run = runid, verbose = TRUE, sleep = 1)
            }, BPPARAM = bp),
        use.names = FALSE)
    }
}

## Save the metadata (backup with geo accession ids)
save(metadata, file = paste0('metadata_', opt$project, '.Rdata'))


## Find some sample information from geo
dir.create(paste0('geo_info_', opt$project), showWarnings = FALSE)
extract_geo <- function(id) {
    if(is.na(id)) {
        res <- DataFrame('title' = NA, 'characteristics' = CharacterList(NA))
    } else {
        info <- runMyFun(geo_info, geoid = id, verbose = TRUE,
            destdir = paste0('geo_info_', opt$project), sleep = 1)
        res <- DataFrame('title' = info$title,
            'characteristics' = info$characteristics)
    }
    return(res)
}
geo <- do.call(rbind, mclapply(metadata$geo_accession, extract_geo, mc.cores = 3))

## Combine information (metadata will now be a DataFrame object)
metadata <- cbind(metadata, geo)

## Save the metadata (final version)
save(metadata, file = paste0('metadata_', opt$project, '.Rdata'))

## Explore metadata
print('Number of dimensions')
dim(metadata)

print('Number of unique project IDs')
length(unique(metadata$project))

print('Number of unique run IDs')
length(unique(metadata$run))

print('First couple of rows')
head(metadata)

print('Number of NAs per column')
sapply(metadata, function(x) { sum(sum(is.na(x))) })

print('Percent of NAs per column')
sapply(metadata, function(x) { sum(sum(is.na(x))) }) / nrow(metadata) * 100

## Save a file per project ID
dir.create(paste0('project_metadata_', opt$project), showWarnings = FALSE)

metadata_clean <- metadata[, !colnames(metadata) %in% c('bigwig_path',
    'tsv_path')]
save(metadata_clean, file = paste0('metadata_clean_', opt$project, '.Rdata'))

meta <- split(metadata_clean, metadata_clean$project)
stopifnot(length(meta) == length(unique(metadata_clean$project)))

xx <- sapply(unique(metadata_clean$project), function(project) {
    project_metadata <- meta[[project]]
    write.table(project_metadata, file.path(paste0('project_metadata_',
        opt$project), paste0(project, '.tsv')), sep = '\t', row.names = FALSE,
        quote = FALSE)
})

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
