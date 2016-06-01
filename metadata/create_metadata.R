## Prepare metadata
# module load R/3.3
# mkdir -p logs
# Rscript create_metadata.R -p "sra" > logs/create_metadata_sra_log.txt 2>&1
# Rscript create_metadata.R -p "gtex" > logs/create_metadata_gtex_log.txt 2>&1

library('getopt')

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
    metadata <- read.table('/dcl01/leek/data/gtex_work/runs/sra/v2/recount2_metadata.tsv',
    header = TRUE, sep = '\t', stringsAsFactors = FALSE
    )
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
    pheno$paired_end <- pheno$librarylayout == 'PAIRED'
    pheno$project <- as.character(pheno$srastudy)
    pheno$sample <- as.character(pheno$sample)
    pheno$experiment <- as.character(pheno$experiment)
    pheno$proportion_of_reads_reported_by_sra_aligned <- pheno$smmaprt
    
    
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
    
    ## Find number of reads (for paired-end samples, that's 2x)
    counts$n_reads <- as.numeric(sapply(strsplit(counts$total.reads, ','),
        '[[', 1))
    
    pheno$read_count_as_reported_by_sra <- counts$n_reads
    pheno$reads_aligned <- pheno$read_count_as_reported_by_sra * pheno$proportion_of_reads_reported_by_sra_aligned
    pheno$mapped_read_count <- counts$totalMapped
    
    ## Mis-reported?
    ## See https://github.com/nellore/runs/blob/93c80b34f9e09c84831d4ffd652c3742dd804487/sra/v2/recount2_metadata.py#L134-L148
    ratio <- pheno$mapped_read_count / pheno$read_count_as_reported_by_sra
    summary(ratio)
    pheno$sra_misreported_paired_end <- ratio == 0.5 | ratio > 1
    
    
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
    
} else {
    stop("Invalid 'project' choice. Use gtex or sra")
}

## Change sharq_* to sharq_beta_*
colnames(metadata)[colnames(metadata) == 'sharq_tissue'] <- 'sharq_beta_tissue'
colnames(metadata)[colnames(metadata) == 'sharq_cell_type'] <- 'sharq_beta_cell_type'

## Find GEO number
find_geo <- function(run) {
    geo <- system(paste0("curl \"http://www.ncbi.nlm.nih.gov/gds?LinkName=sra_gds&from_uid=$(curl \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=", run ,"\" | sed -n 's|[^<]*<Id>\\([^<]*\\)</Id>[^<]*|\\1|gp')\" | grep Series | awk -F 'acc\\\\=GSM' '{print \"GSM\" $2}' | grep -vFx GSM | awk -F '\"' '{print $1}'"), intern = TRUE)
    if(length(geo) == 0) {
        return(NA)
    } else {
        return(geo)
    }
}

## Not all cases have GEO id's, like:
# find_geo('DRR000897')
metadata$geo_accession <- sapply(metadata$run, find_geo)

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

## Save the metadata
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
sapply(metadata, function(x) { sum(is.na(x)) })

print('Percent of NAs per column')
sapply(metadata, function(x) { sum(is.na(x)) }) / nrow(metadata) * 100

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

## Save project ids in a file
write.table(unique(metadata_clean$project), file = paste0('project_ids_',
    opt$project, '.txt'), sep = '\t', row.names = FALSE, quote = FALSE,
    col.names = FALSE)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
