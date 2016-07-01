## Load libraries
library('getopt')
suppressPackageStartupMessages(library('GenomicRanges'))
suppressPackageStartupMessages(library('SummarizedExperiment'))
suppressPackageStartupMessages(library('rtracklayer'))
suppressPackageStartupMessages(library('TxDb.Hsapiens.UCSC.hg38.knownGene'))

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
    opt <- list(project = 'sra', 'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata',
        projectid = 'DRP000499')
        
    ## Debugging
    opt <- list(project = 'sra', 'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata',
            projectid = 'DRP000366')
    opt <- list(project = 'sra', 'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata',
            projectid = 'DRP000987')
    ## Largest one, to find memory needed
    opt <- list(project = 'sra', 'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata',
        projectid = 'SRP025982')
}

## Create output dir
dir.create(paste0('rse_', opt$project), showWarnings = FALSE)

## Load GRanges and metadata
load('/dcl01/leek/data/recount-website/genes/ucsc-knowngene-hg38-genes-bp-length.Rdata')
load('/dcl01/leek/data/recount-website/genes/ucsc-knowngene-hg38-exons.Rdata')
load('/dcl01/leek/data/recount-website/genes/count_groups.Rdata')
load(opt$metadata)

## Subset to project of interest
metadata <- subset(metadata, project == opt$projectid)
if(nrow(metadata) == 0) stop(paste('Invalid project id', opt$projectid))

## Subset to only use the samples that have tsv files
metadata <- metadata[!is.na(metadata$tsv_path), ]
if(nrow(metadata) == 0) stop(paste('No samples have bwtool tsv files for project', opt$projectid))
rownames(metadata) <- NULL

## Create output dir for the project
outdir <- paste0('rse_', opt$project, '/', opt$projectid)
dir.create(outdir, showWarnings = FALSE)

## Read counts from bwtool tsv output files
counts <- mapply(function(tsvFile, sampleName) {
    message(paste(Sys.time(), 'reading file', tsvFile))
    res <- read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL, NULL, 'numeric'))
    colnames(res) <- sampleName
    return(as.matrix(res))
}, metadata$tsv_path, metadata$run, SIMPLIFY = FALSE)
counts <- do.call(cbind, counts)

## Memory used by counts
print('Memory used by exon counts')
print(object.size(counts), units = 'Mb')
save(counts, file = file.path(outdir, 'counts_exon.Rdata'))

## Save exon counts
message(paste(Sys.time(), 'writing file', file.path(outdir, 'counts_exon.tsv')))
write.table(as.data.frame(counts), file = file.path(outdir, 'counts_exon.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)
system(paste('gzip', file.path(outdir, 'counts_exon.tsv')))

## Remove bigwig and tsv file paths
metadata_clean <- metadata[, !colnames(metadata) %in% c('bigwig_path',
    'tsv_path')]

## Create exon level rse
exons_all <- unlist(exons)
rse_exon <- SummarizedExperiment(assays = list('counts' = counts),
    colData = DataFrame(metadata_clean), rowRanges = exons_all)
message(paste(Sys.time(), 'writing file', file.path(outdir, 'rse_exon.Rdata')))
save(rse_exon, file = file.path(outdir, 'rse_exon.Rdata'))

## Summarize counts at gene level
counts_gene <- lapply(split(as.data.frame(counts), count_groups), colSums)
counts_gene <- do.call(rbind, counts_gene)
rownames(counts_gene) <- names(genes)

## Memory used by counts at gene level
print('Memory used by gene counts')
print(object.size(counts_gene), units = 'Mb')
save(counts_gene, file = file.path(outdir, 'counts_gene.Rdata'))

## Save gene counts
message(paste(Sys.time(), 'writing file', file.path(outdir, 'counts_gene.tsv')))
write.table(as.data.frame(counts_gene), file = file.path(outdir,
    'counts_gene.tsv'), sep = '\t', row.names = FALSE, quote = FALSE,
    col.names = TRUE)
system(paste('gzip', file.path(outdir, 'counts_gene.tsv')))

## Create gene level rse
rse_gene <- SummarizedExperiment(assays = list('counts' = counts_gene),
    colData = DataFrame(metadata_clean), rowRanges = genes)
message(paste(Sys.time(), 'writing file', file.path(outdir, 'rse_gene.Rdata')))
save(rse_gene, file = file.path(outdir, 'rse_gene.Rdata'))


## Load junctions sample information
message(paste(Sys.time(), 'loading junctions sample information'))
jx_samples <- read.table('/dcl01/leek/data/recount_junctions/sample_ids.tsv',
    sep = '\t', col.names = c('sample_id', 'project', 'run'),
    stringsAsFactors = FALSE, colClasses = 'character')

message(paste(Sys.time(), 'processing project junctions read information'))
## Load project junctions info
jx_project <- read.table(file.path('/dcl01/leek/data/recount_junctions',
    paste0(opt$projectid, '.junction_coverage.tsv.gz')), sep = '\t',
    col.names = c('jx_id', 'sample_ids', 'reads'), stringsAsFactors = FALSE,
    colClasses = 'character')

## Create a table with 1 row per sample for a given junction
jx_project_samples <- strsplit(jx_project$sample_ids, ',')
jx_project_reads <- strsplit(jx_project$reads, ',')
stopifnot(identical(elementNROWS(jx_project_samples),
    elementNROWS(jx_project_reads)))
jx_project_tab <- data.frame(
    jx_id = rep(jx_project$jx_id, elementNROWS(jx_project_samples)),
    sample_id = unlist(jx_project_samples),
    reads = as.numeric(unlist(jx_project_reads))
)
rm(jx_project_samples, jx_project_reads)


message(paste(Sys.time(), 'creating junction counts table'))
## Create junction counts table
jx_counts <- matrix(0, ncol = nrow(metadata_clean), nrow = nrow(jx_project))
colnames(jx_counts) <- metadata_clean$run

## Fill in table
for(run in metadata_clean$run) {
    sample <- jx_samples$sample_id[jx_samples$run == run]
    sample_reads <- subset(jx_project_tab, sample_id == sample)
    if(nrow(sample_reads) == 0)  {
        message(paste(Sys.time(), 'found no junction counts for run', run))
        next
    }
    jx_map <- match(jx_project$jx_id, sample_reads$jx_id)
    jx_counts[!is.na(jx_map), run] <- sample_reads$reads[jx_map[!is.na(jx_map)]]
}
rm(sample, sample_reads, jx_map)

print('Memory used by junction counts')
print(object.size(jx_counts), units = 'Mb')

## Save junction counts
message(paste(Sys.time(), 'writing file', file.path(outdir, 'counts_jx.tsv')))
write.table(as.data.frame(jx_counts), file = file.path(outdir,
    'counts_jx.tsv'), sep = '\t', row.names = FALSE, quote = FALSE,
    col.names = TRUE)
system(paste('gzip', file.path(outdir, 'counts_jx.tsv')))


## Load junction bed file
message(paste(Sys.time(), 'loading the junctions bed file'))
jx_bed <- import.bed(file.path('/dcl01/leek/data/recount_junctions',
    paste0(opt$projectid, '.junction_id_with_transcripts.bed.gz')))

message(paste(Sys.time(), 'parsing the junctions information'))
jx_bed_name <- strsplit(jx_bed$name, '\\|')
jx_bed$junction_id <- sapply(jx_bed_name, '[[', 1)
stopifnot(identical(jx_bed$junction_id, as.character(jx_project$jx_id)))
rm(jx_project)

parse_bed_name <- function(pattern = 'D:', slot = 2) {
    CharacterList(
        lapply(
            strsplit(gsub(pattern, '', sapply(jx_bed_name, '[[', slot)), ';'),
                function(y) { 
                    if(y[1] == 'NA') return(NA) else return(y)
                }
        )
    )
}
jx_bed$found_donor <- parse_bed_name('D:', slot = 2)
jx_bed$found_acceptor <- parse_bed_name('A:', slot = 3)
jx_bed$found_junction <- parse_bed_name('J:', slot = 4)

mcols(jx_bed) <- mcols(jx_bed)[, c('junction_id', 'found_donor',
    'found_acceptor', 'found_junction')]

## Fix seqlengths, have to use data from web for chrEBV
chr_info <- read.table(
    'https://raw.githubusercontent.com/nellore/runs/master/gtex/hg38.sizes',
    sep = '\t', col.names = c('chr', 'length'), stringsAsFactors = FALSE)
chrs <- chr_info$length
names(chrs) <- chr_info$chr
seqlengths(jx_bed) <- chrs[names(seqlengths(jx_bed))]

## Find all transcripts
message(paste(Sys.time(), 'setup for identifying gene ids for transcripts'))
transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene,
    columns = c('tx_name', 'gene_id'))

## Make a table with all the transcript names
trans_names <- DataFrame(
    jx_id = c(
        rep(jx_bed$junction_id, elementNROWS(jx_bed$found_donor)),
        rep(jx_bed$junction_id, elementNROWS(jx_bed$found_acceptor)),
        rep(jx_bed$junction_id, elementNROWS(jx_bed$found_junction))
    ),
    name = c(
        unlist(jx_bed$found_donor),
        unlist(jx_bed$found_acceptor),
        unlist(jx_bed$found_junction)
    ),
    gene_id = CharacterList(NA)
)

## Find ids for unique names
unique_names <- DataFrame(
    name = unique(trans_names$name[!is.na(trans_names$name)]),
    gene_id = CharacterList(NA)
)
map_gene <- match(unique_names$name, transcripts$tx_name)
unique_names$gene_id[!is.na(map_gene)] <- transcripts$gene_id[map_gene[!is.na(map_gene)]]
unique_names$gene_id[elementNROWS(unique_names$gene_id) == 0] <- CharacterList(NA)

## Merge back results to large table
map_names <- match(trans_names$name[!is.na(trans_names$name)],
    unique_names$name)
trans_names$gene_id[!is.na(trans_names$name)] <- unique_names$gene_id[map_names]

## Make table smaller by removing NAs
trans_names <- trans_names[any(!is.na(trans_names$gene_id)), ]

## Initialize the gene ids
jx_bed$gene_ids <- CharacterList(NA)

## Find gene ids
find_gene <- function(jx_id) {
    gene_ids <- trans_names$gene_id[trans_names$jx_id == jx_id]
    gene_ids <- gene_ids[sapply(gene_ids, function(x) !is.na(x))]
    if(length(gene_ids) == 0) {
        res <- CharacterList(NA)
    } else {
        res <- CharacterList(unique(do.call(c, gene_ids)))
    }
    return(res)
}
message(paste(Sys.time(), 'finding the gene ids for each transcript'))
map_jx <- match(jx_bed$junction_id, trans_names$jx_id)
system.time( jx_bed$gene_ids[!is.na(map_jx)] <- do.call(c,
    lapply(jx_bed$junction_id[!is.na(map_jx)], find_gene)) )

## Create the junctions level rse
rse_jx <- SummarizedExperiment(assays = list('counts' = jx_counts),
    colData = DataFrame(metadata_clean), rowRanges = jx_bed)
message(paste(Sys.time(), 'writing file', file.path(outdir, 'rse_jx.Rdata')))
save(rse_jx, file = file.path(outdir, 'rse_jx.Rdata'))


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
