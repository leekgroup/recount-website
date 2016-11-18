library('jsonlite')
library('recount')
library('BiocParallel')
library('devtools')

## Load TCGA metadata from the JSON file
metaraw <- fromJSON(
    '/dcl01/leek/data/gtex_work/runs/tcga/tcga_metadata.json.gz',
    flatten = TRUE)
    
## Parse out
parse <- function(var, df = metaraw$data$hits) {
    ## Find which are NULL
    i.null <- which(sapply(df[, var], is.null))
    
    if(length(i.null) == 0) {
        i.null <- which(sapply(df[, var], is.logical))
    }
    
    if(length(i.null) > 0) {
        ## Make NAs the ones that are missing
        i.df <- df[, var][which(!sapply(df[, var], is.null))[1]]
        i.df[[1]][1, ] <- NA
        i.list <- rep(i.df, length(i.null))
        
        ## Make a list with the complete data
        ori.list <- df[, var]
        ori.list[i.null] <- i.list
        
        ## rbind
        res <- rbind.pages(ori.list)
    } else {
        res <- rbind.pages(df[, var])
    }        
    
    colnames(res) <- paste0(var, '.', colnames(res))
    return(res)
}

parse_nested <- function(var, df) {
    new <- parse(var, df)
    res <- cbind(df[, colnames(df) != var], new)
}

meta.special <- c('metadata_files', 'cases', 'acl', 'associated_entities')
stopifnot(all(colnames(metaraw$data$hits)[
    sapply(metaraw$data$hits, is.list)] %in% c(meta.special,
        'annotations')))

mfiles <- parse('metadata_files')
mfiles$idvar <- rep(seq_len(nrow(mfiles) / 3), each = 3)
mfiles <- reshape(mfiles, timevar = 'metadata_files.data_type',
    direction = 'wide', idvar = 'idvar')[, -1]
colnames(mfiles) <- gsub(' ', '_', tolower(colnames(mfiles)))
colnames(mfiles) <- gsub('_metadata', '', colnames(mfiles))

cases <- parse('cases')
cases.special <- c('cases.diagnoses', 'cases.samples', 'cases.exposures')

cases.vars <- lapply(cases.special, parse, df = cases)
names(cases.vars) <- cases.special

## Fix nested case for diagnoses
cases.vars$cases.diagnoses <- parse_nested('cases.diagnoses.treatments',
    cases.vars$cases.diagnoses)

## Fix nested cases for samples
cases.vars$cases.samples <- parse_nested('cases.samples.portions',
    cases.vars$cases.samples)
cases.vars$cases.samples <- parse_nested('cases.samples.portions.analytes',
    cases.vars$cases.samples)
cases.vars$cases.samples <- parse_nested(
    'cases.samples.portions.analytes.aliquots', cases.vars$cases.samples)

cases.vars$cases.samples <- parse_nested('cases.samples.portions.annotations',
    cases.vars$cases.samples)

cases.vars$cases.samples <- parse_nested('cases.samples.annotations',
    cases.vars$cases.samples)


stopifnot(all(sapply(cases.vars, nrow) == nrow(cases)))

cases <- cases[, !colnames(cases) %in% cases.special]
## Avoid longer names
names(cases.vars) <- NULL
cases <- cbind(cases, cases.vars)

acl <- data.frame('acl' = unlist(metaraw$data$hits$acl), stringsAsFactors = FALSE)
assoc <- parse('associated_entities')
    
meta <- metaraw$data$hits[, !colnames(metaraw$data$hits) %in% meta.special]
meta <- cbind(meta, mfiles, cases, acl, assoc)

## Clean up workspace
rm(mfiles, cases, acl, assoc, cases.vars, cases.special, metaraw, meta.special)

are.list <- function(df) { colnames(df)[sapply(df, is.list)] }
print('The following variables are still lists but do not match the 11285 rows expected and will be left as lists')
lapply(are.list(meta), function(var) { dim(parse(var, meta)) })

sapply(are.list(meta), function(var) { length(meta[, var]) })




## Load clinical data obtained with TCGAbiolinks
load('/dcl01/leek/data/recount-website/metadata/tcga_prep/clin_all.Rdata')



## Load metadata from CGC queries
cgc <- read.table(
    '/dcl01/leek/data/gtex_work/runs/tcga/all_cgc_metadata.tsv.gz',
    header = TRUE, sep = '\t', comment.char = '', quote = '',
    stringsAsFactors = FALSE)

## Set column names depending on the data source
colnames(meta) <- paste0('gdc_', tolower(colnames(meta)))
colnames(clin_all) <- paste0('xml_', tolower(colnames(clin_all)))
colnames(cgc) <- paste0('cgc_', tolower(colnames(cgc)))




## Find the variables that are unique in all given data.frames and
## match the number of rows.
uniq <- function(df) {
    y <- sapply(df, function(x) { length(unique(x))})
    y[y == nrow(df)]
}
print('Unique columns by data source')
uniq(meta)
uniq(clin_all)
uniq(cgc)

head(meta[, names(uniq(meta))])
head(clin_all[, names(uniq(clin_all))])
head(cgc[, names(uniq(cgc))])

## Find the column containing the given key
find_col <- function(key) {
    y <- sapply(meta, function(x) { match(tolower(key), tolower(x))})
    y[!is.na(y)]
}

## Find columns to match against
find_col('1BBD421E-0B42-4760-9E9D-4819474ECF15')
find_col('5452FFD8-0408-4A78-B2FD-F1E0B0D05276')


## Merge the data
cgc$cgc_gdc_uuid <- tolower(cgc$cgc_gdc_uuid)
metadata <- merge(meta, cgc, by.x = 'gdc_file_id', by.y = 'cgc_gdc_uuid',
    all.x = TRUE)
clin_all$xml_bcr_patient_uuid <- tolower(clin_all$xml_bcr_patient_uuid)
metadata <- merge(metadata, clin_all, by.x = 'gdc_cases.case_id',
    by.y = 'xml_bcr_patient_uuid', all.x = TRUE)
print('Dimensions of the merged data')
dim(metadata)

## Remove columns that are all NA
nas <- sapply(metadata, function(x) { sum(is.na(x)) == nrow(metadata) })
print('Number of columns that are all NAs')
sum(nas)

metadata <- metadata[, nas == FALSE]
print('Dimensions of the data without variables made up only of NAs')
dim(metadata)


## Find the columns that are identical
ident <- function(df) {
    result <- lapply(seq_len(ncol(df)), function(i) {
        res <- lapply(seq_len(ncol(df) - i) + i, function(j) {
            data.frame(i, j, ident = identical(df[, i], df[, j]))
        })
        do.call(rbind, res)
    })
    final <- do.call(rbind, result)
    subset(final, ident == TRUE)
}

## Search only within the same class
cl <- sapply(metadata, class)
print('Removing columns that are duplicated')
table(cl)
print('logical columns')
ident(metadata[, cl == 'logical'])
colnames(metadata[, cl == 'logical'][c(2, 3)])
remove <- c('gdc_cases.samples.portions.is_ffpe')

print('numerical columns')
ident(metadata[, cl == 'numeric'])
colnames(metadata[, which(cl == 'numeric')[c(12, 27, 15, 28)] ])
remove <- c(remove, 'xml_weight', 'xml_height')

print('integer, list columns')
ident(metadata[, cl == 'integer'])
ident(metadata[, cl == 'list'])

print('character columns')
m.char <- metadata[, cl == 'character']
for(i in seq_len(ncol(m.char))) m.char[, i] <- tolower(m.char[, i])
i.char <- ident(m.char)
i.char
print('Character columns that are repeated (see other list below)')
colnames(m.char[, unique(i.char$j)])
print('Character columns that are duplicates of the previous list')
colnames(m.char[, unique(i.char$i[!i.char$i %in% i.char$j]) ])

remove <- c(remove, colnames(m.char[, unique(i.char$j)]))

## Final merged metadata
metadata <- metadata[, !colnames(metadata) %in% remove]
print('Dimensions of the final merged metadata')
dim(metadata)

## Clean up workspace
rm(cgc, clin_all, meta, cl, remove, i.char, m.char, nas, i)


## Find count files
counts_files <- file.path(dir('/dcl01/leek/data/tcga/v1', pattern = 'batch',
    full.names = TRUE), 'cross_sample_results', 'counts.tsv.gz')
names(counts_files) <- dir('/dcl01/leek/data/tcga/v1', pattern = 'batch')

## Read in counts info
counts <- lapply(counts_files, read.table, header = TRUE, sep = '\t',
    stringsAsFactors = FALSE, fill = TRUE)
counts <- do.call(rbind, counts)
counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads,
    ','), '[[', 1))
## Fix the one weird sample that did not download
counts$total.reads[counts$total.reads == ''] <- '0,0'
counts$totalReads <- as.integer(sapply(strsplit(counts$total.reads,
    ','), '[[', 1))
## Re-order accordingly
counts <- counts[match(metadata$gdc_file_id, tolower(counts$X)), ]
counts$bigwig_path <- paste0(sapply(strsplit(rownames(counts), '\\.'),
    function(x) {
        file.path('/dcl01/leek/data/tcga/v1', x[1], 'coverage_bigwigs/')
    }), counts$X, '.bw')
## Only one bigwig doesn't exist (again, that weird sample)
file.exists(counts$bigwig_path[counts$totalReads == 0])
counts$bigwig_path[counts$totalReads == 0] <- NA


## Create the reads metadata table for TCGA
readmeta <- matrix(NA, ncol = ncol(colData(rse_gene_SRP009615)), nrow = nrow(metadata))
readmeta <- as.data.frame(readmeta)
colnames(readmeta) <- colnames(colData(rse_gene_SRP009615))

readmeta$project <- 'TCGA'
readmeta$paired_end <- TRUE
readmeta$reads_downloaded <- counts$totalReads
readmeta$mapped_read_count <- counts$totalMapped
readmeta$bigwig_path <- counts$bigwig_path
readmeta$bigwig_file <- gsub('.*coverage_bigwigs/', '', readmeta$bigwig_path)


## Calculate AUCs
calculate_auc <- function(bw) {
    ## Choose name for temporary file
    auc_file <- file.path(tempdir(), paste0(gsub('.bw', '', basename(bw)),
        '.auc'))
    system(paste('wiggletools', 'AUC', auc_file, bw))
    res <- as.numeric(readLines(auc_file))
    ## Clean up
    unlink(auc_file)
    
    return(res)
}

readmeta$auc[!is.na(readmeta$bigwig_path)] <- unlist(
    bplapply(readmeta$bigwig_path[!is.na(readmeta$bigwig_path)], calculate_auc,
    BPPARAM = MulticoreParam(workers = 10,
        outfile = Sys.getenv('SGE_STDERR_PATH'))))

## Final merge
metadata <- cbind(readmeta, metadata)

## Remove the weird sample
metadata <- metadata[!is.na(metadata$bigwig_path), ]
rownames(metadata) <- NULL

print('Final dimensions of the metadata')
dim(metadata)

## Save final version
save(metadata, file = 'metadata.Rdata')

## Reproducibility info
Sys.time()
options(width = 120)
session_info()
