## Merge the RSE exon and gene objects and create the jx RSE object

## For specifying parameters
library('getopt')

## Specify parameters
spec <- matrix(c(
    'bigwig_path', 'b', 1, 'character',
    'Path to the directory with the bigwig files',
	'jx_file', 'j', 1, 'character',
    'Path to the first_pass_junctions.tsv.gz file from the cross-sample results',
    'manifest_file', 'm', 1, 'character',
    'Path to the manifest file used',
    'wiggletools', 'w', 2, 'character',
    "Path to wiggletools. If not provided, it's assumed that it is on the $PATH",
    'wigToBigWig', 't', 2, 'character',
    "Path to wigToBigWig. If not provided, it's assumed that it is on the $PATH",
    'calculate_mean', 'c', 2, 'logical', 'Whether to calculate the AUC',
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
        'bigwig_path' = '/dcl01/leek/data/sunghee_analysis/processed/coverage_bigwigs',
        'jx_file' = '/dcl01/leek/data/sunghee_analysis/processed/cross_sample_results/first_pass_junctions.tsv.gz',
        'manifest_file' = '/dcl01/leek/data/sunghee/all_s3.manifest',
        'wiggletools' = 'wiggletools',
        'wigToBigWig' = 'wigToBigWig',
        'calculate_mean' = TRUE
    )
}

## Check inputs
stopifnot(dir.exists(opt$bigwig_path))
stopifnot(file.exists(opt$jx_file))
stopifnot(file.exists(opt$manifest_file))
stopifnot(file.exists('introns_unique.Rdata'))
stopifnot(file.exists('hg38.sizes'))

## Check that outputs don't exist, to avoid overwriting
stopifnot(!file.exists('rse_gene.Rdata'))
stopifnot(!file.exists('rse_exon.Rdata'))
stopifnot(!file.exists('rse_jx.Rdata'))

## Are we on JHPCE? Print some helpful info
jhpce <- grepl('compute-', Sys.info()['nodename'])
if(jhpce & is.null(opt$wiggletools)) {
    message(paste(Sys.time(), 'Note that you can use wiggletools with:
    module load wiggletools/default
    '))
}
if(jhpce & is.null(opt$wigToBigWig)) {
    message(paste(Sys.time(), 'Note that you can use wigToBigWig with:
    module load ucsctools
    '))
}


## Set some defaults
if(is.null(opt$calculate_mean)) opt$calculate_mean <- TRUE
if(is.null(opt$wiggletools)) opt$wiggletools <- 'wiggletools'
if(is.null(opt$wigToBigWig)) opt$wigToBigWig <- 'wigToBigWig'
    
## Print options used
message(paste(Sys.time(), 'options used:'))
print(opt)

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

## Read manifest info
message(paste(Sys.time(), 'reading', opt$manifest_file))
manifest <- read.table(opt$manifest_file, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE, fill = TRUE)

## Get sample names from the manifest file. Note that a manifest file can
## have both paired-end and single-end data
manifest_cols <- apply(manifest, 1, function(x) { sum(x != '') })
manifest_samples <- sapply(seq_len(nrow(manifest)), function(i) {
    manifest[i, manifest_cols[i]]
})
manifest$paired <- ifelse(manifest_cols == 5, TRUE, FALSE)

## Locate exon rse objects, load them, merge them and save results
exon_files <- dir('rse_temp', 'rse_exon_', full.names = TRUE)
rse_exon <- do.call(cbind, lapply(exon_files, load_rse))

## Assign paired-end info
colData(rse_exon)$paired <- manifest$paired[match(rownames(colData(rse_exon)),
    manifest_samples)]

message(paste(Sys.time(), 'saving rse_exon.Rdata'))
save(rse_exon, file = 'rse_exon.Rdata')
rm(rse_exon, exon_files)

## Same for gene rse objects
gene_files <- dir('rse_temp', 'rse_gene_', full.names = TRUE)
rse_gene <- do.call(cbind, lapply(gene_files, load_rse, type = 'gene'))
colData(rse_gene)$paired <- manifest$paired[match(rownames(colData(rse_gene)),
    manifest_samples)]
message(paste(Sys.time(), 'saving rse_gene.Rdata'))
save(rse_gene, file = 'rse_gene.Rdata')

## Get metadata information
metadata <- colData(rse_gene)
rm(rse_gene, gene_files)

## Calculate the mean bigwig if necessary
if(opt$calculate_mean) {
    ## Check that outputs don't exist
    stopifnot(!file.exists('bw/mean.bw'))
    stopifnot(!file.exists('bw/mean.wig'))
    
    dir.create('bw', showWarnings = FALSE)
    
    ## Name resulting mean.bw file
    outbw <- 'bw/mean.bw'
    outwig <- 'bw/mean.wig'
    

    scaleWig <- function(m) {
        print(m)
        paste(paste('scale', round(1e6*100*40 / m$auc, digits = 17),
            file.path(opt$bigwig_path, m$bigwig_file)), collapse = ' ')
    }
    runCmd <- function(cmd, i = NULL) {
        if(is.null(i)) {
            shell_name <- '.createWig.sh'
        } else {
            shell_name <- paste0('.createWig_part', i, '.sh')
        }
        message(paste(Sys.time(), 'command used:', cmd))
        cat(cmd, file = shell_name)
        system(paste('sh', shell_name))
    }

    ## Calculate mean bigwig
    if(nrow(metadata) < 100) {
        ## Scale commands
        cmd <- scaleWig(metadata)
        ## Calculate mean wig file
        message(paste(Sys.time(), 'creating file', outwig))
        cmd <- paste(opt$wiggletools, 'write', outwig, 'mean', cmd)
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
            cmd <- paste(opt$wiggletools, 'write', tmpwig, 'sum', cmd)
            runCmd(cmd, i)
            return(tmpwig)
        }, meta, names(meta)) )
    
        ## Calculate final mean
        cmd <- paste(opt$wiggletools, 'write', outwig, 'scale',
            1/nrow(metadata), 'sum', paste(tmpfiles, collapse = ' '))
        system.time( runCmd(cmd) )
        
        ## Clean up
        sapply(tmpfiles, unlink)
    }

    ## Transform to bigwig file
    message(paste(Sys.time(), 'creating file', outbw))
    cmd2 <- paste(opt$wigToBigWig, outwig, 'hg38.sizes', outbw)
    system.time( system(cmd2) )
}


## Code for creating rse_jx
message(paste(Sys.time(), 'reading', opt$jx_file))
jx_info <- read.table(opt$jx_file, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE, check.names = FALSE)
colnames(jx_info) <- c('chr', 'start', 'end', 'sample_ids', 'reads')
    
## Create the counts matrix
message(paste(Sys.time(), 'processing count information'))
jx_info_samples <- strsplit(jx_info$sample_ids, ',')
jx_info_reads <- strsplit(jx_info$reads, ',')
stopifnot(identical(elementNROWS(jx_info_samples), elementNROWS(jx_info_reads)))
jx_info_tab <- data.frame(
    jx_id = rep(seq_len(nrow(jx_info)), elementNROWS(jx_info_samples)),
    sample_id = unlist(jx_info_samples),
    reads = as.numeric(unlist(jx_info_reads)), stringsAsFactors = FALSE
)
rm(jx_info_samples, jx_info_reads)

message(paste(Sys.time(), 'creating junction counts table'))
## Create junction counts table
jx_counts <- matrix(0, ncol = nrow(metadata), nrow = nrow(jx_info))
colnames(jx_counts) <- rownames(metadata)

## Fill in table
for(sample in rownames(metadata)) {
    sampleId <- as.character(which(manifest_samples == sample))
    sample_reads <- subset(jx_info_tab, sample_id == sampleId)
    if(nrow(sample_reads) == 0)  {
        message(paste(Sys.time(), 'found no junction counts for sample',
            sample))
        next
    }
    jx_map <- match(jx_info$jx_id, sample_reads$jx_id)
    jx_counts[!is.na(jx_map), sample] <- sample_reads$reads[jx_map[!is.na(jx_map)]]
}
rm(sample, sampleId, sample_reads, jx_map, jx_info_tab)

## Create a GRanges object for the exon-exon junctions
message(paste(Sys.time(),
    'forming GRanges object with exon-exon junctions information'))
jx_gr <- GRanges(seqnames = gsub('\\+|-', '', jx_info$chr), IRanges(
    start = jx_info$start, end = jx_info$end), strand = ifelse(grepl('\\+',
        jx_info$chr), '+', ifelse(grepl('-', jx_info$chr), '-', '*')))


chr_info <- read.table('hg38.sizes', sep = '\t',
    col.names = c('chr', 'length'), stringsAsFactors = FALSE)
chrs <- chr_info$length
names(chrs) <- chr_info$chr
seqlengths(jx_gr) <- chrs[names(seqlengths(jx_gr))]

## Add recount columns that are NA here
jx_gr$junction_id <- as.character(NA)
jx_gr$found_junction_gencode_v24 <- rep(CharacterList(NA), length(jx_gr))


message(paste(Sys.time(),
    'finding tx_name and gene_id based on the intron reference set'))
## Now to actual data, add the transcript names and gene ids
load('introns_unique.Rdata')
oo <- findOverlaps(jx_gr, introns_unique, type = 'equal')
stopifnot(length(unique(queryHits(oo))) == length(oo))

jx_gr$symbol <- jx_gr$gene_id <- jx_gr$tx_name <- jx_gr$gene_id_proposed  <-  jx_gr$symbol_proposed <- CharacterList(NA)
left_gene <- left_symbol <- right_gene <- right_symbol <- jx_gr$gene_id_proposed


## Partial overlap
message(paste(Sys.time(),
    'finding gene ids and symbols based on partial matching'))
both <- countOverlaps(jx_gr, introns_unique, type = 'equal') > 0
not_both <- which(!both)

## Left
oo_left <- findOverlaps(jx_gr[not_both], introns_unique, type = 'start')
left_gene[not_both[queryHits(oo_left)]] <- introns_unique$gene_id[subjectHits(oo_left)]
left_symbol[not_both[queryHits(oo_left)]] <- introns_unique$symbol[subjectHits(oo_left)]

## Right
oo_right <- findOverlaps(jx_gr[not_both], introns_unique, type = 'end')
right_gene[not_both[queryHits(oo_right)]] <- introns_unique$gene_id[subjectHits(oo_right)]
right_symbol[not_both[queryHits(oo_right)]] <- introns_unique$symbol[subjectHits(oo_right)]

## Manually combine, since using paste(x, y, sep = '-') won't work
## for cases where x and y have different lengths
manual_c <- function(l, r) {
    tmp <- merge(l, r, all = TRUE)
    tmp2 <- tmp[!is.na(tmp$value), ]
    ## Add back if only it's NAs
    missed <- unique(tmp$group)[!unique(tmp$group) %in% tmp2$group]
    if(length(missed) > 0) {
        tmp2 <- rbind(tmp2, subset(tmp, group %in% missed))
    }
    
    res <- CharacterList(split(as.character(tmp2$value), tmp2$group))
    names(res) <- NULL
    return(res)
}

message(paste(Sys.time(), 'combining left and right results'))
has_hit <- not_both[unique(c(queryHits(oo_left), queryHits(oo_right)))]
ends_hit <- not_both[intersect(queryHits(oo_left), queryHits(oo_right))]
jx_gr$gene_id_proposed[has_hit] <- manual_c(left_gene[has_hit],
    right_gene[has_hit])
jx_gr$symbol_proposed[has_hit] <- manual_c(left_symbol[has_hit],
    right_symbol[has_hit])

## Full match
jx_gr$gene_id_proposed[queryHits(oo)] <- jx_gr$gene_id[queryHits(oo)] <- introns_unique$gene_id[subjectHits(oo)]
jx_gr$symbol_proposed[queryHits(oo)] <- jx_gr$symbol[queryHits(oo)] <- introns_unique$symbol[subjectHits(oo)]
jx_gr$tx_name[queryHits(oo)] <- introns_unique$tx_name[subjectHits(oo)]

## See how junctions matched
print('Exon-exon junctions matching or partial matching')
print(addmargins(table('has proposed gene_id' = !any(is.na(jx_gr$gene_id_proposed)), 'exact match, gene_id' = !any(is.na(jx_gr$gene_id)))))
 
## Assign class
message(paste(Sys.time(), 'assigning class'))
left <- countOverlaps(jx_gr, introns_unique, type = 'start') > 0
right <- countOverlaps(jx_gr, introns_unique, type = 'end') > 0
jx_gr$class <- ifelse(both, 'annotated',
    ifelse(left & right, 'exon_skip',
    ifelse(left | right, 'alternative_end', 'novel')))
    
## Detect fusion
message(paste(Sys.time(), 'detecting gene fusions'))

find_fusions <- function(fu) {
    tmp <- merge(left_gene[fu], right_gene[fu])
    seq_len(length(fu))[!seq_len(length(fu)) %in% tmp$group]
}    
jx_gr$class[ends_hit[find_fusions(ends_hit)]] <- 'fusion'

print('Exon-exon junctions by class')
print(table(jx_gr$class))


## Create the junctions level rse
message(paste(Sys.time(), 'creating rse_jx object'))
rse_jx <- SummarizedExperiment(assays = list('counts' = jx_counts),
    colData = metadata, rowRanges = jx_gr)
message(paste(Sys.time(), 'saving rse_jx.Rdata'))
save(rse_jx, file = 'rse_gene.Rdata')


## Clean up
message(paste(Sys.time(), 'cleaning up temporary files'))
to_clean <- c(outwig, 'ucsc-knowngene-hg38.bed',
    'ucsc-knowngene-hg38-genes-bp-length.Rdata',
    'ucsc-knowngene-hg38-exons.Rdata', 'count_groups.Rdata', 'hg38.sizes',
    'introns_unique.Rdata')
sapply(to_clean, unlink)
unlink('rse_temp', recursive = TRUE)

## Reproducibility info
proc.time()
options(width = 120)
session_info()
