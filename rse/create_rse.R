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
    opt <- list(project = 'sra', 'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata',
            projectid = 'ERP005274')
    opt <- list(project = 'sra', 'metadata' = '/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata',
            projectid = 'SRP002915')
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

## Remove bigwig and tsv file paths
metadata_clean <- metadata[, !colnames(metadata) %in% c('bigwig_path',
    'tsv_path')]

## Create output dir for the project
outdir <- paste0('rse_', opt$project, '/', opt$projectid)
dir.create(outdir, showWarnings = FALSE)



#### Junction level ####

message(paste(Sys.time(), 'loading junction coverage file'))
## Load project junctions info
jx_file <- file.path('/dcl01/leek/data/recount_junctions',
    paste0(opt$projectid, '.junction_coverage.tsv.gz'))
jx_file_bed <- file.path('/dcl01/leek/data/recount_junctions',
    paste0(opt$projectid, '.junction_id_with_transcripts.bed.gz'))
hasJx <- file.exists(jx_file)
if(!hasJx) message(paste('Missing file', jx_file))

if(!file.exists(jx_file_bed)) {
    message(paste('Missing file', jx_file_bed))
    hasJx <- FALSE
}

if(hasJx) {
    ## Load junctions sample information
    message(paste(Sys.time(), 'loading junctions sample information'))
    jx_samples <- read.table(
        '/dcl01/leek/data/recount_junctions/sample_ids.tsv',
        sep = '\t', col.names = c('sample_id', 'project', 'run'),
        stringsAsFactors = FALSE, colClasses = 'character')
        
        jx_project <- read.table(jx_file, sep = '\t',
            col.names = c('jx_id', 'sample_ids', 'reads'),
            stringsAsFactors = FALSE, colClasses = 'character')
        print('jx_project dimensions')
        dim(jx_project)

        message(paste(Sys.time(), 'creating jx_project_tab object'))


    if(opt$project == 'gtex') {
        ## Create a table with 1 row per sample for a given junction
        if(!file.exists(file.path(outdir, 'jx_project_tab.Rdata'))) {
            jx_project.start <- seq(from = 1, to = nrow(jx_project), by = 1e5)
            jx_project.end <- c(
                jx_project.start[2:length(jx_project.start)] - 1,
                nrow(jx_project)
            )
        
            jx_project_tab <- mapply(function(start, end) {
                jx_split <- jx_project[start:end, ]
                jx_project_samples <- strsplit(jx_split$sample_ids, ',')
                jx_project_reads <- strsplit(jx_split$reads, ',')
                stopifnot(identical(elementNROWS(jx_project_samples),
                    elementNROWS(jx_project_reads)))
                res <- data.frame(
                    jx_id = rep(jx_split$jx_id,
                        elementNROWS(jx_project_samples)),
                    sample_id = unlist(jx_project_samples),
                    reads = as.numeric(unlist(jx_project_reads)),
                    stringsAsFactors = FALSE
                )
                return(res)
            }, jx_project.start, jx_project.end, SIMPLIFY = FALSE)
            rm(jx_project.start, jx_project.end)
            message(paste(Sys.time(), 'saving jx_project_tab'))
            save(jx_project_tab, file = file.path(outdir,
                'jx_project_tab.Rdata'))
            rm(jx_project_tab)
        } else {
            message(paste(Sys.time(), 'using previously created jx_project_tab'))
        }
        

        suppressPackageStartupMessages(library('Matrix'))
        suppressPackageStartupMessages(library('BiocParallel'))
        suppressPackageStartupMessages(library('BatchJobs'))
        suppressPackageStartupMessages(library('Hmisc'))
        
        ## Setup parallel envir
        funs <- makeClusterFunctionsSGE('bioc.tmpl')
        param <- BatchJobsParam(36, cluster.functions = funs, cleanup = FALSE)
        
        message(paste(Sys.time(), 'creating junction counts (list)'))
        
        ## Find sample ids
        sample_ids <- lapply(metadata_clean$run, function(run) {
            sample <- jx_samples$sample_id[jx_samples$run == run]
        })
        
        ## Number of junctions
        jx_n <- length(unique(jx_project$jx_id))
        print('Number of junctions in GTEx')
        print(jx_n)
        
        
        bioc_prep <- function(runs, samples, jx_n) {
            
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
            
            library('Matrix')
            
            message(paste(Sys.time(), 'loading jx_project_tab file'))
            load('/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682/jx_project_tab.Rdata')
            
            res_files <- mapply(bioc_run, runs, samples,
                MoreArgs = list(jx_n = jx_n), SIMPLIFY = FALSE)
            
            ## Load results
            message(paste(Sys.time(), 'loading sample results'))
            jx_count <- lapply(res_files, bioc_load)
            
            ## Create junction counts table for subset
            message(paste(Sys.time(), 'running cbind on junction counts'))
            jx_count <- do.call(cbind, jx_count)
            
            print('jx_count size and dimensions')
            print(object.size(jx_count), units = 'Mb')
            print(dim(jx_count))
            
            return(jx_count)
        }
        
        ## Create groups of samples to analyze at a time
        i.groups <- cut2(seq_len(length(metadata_clean$run)), g = 36)
        run.groups <- split(metadata_clean$run, i.groups)
        sample.groups <- split(sample_ids, i.groups)
        
        ## Extract information by sample
        jx_counts <- bpmapply(bioc_prep, run.groups, sample.groups,
            MoreArgs = list(jx_n = jx_n), BPPARAM = bpparam, SIMPLIFY = FALSE)
        rm(i.groups, run.groups, sample.groups, jx_n, sample_ids)

        ## Create junction counts table
        message(paste(Sys.time(), 'running cbind on junction counts subsets'))
        jx_counts <- do.call(cbind, jx_counts)

        message(paste(Sys.time(), 'saving junction counts'))
        save(jx_counts, file = file.path(outdir, 'jx_counts.Rdata'))
    } else {

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
        jx_counts <- matrix(0, ncol = nrow(metadata_clean),
            nrow = nrow(jx_project))
        colnames(jx_counts) <- metadata_clean$run

        ## Fill in table
        for(run in metadata_clean$run) {
            sample <- jx_samples$sample_id[jx_samples$run == run]
            sample_reads <- subset(jx_project_tab, sample_id == sample)
            if(nrow(sample_reads) == 0)  {
                message(paste(Sys.time(), 'found no junction counts for run',
                    run))
                next
            }
            jx_map <- match(jx_project$jx_id, sample_reads$jx_id)
            jx_counts[!is.na(jx_map), run] <- sample_reads$reads[jx_map[!is.na(jx_map)]]
        }
        rm(sample, sample_reads, jx_map, run, jx_project_tab)
    }

    print('Memory used by junction counts')
    print(object.size(jx_counts), units = 'Mb')

    ## Save junction counts
    if(opt$project != 'gtex') {
        message(paste(Sys.time(), 'writing file', file.path(outdir,
            'counts_jx.tsv')))
        write.table(as.data.frame(jx_counts), file = file.path(outdir,
            'counts_jx.tsv'), sep = '\t', row.names = FALSE, quote = FALSE,
            col.names = TRUE)
        system(paste('gzip', file.path(outdir, 'counts_jx.tsv')))
    }



    ## Load junction bed file
    message(paste(Sys.time(), 'loading the junctions bed file'))
    jx_bed <- import.bed(jx_file_bed)

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
        '/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes',
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
    hasIds <- any(!is.na(trans_names$name))
    if(hasIds) {
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

    }

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

    if(hasIds) {
        message(paste(Sys.time(), 'finding the gene ids for each transcript'))
        map_jx <- match(jx_bed$junction_id, trans_names$jx_id)
        system.time( jx_bed$gene_ids[!is.na(map_jx)] <- do.call(c,
            lapply(jx_bed$junction_id[!is.na(map_jx)], find_gene)) )
    }


    ## Create the junctions level rse
    rse_jx <- SummarizedExperiment(assays = list('counts' = jx_counts),
        colData = DataFrame(metadata_clean), rowRanges = jx_bed)
    message(paste(Sys.time(), 'writing file', file.path(outdir, 'rse_jx.Rdata')))
    save(rse_jx, file = file.path(outdir, 'rse_jx.Rdata'))

    rm(rse_jx, jx_bed, jx_counts, trans_names)
    if(hasIds) rm(map_jx, map_gene, unique_names)
} else {
    message('Skipping junctions since files are missing')
}





#### Gene and exon level ####
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
rm(counts, counts_gene, rse_exon, rse_gene)




## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
