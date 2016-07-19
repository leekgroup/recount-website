## Prepare metadata
# module load R/3.3
# mkdir -p logs
# Rscript create_meta_web.R -p "sra" > logs/create_meta_web_sra_log.txt 2>&1
# Rscript create_meta_web.R -p "gtex" > logs/create_meta_web_gtex_log.txt 2>&1

library('getopt')
library('GenomicRanges')

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

## For testing
if(FALSE) {
    opt <- list(project = 'sra')
}

## Load metadata
file <- file.path('/dcl01/leek/data/recount-website/metadata/',
    paste0('metadata_', opt$project, '.Rdata'))
stopifnot(file.exists(file))
load(file)

## Load abstracts
file <- file.path('/dcl01/leek/data/recount-website/metadata_abstract/', 
    paste0('abstracts_', opt$project, '.Rdata'))
stopifnot(file.exists(file))
load(file)

projects <- unique(metadata$project)
meta_web <- data.frame(
    accession = paste0(
        '<a href="http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=',
        projects, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'srainfo\', 1)">', projects, '</a>'),
    number_samples = sapply(projects, function(project) {
        sum(metadata$project == project)}),
    species = 'human',
    abstract = abstracts$study_abstract[match(projects,
        abstracts$study_accession)],
    gene = NA,
    exon = NA,
    junctions = NA,
    phenotype = NA,
    genes = '<a href="https://jhubiostatistics.shinyapps.io/recount/ucsc-knowngene-hg38-genes-bp-length.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'genes\', 1)">link</a>',
    exons = '<a href="https://jhubiostatistics.shinyapps.io/recount/ucsc-knowngene-hg38-exons.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'exons\', 1)">link</a>',
    files_info = NA,
    stringsAsFactors = FALSE
)
rownames(meta_web) <- NULL
for(project in projects) {
    jx_files <- c(
        file.path('/dcl01/leek/data/recount-website/rse',
            paste0('rse_', opt$project), project, 'rse_jx.Rdata'),
        file.path('/dcl01/leek/data/recount_junctions',
            paste0(project, '.junction_id_with_transcripts.bed.gz')),
        file.path('/dcl01/leek/data/recount_junctions',
            paste0(project, '.junction_coverage.tsv.gz')),
        file.path('/dcl01/leek/data/recount-website/rse',
            paste0('rse_', opt$project), project, 'counts_jx.tsv.gz')
    )
    exon_files <- c(
        file.path('/dcl01/leek/data/recount-website/rse/',
            paste0('rse_', opt$project), project, 'rse_exon.Rdata'),
        file.path('/dcl01/leek/data/recount-website/rse/',
            paste0('rse_', opt$project), project, 'counts_exon.tsv.gz')
    )
    gene_files <- c(
        file.path('/dcl01/leek/data/recount-website/rse/',
            paste0('rse_', opt$project), project, 'rse_gene.Rdata'),
        file.path('/dcl01/leek/data/recount-website/rse/',
            paste0('rse_', opt$project), project, 'counts_gene.tsv.gz')
    )
    
    meta_web$gene[projects == project] <- paste(c(
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_gene.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-gene\', 1)">RSE</a>'
        ), 
        paste0(' <a href="http://duffel.rail.bio/recount/', project,
        '/counts_gene.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-counts-gene\', 1)">counts</a>')
        )[file.exists(gene_files)], collapse = ' ')
    meta_web$exon[projects == project] <- paste(c(
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_exon.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-exon\', 1)">RSE</a>'
        ), 
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
        '/counts_exon.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-counts-exon\', 1)">counts</a>')
        )[file.exists(exon_files)], collapse = ' ')
    if(opt$project == 'gtex') {      
        ## Add tissue files
        extra_exon <- dir(
            '/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682',
            'rse_exon_')
        names(extra_exon) <- gsub('_', ' ', gsub('rse_exon_|.Rdata', '',
            extra_exon))
        extra_gene <- dir(
            '/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682',
            'rse_gene_')
        names(extra_gene) <- gsub('_', ' ', gsub('rse_gene_|.Rdata', '',
            extra_gene))
        
        meta_web$gene[projects == project] <- paste(
            meta_web$gene[projects == project], 'RSE by tissue:',
            paste0(
                '<a href="http://duffel.rail.bio/recount/', project,
                '/', extra_gene, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-gene\', 1)">', names(extra_gene),'</a>'
            , collapse = ' '), collapse = ' ')
        
        meta_web$exon[projects == project] <- paste(
            meta_web$exon[projects == project], 'RSE by tissue:',
            paste0(
                '<a href="http://duffel.rail.bio/recount/', project,
                '/', extra_exon, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-exon\', 1)">', names(extra_exon),'</a>'
            , collapse = ' '), collapse = ' ')
    }
    
    meta_web$junctions[projects == project] <- paste(c(
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_jx.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-jx\', 1)">RSE</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project, '/',
            project, '.junction_id_with_transcripts.bed.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-jx-bed\', 1)">jx_bed</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project, '/',
        project, '.junction_coverage.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-jx-cov\', 1)">jx_cov</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
        '/counts_jx.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-counts-jx\', 1)">counts</a>'
        ))[file.exists(jx_files)], collapse = ' ')
    if(file.exists(file.path('/dcl01/leek/data/recount-website/metadata/', paste0('project_metadata_', opt$project), paste0(project, '.tsv')))) {   
        meta_web$phenotype[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project, '/', project,
            '.tsv" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-phenotype\', 1)">link</a>')
    }
    if(file.exists(file.path('/dcl01/leek/data/recount-website/fileinfo/', paste0('fileinfo_', opt$project), project, 'files_info.tsv'))) {
        meta_web$files_info[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/files_info.tsv" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-filesinfo\', 1)">link</a>')
    }
}

## Explore number of samples per project
message(paste(Sys.time(), "number of samples per project"))
summary(meta_web$number_samples)

## How much of it is missing?
message(paste(Sys.time(), "number missing"))
sapply(meta_web, function(x) sum(is.na(x) | x == ''))

## Save info
save(meta_web, file = paste0('meta_web_', opt$project, '.Rdata'))
write.table(meta_web, file = paste0('meta_web_', opt$project, '.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
