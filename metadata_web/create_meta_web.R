## Prepare metadata
# module load conda_R/3.4.x
# mkdir -p logs
# Rscript create_meta_web.R -p "sra" > logs/create_meta_web_sra_log.txt 2>&1
# Rscript create_meta_web.R -p "gtex" > logs/create_meta_web_gtex_log.txt 2>&1
# Rscript create_meta_web.R -p "tcga" > logs/create_meta_web_tcga_log.txt 2>&1

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
    transcripts = NA,
    phenotype = NA,
    files_info = NA,
    stringsAsFactors = FALSE
)
rownames(meta_web) <- NULL

if(opt$project == 'tcga') meta_web$accession[1] <- '<a href="https://cancergenome.nih.gov/" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'srainfo\', 1)">TCGA</a>'

for(project in projects) {
    if(project != 'TCGA') {
        jx_files <- c(
            file.path('/dcl01/leek/data/recount-website/rse',
                paste0('rse_', opt$project), project, 'rse_jx.Rdata'),
            file.path('/dcl01/leek/data/recount_junctions_2',
                paste0(project, '.junction_id_with_transcripts.bed.gz')),
            file.path('/dcl01/leek/data/recount_junctions_2',
                paste0(project, '.junction_coverage.tsv.gz')),
            file.path('/dcl01/leek/data/recount-website/rse',
                paste0('rse_', opt$project), project, 'counts_jx.tsv.gz')
        )
    } else {
        jx_files <- c(
            file.path('/dcl01/leek/data/recount-website/rse',
                paste0('rse_', opt$project), project, 'rse_jx.Rdata'),
        '/dcl01/leek/data/tcga_work/tcga_recount_junctions/TCGA.junction_id_with_transcripts.bed.gz',
            '/dcl01/leek/data/tcga_work/tcga_recount_junctions/TCGA.junction_coverage.tsv.gz',
            file.path('/dcl01/leek/data/recount-website/rse',
                paste0('rse_', opt$project), project, 'counts_jx.tsv.gz')
        )
    }
    
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
    tx_files <- paste0('/dcl01/leek/data/ta_poc/recount_out/rse_new/',
        project, '/rse_tx.RData')
    
    meta_web$gene[projects == project] <- paste(c(
        paste0(
            '<a href="http://duffel.rail.bio/recount/v2/', project,
            '/rse_gene.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-gene\', 1)">RSE v2</a>'
        ),
        paste0(' <a href="http://duffel.rail.bio/recount/v2/', project,
        '/counts_gene.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-counts-gene\', 1)">counts v2</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_gene.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-gene\', 1)">RSE v1</a>'
        ),
        paste0(' <a href="http://duffel.rail.bio/recount/', project,
        '/counts_gene.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-counts-gene\', 1)">counts v1</a>')
        )[rep(file.exists(gene_files), each = 2)], collapse = ' ')
    meta_web$exon[projects == project] <- paste(c(
        paste0(
            '<a href="http://duffel.rail.bio/recount/v2/', project,
            '/rse_exon.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-exon\', 1)">RSE v2</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/v2/', project,
        '/counts_exon.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-counts-exon\', 1)">counts v2</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_exon.Rdata" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-exon\', 1)">RSE v1</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
        '/counts_exon.tsv.gz" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-counts-exon\', 1)">counts v1</a>')
        )[rep(file.exists(exon_files), each = 2)], collapse = ' ')    
    meta_web$transcripts[projects == project] <- paste(c(
        paste0(
            '<a href="http://duffel.rail.bio/recount/v2/', project,
            '/rse_tx.RData" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-tx\', 1)">RSE v2</a>'
        ),
        paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_tx.RData" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-tx\', 1)">RSE v1</a>'
        ))[rep(file.exists(tx_files), each = 2)], collapse = ' ')
    
    if(opt$project != 'sra') {      
        ## Add tissue files
        rse_dir <- ifelse(opt$project == 'gtex',
            '/dcl01/leek/data/recount-website/rse/rse_gtex/SRP012682',
            '/dcl01/leek/data/recount-website/rse/rse_tcga/TCGA')
        extra_exon <- dir(rse_dir, 'rse_exon_')
        names(extra_exon) <- gsub('_', ' ', gsub('rse_exon_|.Rdata', '',
            extra_exon))
        extra_gene <- dir(rse_dir, 'rse_gene_')
        names(extra_gene) <- gsub('_', ' ', gsub('rse_gene_|.Rdata', '',
            extra_gene))
        rse_tx_dir <- ifelse(opt$project == 'gtex',
            '/dcl01/leek/data/ta_poc/recount_out/rse_new/SRP012682',
            '/dcl01/leek/data/ta_poc/recount_out/rse_new/TCGA')
        extra_tx <- dir(rse_tx_dir, 'rse_tx_.*R[d|D]ata')
        names(extra_tx) <- gsub('_', ' ', gsub('rse_tx_|.R[d|D]ata', '',
            extra_tx))
        
        meta_web$gene[projects == project] <- paste(
            meta_web$gene[projects == project], 'RSE by tissue (version 2):',
            paste0(
                '<a href="http://duffel.rail.bio/recount/v2/', project,
                '/', extra_gene, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-gene\', 1)">', names(extra_gene),'</a>'
            , collapse = ' '), 'RSE by tissue (version 1):',
            paste0(
                '<a href="http://duffel.rail.bio/recount/', project,
                '/', extra_gene, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-gene\', 1)">', names(extra_gene),'</a>'
            , collapse = ' '), collapse = ' ')
        
        meta_web$exon[projects == project] <- paste(
            meta_web$exon[projects == project], 'RSE by tissue (version 2):',
            paste0(
                '<a href="http://duffel.rail.bio/recount/v2/', project,
                '/', extra_exon, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-exon\', 1)">', names(extra_exon),'</a>'
            , collapse = ' '), 'RSE by tissue (version 1):',
            paste0(
                '<a href="http://duffel.rail.bio/recount/', project,
                '/', extra_exon, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-exon\', 1)">', names(extra_exon),'</a>'
            , collapse = ' '), collapse = ' ')
            
        meta_web$transcripts[projects == project] <- paste(
            meta_web$tx[projects == project], 'RSE by tissue (version 2):',
            paste0(
                '<a href="http://duffel.rail.bio/recount/v2/', project,
                '/', extra_tx, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-tx\', 1)">', names(extra_tx),'</a>'
            , collapse = ' '), 'RSE by tissue (version 1):',
            paste0(
                '<a href="http://duffel.rail.bio/recount/', project,
                '/', extra_tx, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-tx\', 1)">', names(extra_tx),'</a>'
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
        meta_web$files_info[projects == project] <- paste(c(
            paste0(
            '<a href="http://duffel.rail.bio/recount/v2/', project,
            '/files_info.tsv" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-filesinfo\', 1)">v2</a>'),
            paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/files_info.tsv" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-filesinfo\', 1)">v1</a>')
        ), collapse = ' ')
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
