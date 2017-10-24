## Usage
# mkdir -p logs
# module load R/3.3.x
# Rscript upload_table.R > logs/upload_table_log.txt 2>&1

## Identify files to upload
upload <- c(
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_sra',
        full.names = TRUE),
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_gtex',
        full.names = TRUE),
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_tcga',
        full.names = TRUE)
)
upload <- file.path(upload, 'upload_files.Rdata')
names(upload) <-  c(
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_sra'),
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_gtex'),
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_tcga')
)

## Find all the info to upload
upload_table <- mapply(function(rda, project) {
    load(rda)
    data.frame(path = upload_files, file_name = names(upload_files),
        project = project, stringsAsFactors = FALSE)
}, upload, names(upload), SIMPLIFY = FALSE)
upload_table <- do.call(rbind, upload_table)

## Add files that have all of SRA
# upload_table <- rbind(upload_table,
#     data.frame(
#         path = c(
#             '/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_exon.Rdata',
#             '/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_gene.Rdata'
#         ),
#         file_name = c('rse_exon.Rdata', 'rse_gene.Rdata'),
#         project = 'sra', stringsAsFactors = FALSE
#     )
# )
## rse_exon.Rdata failed in the last run (it's too large...)
upload_table <- rbind(upload_table,
    data.frame(
        path = c(
            '/dcl01/leek/data/recount-website/rse/rse_sra/all/rse_gene.Rdata'
        ),
        file_name = c('rse_gene.Rdata'),
        project = 'sra', stringsAsFactors = FALSE
    )
)

rownames(upload_table) <- NULL

## Internal table with all paths for recount.bwtool
#local_url <- upload_table
#is.bw <- grepl('[.]bw$', local_url$file_name)
#local_url$url <- NA
#local_url$url[!is.bw] <- paste0('http://duffel.rail.bio/recount/',
#    local_url$project[!is.bw], '/', local_url$file_name[!is.bw])
#local_url$url[is.bw] <- paste0('http://duffel.rail.bio/recount/',
#    local_url$project[is.bw], '/bw/', local_url$file_name[is.bw])
#save(local_url, file = 'local_url.RData', compress = 'xz')

## Remove GTEx bigWigs for now
#print('Removing GTEx bigWigs')
#dim(upload_table)
#gtex <- upload_table$project == 'SRP012682' & grepl('bw',
#    upload_table$file_name)
#table(gtex)
#upload_table <- upload_table[!gtex, ]
#rownames(upload_table) <- NULL

## Explore table a little bit
print('Final table: minor exploration')
head(upload_table)
dim(upload_table)

## Save info
save(upload_table, file = 'upload_table.Rdata')

write.table(upload_table, file = 'upload_table.tsv', sep = '\t',
    row.names = FALSE, quote = FALSE, col.names = FALSE)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
