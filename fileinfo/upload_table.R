## Usage
# mkdir -p logs
# module load R/3.3
# Rscript upload_table.R > logs/upload_table_log.txt 2>&1

## Identify files to upload
upload <- c(
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_sra',
        full.names = TRUE),
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_gtex',
        full.names = TRUE)
)
upload <- file.path(upload, 'upload_files.Rdata')
names(upload) <-  c(
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_sra'),
    dir('/dcl01/leek/data/recount-website/fileinfo/fileinfo_gtex')
)

## Find all the info to upload
upload_table <- mapply(function(rda, project) {
    load(rda)
    data.frame(path = upload_files, file_name = names(upload_files),
        project = project, stringsAsFactors = FALSE)
}, upload, names(upload), SIMPLIFY = FALSE)
upload_table <- do.call(rbind, upload_table)
rownames(upload_table) <- NULL

## Remove GTEx bigWigs for now
print('Removing GTEx bigWigs')
dim(upload_table)
gtex <- upload_table$project == 'SRP012682' & grepl('bw',
    upload_table$file_name)
table(gtex)
upload_table <- upload_table[!gtex, ]
rownames(upload_table) <- NULL

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
