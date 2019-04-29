library('tools')
library('sessioninfo')

## Locate files
upload_files <- dir('/dcl01/lmarchio/data/FANTOM6/CAT_djn_analysis/RSE', full.names = TRUE)
names(upload_files) <- dir('/dcl01/lmarchio/data/FANTOM6/CAT_djn_analysis/RSE')

## Compute md5sums
print('Time spent calculating md5sums')
system.time( files_md5sum <- md5sum(upload_files) )

## Calculate total file size
print('Time spent calculating file sizes')
fileSize <- function(files) {
    as.numeric(system(paste('du -l', paste(files, collapse = ' '), 
        '| cut -f1'), intern = TRUE))
}

system.time( file_size <- sapply(upload_files, fileSize) )
project_size <- sum(file_size / 1024)
if(project_size < 1024) {
    print(paste('The total file size for fc_rc is', round(project_size, digits = 3), 'Mb.'))
} else {
    print(paste('The total file size for fc_rc is', round(project_size / 1024, digits = 3), 'Gb.'))
}

## Format consistent with the rest of the recount files_info.tsv files
write.table(data.frame(file = names(upload_files), md5sum = files_md5sum,
    size = file_size, stringsAsFactors = FALSE),
    file = 'fc_rc_files_info.tsv',
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)  

## For adding to recount_url
fc_rc_url_table <- data.frame(
    path = upload_files,
    file_name = names(upload_files),
    project = gsub('_.*', '', gsub('rse_fc_|\\.Rdata', '', names(upload_files))),
    version1 = NA,
    version2 = NA,
    url = paste0('http://idies.jhu.edu/recount/data/fc_rc/', names(upload_files)),
    stringsAsFactors = FALSE)
rownames(fc_rc_url_table) <- NULL

## save for later
save(fc_rc_url_table, file = 'fc_rc_url_table.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
