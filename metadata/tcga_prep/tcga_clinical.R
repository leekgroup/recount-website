## Usage
# module load R/devel
# mkdir -p logs
# Rscript tcga_clinical.R > logs/tcga_clinical_log.txt 2>&1

## Based on https://support.bioconductor.org/p/89315/#89375

library('TCGAbiolinks')
library('plyr')
library('devtools')
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]

clin <- lapply(projects, function(p) {
    message(paste(Sys.time(), 'processing project', p))
    result <- tryCatch({
        query <- GDCquery(project = p, data.category = 'Clinical')
        GDCdownload(query)
        GDCprepare_clinic(query, clinical.info = 'patient')
    }, error = function(e) {
        message(paste0('Error clinical: ', p))
        return(NULL)
    })
    return(result)
})
names(clin) <- projects

## Merge all
clin_all <- rbind.fill(clin)

## Fix columns that have '' that should be NAs
for(j in seq_len(ncol(clin_all))) {
    i <- which(clin_all[, j] == '')
    if(length(i) > 0) clin_all[i, j] <- NA
}

save(clin_all, file = 'clin_all.Rdata')

write.table(clin_all, file = 'clin_all.tsv', quote = FALSE, row.names = FALSE,
    sep = '\t')

## Reproducibility info
Sys.time()
options(width = 120)
session_info()
