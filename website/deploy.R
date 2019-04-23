## Can now be deployed with R 3.3.0

## Load metadata information
load('../metadata_web/meta_web_sra.Rdata')
if(file.exists('../metadata_web/meta_web_gtex.Rdata')) {
    m <- meta_web
    load('../metadata_web/meta_web_gtex.Rdata')
    meta_web <- rbind(m, meta_web)
}
if(file.exists('../metadata_web/meta_web_tcga.Rdata')) {
    m <- meta_web
    load('../metadata_web/meta_web_tcga.Rdata')
    meta_web <- rbind(m, meta_web)
}
save(meta_web, file = 'meta_web.Rdata')

## Add FC-RC files
# system('cp meta_web.Rdata meta_web_original.Rdata')
load('meta_web_original.Rdata')
load('../fileinfo/fc_rc_url_table.Rdata')

meta_web$'FANTOM-CAT' <- NA
for(i in seq_len(nrow(meta_web))) {
    proj <- gsub('.*">|</a>', '', meta_web$accession[i])
    if(proj %in% fc_rc_url_table$project) {
        
        proj_file <- paste0('rse_fc_', proj, '.rda')
        
        res <- paste0(
            '<a href="', fc_rc_url_table$url[fc_rc_url_table$file_name == proj_file],'" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-fc-rc\', 1)">RSE</a>'
        )
        if(proj %in% c('SRP012682', 'TCGA')) {
            ## Find tissue files (ignore the main one that we added already)
            proj_tab <- subset(fc_rc_url_table, project == proj & file_name != proj_file)
            
            ## Add tissues
            res <- paste0(res, ' RSE by tissue: ' , paste(paste0(
                '<a href="', proj_tab$url, 
                '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'data-rse-fc-rc\', 1)">',
                gsub('_', ' ', gsub(paste0('.*_', proj, '_|.rda'), '', proj_tab$file_name)),
                '</a>'
            ), collapse = ' '))
        }
        meta_web$'FANTOM-CAT'[i] <- res
    }
}
save(meta_web, file = 'meta_web.Rdata')

dir.create('www', showWarnings = FALSE)
system('scp e:/dcl01/leek/data/tcga_work/tcga_recount_junctions/sample_ids.tsv www/')


library('rsconnect')
load('.deploy_info.Rdata')
rsconnect::setAccountInfo(name=deploy_info$name, token=deploy_info$token,
    secret=deploy_info$secret)
deployApp(appFiles = c('ui.R', 'server.R', 'meta_web.Rdata', 'download.md',
    'google-analytics.js', 'www/sample_ids.tsv',
    'www/highlight.pack.js', 'www/shiny-showcase.js', 'www/rstudio.css',
    'www/LICENSE.txt', 'reducedexons.md'),
    appName = 'recount')#, account = deploy_info$name)
Y
