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
    appName = 'recount', account = deploy_info$name)
