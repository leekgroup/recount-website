## Currently can only be deployed with R 3.2.5

## Load metadata information
load('../metadata_web/meta_web_sra.Rdata')
if(file.exists('../metadata_web/meta_web_gtex.Rdata')) {
    m <- meta_web
    load('../metadata_web/meta_web_gtex.Rdata')
    meta_web <- rbind(m, meta_web)
}
save(meta_web, file = 'meta_web.Rdata')


library('rsconnect')
load('.deploy_info.Rdata')
rsconnect::setAccountInfo(name=deploy_info$name, token=deploy_info$token,
    secret=deploy_info$secret)
deployApp(appFiles = c('ui.R', 'server.R', 'meta_web.Rdata', 'download.md',
    'www/ucsc-knowngene-hg38-genes-bp-length.Rdata',
    'www/ucsc-knowngene-hg38-exons.Rdata'), appName = 'recount')
