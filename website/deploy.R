## Can now be deployed with R 3.3.0

## Load metadata information
load('../metadata_web/meta_web_sra.Rdata')
if(file.exists('../metadata_web/meta_web_gtex.Rdata')) {
    m <- meta_web
    load('../metadata_web/meta_web_gtex.Rdata')
    meta_web <- rbind(m, meta_web)
}
save(meta_web, file = 'meta_web.Rdata')


dir.create('www', showWarnings = FALSE)
system('cp ../genes/ucsc-knowngene-hg38-exons.Rdata www/')
system('cp ../genes/ucsc-knowngene-hg38-genes-bp-length.Rdata www/')


library('rsconnect')
load('.deploy_info.Rdata')
rsconnect::setAccountInfo(name=deploy_info$name, token=deploy_info$token,
    secret=deploy_info$secret)
deployApp(appFiles = c('ui.R', 'server.R', 'meta_web.Rdata', 'download.md',
    'www/ucsc-knowngene-hg38-genes-bp-length.Rdata',
    'www/ucsc-knowngene-hg38-exons.Rdata', 'google-analytics.js', 
    'www/TxDb.Hsapiens.UCSC.hg38.knownGene_3.1.3.tar.gz'),
    appName = 'recount', account = deploy_info$name)
