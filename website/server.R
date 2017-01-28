library('shiny')
library('DT')
library('markdown')
library('shinyBS')

load('meta_web.Rdata')

not_massive <- meta_web$number_samples < 9662
popular_i <- meta_web$number_samples > 400 & not_massive
gtex_i <- meta_web$number_samples == 9662
tcga_i <- meta_web$number_samples == 11284
meta_web$species <- factor(meta_web$species)
colnames(meta_web)[colnames(meta_web) == 'number_samples'] <- 'number of samples'
colnames(meta_web)[colnames(meta_web) == 'files_info'] <- 'files info'



shinyServer(function(input, output, session) {
    createAlert(session, 'annotationalert', 'gencode', 'Annotation change: now using Gencode v25', 'Gene and exon counts are now based on Gencode v25 CHR annotation. This adds non-protein coding genes to recount! We will no longer support the UCSC hg38 knownGene annotation.')
    
    output$metadata <- DT::renderDataTable(
        meta_web[not_massive, ],
        escape = which(colnames(meta_web) %in% c('number of samples', 'species',
            'abstract')),
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 10,
            lengthMenu = c(5, 10, 25, 50, 100, nrow(meta_web)),
            order = list(list(1, 'desc'))
        )
    )
    output$popular <- DT::renderDataTable(
        meta_web[popular_i, ],
        escape = which(colnames(meta_web) %in% c('number of samples', 'species',
            'abstract')),
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 10,
            lengthMenu = c(5, 10, 25),
            order = list(list(1, 'desc'))
        )
    )
    output$gtex <- DT::renderDataTable(
        meta_web[gtex_i, ],
        escape = which(colnames(meta_web) %in% c('number of samples', 'species',
            'abstract')),
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 1,
            lengthMenu = c(1, 5),
            order = list(list(1, 'desc'))
        )
    )
    output$tcga <- DT::renderDataTable(
        meta_web[tcga_i, ],
        escape = which(colnames(meta_web) %in% c('number of samples', 'species',
            'abstract')),
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 1,
            lengthMenu = c(1, 5),
            order = list(list(1, 'desc'))
        )
    )
    
    output$downloadData <- downloadHandler(
        filename = function() { paste0('recount_selection_', Sys.time(),
            '.csv') },
        content = function(file) {
            current <- meta_web[not_massive, c('accession',
                'number of samples', 'species',
                'abstract')][input$metadata_rows_all, ]
            current$accession <- gsub('</a>', '', gsub('.*">', '',
                current$accession))
            write.csv(current, file, row.names = FALSE)
        }
    )
})
