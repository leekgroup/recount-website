## Load libraries
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading introns_unique.Rdata'))
system.time( load('introns_unique.Rdata') )
## Get junctions in the data
message(paste(Sys.time(), 'creating jx.tsv file'))
jx_df <- data.frame(chr = seqnames(introns_unique),
    start = start(introns_unique), end = end(introns_unique),
    strand = strand(introns_unique))

## Write to disk
write.table(jx_df, file = 'jx.tsv', sep = '\t', quote = FALSE,
    row.names = FALSE, col.names = FALSE)

## Reproducibility info
proc.time()
options(width = 120)
session_info()
