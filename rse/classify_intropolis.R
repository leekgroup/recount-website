## Load libraries
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading introns_unique.Rdata'))
system.time( load('introns_unique.Rdata') )

## Get junctions in the data
message(paste(Sys.time(), 'reading jx_in_intropolis.tsv file'))
jx_intropolis_raw <- read.table('jx_in_intropolis.tsv',
    col.names = c('count', 'chr', 'start', 'end', 'strand'),
    colClasses = c('integer', 'character', 'integer', 'integer', 'character'))

## Find which are in intropolis or not
jx_intropolis <- GRanges(seqnames = jx_intropolis_raw$chr,
    IRanges(jx_intropolis_raw$start, jx_intropolis_raw$end),
    strand = jx_intropolis_raw$strand)

introns_unique$in_intropolis <- countOverlaps(introns_unique,
    jx_intropolis, type = 'equal') > 0

print('Exon-exon junctions by presence in Intropolis')
table(introns_unique$in_intropolis)
round(table(introns_unique$in_intropolis) / length(introns_unique) * 100, 2)

## Save results
message(paste(Sys.time(),
    'saving results in introns_unique_with_intropolis.Rdata'))
save(introns_unique, file = 'introns_unique_with_intropolis.Rdata')

## Compress original data
system('gzip jx.tsv')
system('gzip jx_in_intropolis.tsv')

## Reproducibility info
proc.time()
options(width = 120)
session_info()
