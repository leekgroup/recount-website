# qrsh -l mem_free=2G,h_vmem=3G -pe local 10
# module load R/3.3.x

library('R.utils')
library('BiocParallel')

## Parallel environment
bp <- MulticoreParam(workers = 10, outfile = Sys.getenv('SGE_STDERR_PATH'))
tsv <- dir('/dcl01/leek/data/recount2/coverage', pattern = 'tsv', full.names = TRUE)
names(tsv) <- gsub('.sum.tsv', '', dir('/dcl01/leek/data/recount2/coverage', pattern = 'tsv'))

system.time( tsv_lines <- bplapply(tsv, countLines, BPPARAM = bp) )
all(tsv_lines == 603139)

save(tsv_lines, file = 'tsv_lines.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
