## Usage
# mkdir -p logs
# Rscript create_bed.R > logs/create_bed.txt 2>&1

## Extract exons for each gene
library('recount')
library('rtracklayer')

## Export exons as a BED file
export(unlist(recount_exons), con = 'Gencode-v25.bed', format='BED')

## Save how the exons are related, for speeding up the tsv -> count matrix step
## Group counts by gene
n <- elementNROWS(recount_exons)
count_groups <- rep(seq_len(length(n)), n)
save(count_groups, file = 'count_groups.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
