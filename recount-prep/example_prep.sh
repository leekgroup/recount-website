#!/bin/sh

## Usage:
# bash example_prep.sh

## Example in JHPCE
## Have to run this manually, for now??
module load R/3.3.x
module load ucsctools
module load wiggletools/default

## Some common variables
DATADIR="/path/to/rail-rna/output-dir"
MANIFEST="/path/to/rail-rna/manifest/file"
BWTOOL="/path/to/bwtool"
WIGGLE="/path/to/wiggletools"
WIGTOBIGWIG="/path/to/wigToBigWig"

COUNTS=$DATADIR/cross_sample_results/counts.tsv.gz
## Download some required files
Rscript prep_setup.R

## Display help info on how to run prep_sample.R
Rscript prep_sample.R -h
## Note that if you setup the environment variable TMPDIR and export it, this
## will control where R stores the temporary files. See ?tempdir for more info.

## Process one sample individually:
# Rscript prep_sample.R -f ${DATADIR}/JH-13_GGCTAC_L006.bw -c ${COUNTS} -b ${BWTOOL} -w ${WIGGLE} -a TRUE
# Rscript prep_sample.R -f ${DATADIR}/JH-30_GTGAAA_L004.bw -c ${COUNTS} -b ${BWTOOL} -w ${WIGGLE} -a TRUE

## Use full arguments for another sample:
# Rscript prep_sample.R --bigwig_file ${DATADIR}/JH-11_GGCTAC_L003.bw --counts_file ${COUNTS} --bwtool ${BWTOOL} --wiggletools ${WIGGLE} --calculate_auc TRUE

## Process two samples at a time with parallel --jobs 2
## Note that you'll need to request 2 cores
parallel --jobs 2 Rscript prep_sample.R --bigwig_file ${DATADIR}/{} --counts_file ${COUNTS} --bwtool ${BWTOOL} --wiggletools ${WIGGLE} --calculate_auc TRUE ::: JH-13_GGCTAC_L006.bw JH-11_GGCTAC_L003.bw

## Now merge results
paste rse_temp/counts_exon_* > counts_exon.tsv
gzip counts_exon.tsv
paste rse_temp/counts_gene_* > counts_gene.tsv
gzip counts_gene.tsv

## Display help info on how to run prep_merge.R
Rscript prep_merge.R -h

## Merge rse objects and create junction rse object
BWDIR=$DATADIR/coverage_bigwigs
JUNCTIONS=$DATADIR/cross_sample_results/first_pass_junctions.tsv.gz
Rscript prep_merge.R -b ${BWDIR} -j ${JUNCTIONS} -m ${MANIFEST} -w ${WIGGLE} -t ${WIGTOBIGWIG} -c TRUE
