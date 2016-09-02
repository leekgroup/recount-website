#!/bin/sh

## Example in JHPCE
## Have to run this manually, for now??
#module load R/3.3.x
#module load ucsctools
#module load wiggletools/default

## Some common variables
DATADIR="/dcl01/leek/data/sunghee_analysis/processed/coverage_bigwigs/"
COUNTS="/dcl01/leek/data/sunghee_analysis/processed/cross_sample_results/counts.tsv.gz"
BWTOOL="/dcl01/leek/data/bwtool/bwtool-1.0/bwtool"
WIGGLE="wiggletools"

## Download some required files
Rscript prep_setup.R

## Display help info on how to run prep_sample.R
Rscript prep_sample.R -h TRUE

## Process a couple of samples (labeled as single-end)
Rscript prep_sample.R -f ${DATADIR}/JH-13_GGCTAC_L006.bw -c ${COUNTS} -b ${BWTOOL} -w ${WIGGLE} -p FALSE -a "TRUE"
#Rscript prep_sample.R -f ${DATADIR}/JH-11_GGCTAC_L003.bw -c ${COUNTS} -b ${BWTOOL} -w ${WIGGLE} -p FALSE -a "TRUE"
#Rscript prep_sample.R -f ${DATADIR}/JH-30_GTGAAA_L004.bw -c ${COUNTS} -b ${BWTOOL} -w ${WIGGLE} -p FALSE -a "TRUE"
