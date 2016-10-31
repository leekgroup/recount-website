#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=3G
#$ -pe local 25
#$ -N check-tsv-tcga
#$ -hold_jid tcga.bwtool

echo '**** Job starts ****'
date

## Make logs dir
mkdir -p logs

## Annotate regions
module load R/3.3
Rscript check_tsv_tcga.R

## Move log files into the logs directory
mv check-tsv-tcga.* logs/

echo '**** Job ends ****'
date
