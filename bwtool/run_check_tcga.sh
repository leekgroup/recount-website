#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=3G
#$ -pe local 10
#$ -N check-tsv-tcga
#$ -hold_jid tcga.bwtool
#$ -o ./logs/check-tsv-tcga.txt
#$ -e ./logs/check-tsv-tcga.txt

echo '**** Job starts ****'
date

## Make logs dir
mkdir -p logs

## Annotate regions
module load R/3.3.x
Rscript check_tsv_tcga.R

echo '**** Job ends ****'
date
