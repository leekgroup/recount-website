#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=3G
#$ -pe local 10
#$ -N check-tsv-gtex
#$ -hold_jid gtex.bwtool
#$ -o ./logs/check-tsv-gtex.txt
#$ -e ./logs/check-tsv-gtex.txt

echo '**** Job starts ****'
date

## Make logs dir
mkdir -p logs

## Annotate regions
module load R/3.3.x
Rscript check_tsv_gtex.R

echo '**** Job ends ****'
date
