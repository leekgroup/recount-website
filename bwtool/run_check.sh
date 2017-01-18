#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=3G
#$ -pe local 25
#$ -N check-tsv
#$ -hold_jid sra.bwtool
#$ -o ./logs/check-tsv.o.txt
#$ -e ./logs/check-tsv.e..txt

echo '**** Job starts ****'
date

## Make logs dir
mkdir -p logs

## Annotate regions
module load R/3.3.x
Rscript check_tsv.R

echo '**** Job ends ****'
date
