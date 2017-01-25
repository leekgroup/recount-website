#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=220G,h_vmem=300G,h_fsize=100G
#$ -N merge_sra
#$ -hold_jid sra.rse

echo "**** Job starts ****"
date

mkdir -p logs

## Merge all of SRA's RSE objects at the gene and exon levels
module load R/3.3.x
Rscript merge_all.R

echo "**** Job ends ****"
date

## Move log files
mv merge_sra.* logs/
