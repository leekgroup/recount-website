#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=20G,h_vmem=21G
#$ -pe local 10
#$ -N tcga_prep_meta

echo "**** Job starts ****"
date

## Create TCGA metadata
module load R/3.3.x
module load wiggletools/default
Rscript tcga_meta.R

echo "**** Job ends ****"
date

## Move log files
mv tcga_prep_meta.* logs/
