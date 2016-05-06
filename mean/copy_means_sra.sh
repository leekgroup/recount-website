#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=50G,h_vmem=60G,h_fsize=100G
#$ -N copy_means

echo "**** Job starts ****"
date

## Copy means
cd /dcl01/leek/data/recount-website/mean/means_sra
rsync -av /dcl01/leek/data/gtex_work/runs/recount2/mean/means_sra/ .

echo "**** Job ends ****"
date

## Move log files
mv /dcl01/leek/data/recount-website/mean/copy_means.* /dcl01/leek/data/recount-website/mean/logs/
