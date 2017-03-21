#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=12G,h_fsize=100G
#$ -o ./chmod.txt
#$ -e ./chmod.txt

echo "**** Job starts sample \${bwsample} ****"
date

echo "Write-protecting (group)"
chmod g-w -R *
date

echo "Read-protecting (others)"
chmod o-r -R *
date

echo "**** Job ends ****"
date
