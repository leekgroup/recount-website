#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N fc_rc.fileinfo
#$ -o ./logs/fc_rc.fileinfo.txt
#$ -e ./logs/fc_rc.fileinfo.txt

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript get_fileinfo_fc_rc.R

echo "**** Job ends ****"
date