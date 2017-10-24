#!/bin/sh

# Usage
# sh run_fileinfo.sh sra
# sh run_fileinfo.sh gtex
# sh run_fileinfo.sh tcga

# Directories
MAINDIR=/dcl01/leek/data/recount-website
WDIR=${MAINDIR}/fileinfo

# Define variables
PROJECT=$1

# Create log dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell file

if [[ "${PROJECT}" == "sra" ]]
then
    echo "$PROJECT"
    MEM="mem_free=1G,h_vmem=2G,h_fsize=10G"
elif [[ "${PROJECT}" == "gtex" ]]
then
    echo "$PROJECT"
    MEM="mem_free=10G,h_vmem=15G,h_fsize=10G"
elif [[ "${PROJECT}" == "tcga" ]]
then
    echo "$PROJECT"
    MEM="mem_free=10G,h_vmem=15G,h_fsize=10G"
else
    echo "Specify a valid project: gtex, sra, tcga"
fi

# Count how many commands there are
# For testing use: LINES=10
#LINES=10
LINES=$(wc -l ${MAINDIR}/metadata/project_ids_${PROJECT}.txt | cut -f1 -d " ")
METADATA="${MAINDIR}/metadata/metadata_${PROJECT}.Rdata"


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.fileinfo"
    
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m a
#$ -l leek,${MEM}
#$ -N ${sname}
#$ -t 1:${LINES}
#$ -hold_jid copy_means,${PROJECT}.mean,${PROJECT}.rse,split_${PROJECT}
#$ -o ./logs/${PROJECT}.fileinfo.\$TASK_ID.txt
#$ -e ./logs/${PROJECT}.fileinfo.\$TASK_ID.txt

PROJECTNAME=\$(awk "NR==\${SGE_TASK_ID}" ${MAINDIR}/metadata/project_ids_${PROJECT}.txt)

echo "**** Job starts project \${PROJECTNAME} ****"
date

## Run the R script that computes md5sum and file sizes as well as the final
## list of files to upload
module load R/3.3.x
Rscript get_fileinfo.R -p "${PROJECT}" -m "${METADATA}" -i "\${PROJECTNAME}"

echo "**** Job ends ****"
date
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call


