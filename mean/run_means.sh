#!/bin/sh

# Usage
# sh run_means.sh sra
# sh run_means.sh gtex
# sh run_means.sh tcga

# Directories
MAINDIR=/dcl01/leek/data/recount-website
WDIR=${MAINDIR}/mean

# Define variables
PROJECT=$1

# Create log dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell file

if [[ "${PROJECT}" == "sra" ]]
then
    echo "$PROJECT"
    MEM="mem_free=50G,h_vmem=60G,h_fsize=100G"
    CORES=1
elif [[ "${PROJECT}" == "gtex" ]]
then
    echo "$PROJECT"
    MEM="mem_free=100G,h_vmem=150G,h_fsize=100G"
    CORES=2
elif [[ "${PROJECT}" == "tcga" ]]
then
    echo "$PROJECT"
    MEM="mem_free=100G,h_vmem=150G,h_fsize=100G"
    CORES=2
else
    echo "Specify a valid project: gtex, sra, tcga"
fi

# Count how many commands there are
# For testing use: LINES=10
#LINES=10
LINES=$(wc -l ${MAINDIR}/metadata/project_ids_${PROJECT}.txt | cut -f1 -d " ")
METADATA="${MAINDIR}/metadata/metadata_${PROJECT}.Rdata"


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.mean"
    
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m a
#$ -l ${MEM}
#$ -N ${sname}
#$ -t 1:${LINES}
#$ -pe local ${CORES}
#$ -o ./logs/${PROJECT}.mean.o.\${TASK_ID}.txt
#$ -e ./logs/${PROJECT}.mean.e.\${TASK_ID}.txt

PROJECTNAME=\$(awk "NR==\${SGE_TASK_ID}" ${MAINDIR}/metadata/project_ids_${PROJECT}.txt)

echo "**** Job starts project \${PROJECTNAME} ****"
date

## Run the R script that creates mean bigwig files per project
module load R/3.3.x
module load ucsctools
module load wiggletools/default
Rscript calculate_means.R -p "${PROJECT}" -m "${METADATA}" -i "\${PROJECTNAME}"

echo "**** Job ends ****"
date
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call


