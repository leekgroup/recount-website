#!/bin/sh

# Usage
# sh run_rse.sh sra
# sh run_rse.sh gtex
# sh run_rse.sh tcga

# Directories
MAINDIR=/dcl01/leek/data/recount-website
WDIR=${MAINDIR}/rse

# Define variables
PROJECT=$1

# Create log dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell file

if [[ "${PROJECT}" == "sra" ]]
then
    echo "$PROJECT"
    MEM="leek,mem_free=50G,h_vmem=160G,h_fsize=50G"
    EMAIL="a"
elif [[ "${PROJECT}" == "gtex" ]]
then
    echo "$PROJECT"
    MEM="mem_free=100G,h_vmem=150G,h_fsize=50G"
    EMAIL="e"
elif [[ "${PROJECT}" == "tcga" ]]
then
    echo "$PROJECT"
    MEM="mem_free=100G,h_vmem=150G,h_fsize=50G"
    EMAIL="e"
else
    echo "Specify a valid project: gtex, sra, tcga"
fi

# Count how many commands there are
# For testing use: LINES=10
#LINES=10
LINES=$(wc -l ${MAINDIR}/metadata/project_ids_${PROJECT}.txt | cut -f1 -d " ")
METADATA="${MAINDIR}/metadata/metadata_${PROJECT}.Rdata"


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.rse"

cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m ${EMAIL}
#$ -l ${MEM}
#$ -N ${sname}
#$ -t 1:${LINES}
#$ -hold_jid ${PROJECT}.metadata
#$ -o ./logs/${PROJECT}.rse.o.\$TASK_ID.txt
#$ -e ./logs/${PROJECT}.rse.e.\$TASK_ID.txt

PROJECTNAME=\$(awk "NR==\${SGE_TASK_ID}" ${MAINDIR}/metadata/project_ids_${PROJECT}.txt)

echo "**** Job starts project \${PROJECTNAME} ****"
date

## Run the R script that creates the count and rse files for exons and genes
module load R/3.3.x
Rscript create_rse.R -p "${PROJECT}" -m "${METADATA}" -i "\${PROJECTNAME}"

echo "**** Job ends ****"
date
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call


