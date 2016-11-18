#!/bin/sh

# Usage
# sh run_metadata.sh sra
# sh run_metadata.sh gtex
# sh run_metadata.sh tcga

# Directories
MAINDIR=/dcl01/leek/data/recount-website
WDIR=${MAINDIR}/metadata

# Define variables
PROJECT=$1

# Create log dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell file

if [[ "${PROJECT}" == "sra" ]]
then
    echo "${PROJECT}"
elif [[ "${PROJECT}" == "gtex" ]]
then
    echo "${PROJECT}"
elif [[ "${PROJECT}" == "tcga" ]]
then
    echo "${PROJECT}"
else
    echo "Specify a valid project: gtex, sra, tcga"
fi


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.metadata"
    
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l rnet,mem_free=10G,h_vmem=11G
#$ -pe local 3
#$ -N ${sname}
#$ -hold_jid tcga_meta_prep

echo "**** Job starts ****"
date

## Create metadata
module load R/3.3.x
Rscript create_metadata.R -p "${PROJECT}"

echo "**** Job ends ****"
date

## Move log files
mv ${WDIR}/${sname}.* ${WDIR}/logs/
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call

