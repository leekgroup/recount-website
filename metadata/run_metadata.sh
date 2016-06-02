#!/bin/sh

# Usage
# sh run_metadata.sh sra
# sh run_metadata.sh gtex

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
else
    echo "Specify a valid project: gtex, sra"
fi


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.metadata"
    
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=130G,h_vmem=150G
#$ -N ${sname}

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

