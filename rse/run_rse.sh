#!/bin/sh

# Usage
# sh run_rse.sh sra
# sh run_rse.sh gtex

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
    MEM="mem_free=100G,h_vmem=180G,h_fsize=30G"
    CORES=""
elif [[ "${PROJECT}" == "gtex" ]]
then
    echo "$PROJECT"
    MEM="mem_free=84G,h_vmem=117G,h_fsize=50G"
    CORES="#$ -pe local 3"
else
    echo "Specify a valid project: gtex, sra"
fi

# Count how many commands there are
# For testing use: LINES=10
#LINES=10
LINES=$(wc -l ${MAINDIR}/metadata/project_ids_${PROJECT}.txt | cut -f1 -d " ")
METADATA="${MAINDIR}/metadata/metadata_${PROJECT}.Rdata"


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.rse"

for LINES in 95 192 211 221 223 231 233 235 238 239 244 246 256 297 299 300 303 325 346 352 365 375 385 388 413 414 425 442 457 475 489 499 514 519 525 532 533 568 613 617 626 646 700 718 719 820 850 932 943 1025 1027 1158 1312 1357 1405 1485 1491 1592 1725 1733 1995 2034 2035 2036 2038 2039
do
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m a
#$ -l ${MEM}
#$ -N ${sname}
#$ -t ${LINES}
#$ -hold_jid ${PROJECT}.metadata
${CORES}

PROJECTNAME=\$(awk "NR==\${SGE_TASK_ID}" ${MAINDIR}/metadata/project_ids_${PROJECT}.txt)

rm -fr rse_sra/${PROJECTNAME}
echo "**** Job starts project \${PROJECTNAME} ****"
date

## Run the R script that creates the count and rse files for exons and genes
module load R/3.3
Rscript create_rse.R -p "${PROJECT}" -m "${METADATA}" -i "\${PROJECTNAME}"

echo "**** Job ends ****"
date

## Move log files
mv ${WDIR}/${sname}.*.\${SGE_TASK_ID} ${WDIR}/logs/
EOF

    call="qsub ${WDIR}/.${sname}.sh"
    echo $call
    $call
done


