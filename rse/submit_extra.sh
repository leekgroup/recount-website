#!/bin/bash

cd /dcl01/leek/data/recount-website/rse/rse_sra
qsub run_merse.sh
cd /dcl01/leek/data/recount-website/rse/rse_gtex
qsub run_split.sh
cd /dcl01/leek/data/recount-website/rse/rse_tcga
qsub run_split.sh
