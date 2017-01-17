#!/bin/bash

mkdir -p logs
qsub run_check.sh
qsub run_check_gtex.sh
qsub run_check_tcga.sh
