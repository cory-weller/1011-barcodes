#!/usr/bin/env bash

chromosome=${1}
start=${2}
stop=${3}

singularity exec src/singularity.sif Rscript src/query-region.R ${chromosome} ${start} ${stop}