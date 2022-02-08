#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 8G
#SBATCH -t 0-24:00:00
#SBATCH --account sadhumj
#SBATCH -p norm

chromosome=${1}
windowSize=${2}
windowStep=${3}

module load R/3.6.3

Rscript divergence.R ${chromosome} ${windowSize} ${windowStep}
