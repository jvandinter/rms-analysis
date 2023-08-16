#!/bin/bash
#SBATCH --job-name gtex_results
#SBATCH --mem=16G
#SBATCH --time=12:00:00

module load R/4.2.1

Rscript create_evo_res.R
