#!/bin/bash

# Run the R script to make a table with several characteristics of the proteins

mkdir "${outdir}/characteristics_${run}" || { echo "Error: Failed to create the output directory."; exit 1; }

apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch,$TMPDIR:$TMPDIR,/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_4.3.0_libs:/usr/local/lib/R/site-library "${apptainer_dir}/rstudio_4.3.0_bioconductor.sif" \
    Rscript "${scriptdir}/characteristics/calculate_characteristics.R" \
        "${input_fasta}" \
        "${outdir}/characteristics_${run}/characteristics_proteins.tsv"