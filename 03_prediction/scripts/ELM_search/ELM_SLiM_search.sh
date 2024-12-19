#!/bin/bash

mkdir "${outdir}/ELM_search_${run}" || { echo "Error: Failed to create the output directory."; exit 1; }

# script that exectues the SLIM search in fasta and disorder files
apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch,$TMPDIR:$TMPDIR,/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_4.3.0_libs:/usr/local/lib/R/site-library "${apptainer_dir}/rstudio_4.3.0_bioconductor.sif" \
    Rscript "${scriptdir}/ELM_search/ELM_SLiM_search.R" \
        "${input_fasta}" \
        "${outdir}/iupred3_${run}/per_protein" \
        "${elm_classes}" \
        "${outdir}/ELM_search_${run}/ELM_SLiMs_overview.tsv"
