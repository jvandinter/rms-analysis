#!/bin/bash

######################################################################
#
# Authors: l.w.kok-15@prinsesmaximacentrum.nl
# Date: 05-11-2021
#
######################################################################

# Show start time
start=$(date +"%d-%m-%Y %T")
echo "Start : $start"

# Load parameters from main script
allele="deprecated"

# Create output directory
mkdir -p "${outdir}/netMHCpan/${run}"

cd "${outdir}/netMHCpan/${run}"

# Convert the input to the input required by netMHCpan
apptainer exec -B "/hpc:/hpc" --env LC_CTYPE="en_US.UTF-8" ${apptainer_dir}/blast_parser_python-2.0.sif \
    python3 ${scriptdir}/app_prepare_input.py ${input_fasta} \
    "${outdir}/netMHCpan/${run}/peptides.fasta" ${protdir} ${allele}

# Show end time
end=$(date +"%d-%m-%Y %T")
echo "End : $end"
