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
source $1
source $2
output_dir=$3

# Create output directory
mkdir -p ${processed}/${output_dir}

cd ${scriptdir}

# Convert the input to the input required by MHCFlurry
${scriptdir}/env/bin/python3 app_prepare_input.py ${input} \
    ${processed}/${output_dir}/"peptides.csv" ${protdir} ${allele}

# Show end time
end=$(date +"%d-%m-%Y %T")
echo "End : $end"
