#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 24-06-2022
#
######################################################################

# Load parameters from main script
source $1
source $2

# Load correct modules
module load gffcompare/${gffcompare_version}

# Get correct files
get_samples ${wd} ${data_folder}
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Check whether script needs to run
if [[ -f "${wd}/data/processed/gffcompare/${sample_id}/${sample_id}.stats" ]]; then
  echo "`date` GTF stats already present"
  exit 0
fi

# Create output dirs
mkdir -p "${wd}/data/processed/gffcompare/${sample_id}"
cd "${wd}/data/processed/gffcompare/${sample_id}/"

# Run GFFcompare to compare sample GTF to reference GTF
gffcompare --version
gffcompare \
  -r ${reference_gtf} \
  -o "${sample_id}" \
  "${wd}/data/processed/stringtie/${sample_id}/${sample_id}.gtf"


echo "`date` GFFCompare for ${sample_id} finished!"