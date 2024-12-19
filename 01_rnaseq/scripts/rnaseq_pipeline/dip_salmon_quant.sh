#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

# Load parameters from main script
source $1
source $2
cpu=$3

# Load correct modules
module load salmon/${salmon_version}

# Get correct files
get_samples $wd $data_folder

sample_id=${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}

# Create output dirs
mkdir -p "${wd}/data/processed/salmon_quant/${sample_id}/"

# Check whether script needs to run
if [[ $(ls "${wd}/data/processed/salmon_quant/${sample_id}/" | wc -l) -gt 0 ]]; then
  echo "`date` Salmon count files already present"
  exit 0
fi

echo "`date` Running salmon quant for ${sample_id}"

read_1=$(basename ${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]})
read_2=$(echo $read_1 | sed 's/_R1/_R2/')

# Run salmon for transcript counts
salmon quant \
  --libType "A" \
  --validateMappings \
  --gcBias \
  --quiet \
  --numGibbsSamples 30 \
  --threads ${cpu} \
  -i "${wd}/data/processed/salmon_index/${merged_gtf_basename}" \
  -1 "${wd}/data/processed/trimgalore/${sample_id}/${read_1}" \
  -2 "${wd}/data/processed/trimgalore/${sample_id}/${read_2}" \
  --output "${wd}/data/processed/salmon_quant/${sample_id}"

echo "`date` Finished salmon ${sample_id}"