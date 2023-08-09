#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

set -uo pipefail

# Load parameters from main script
source $1
source $2

# Load correct modules
module load howarewestrandedhere

# Get correct files
get_samples $wd $data_folder

full_path_fastq_1="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

full_path_fastq_1="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
full_path_fastq_2=$(echo $full_path_fastq_1 | sed 's/_R1/_R2/')

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Create output dirs
mkdir -p "${wd}/data/processed/check_strandedness"
cd "${wd}/data/processed/check_strandedness"

# Check whether script needs to run
if [[ -s "${wd}/data/processed/check_strandedness/${sample_id}.txt" ]]; then
  echo "`date` ${sample_id} file already present"
  exit 0
fi

# Infer strandedness
check_strandedness -g ${reference_gtf} \
-n 1000000 \
-r1 ${full_path_fastq_1} \
-r2 ${full_path_fastq_2} \
-k "${kallisto_index}" >> "${wd}/data/processed/check_strandedness/${sample_id}.txt"

touch "${wd}/data/processed/check_strandedness/strandedness_all.txt" 
strandedness=$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print $NF}')
printf "%s\t%s\n" "$sample_id" "$strandedness" >> "${wd}/data/processed/check_strandedness/strandedness_all.txt" 

echo "`date` finished mapping ${sample_id}"
