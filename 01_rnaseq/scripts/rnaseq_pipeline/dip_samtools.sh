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
cpu=$3

# Get correct files
get_samples $wd $data_folder

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
bam="${sample_id}.Aligned.out.bam"
new_bam="${sample_id}.Aligned.sortedByCoord.out.bam"

# Load correct modules
module load samtools/${samtools_version}

echo "`date` running samtools for ${sample_id}"

# Check whether script needs to run
if [[ -f "${wd}/data/processed/star/${sample_id}/${new_bam}.bai" ]]; then
  echo "`date` ${sample_id} BAM index already present"
  exit 0
fi

if ! [[ -s "${wd}/data/processed/star/${sample_id}/${bam}" ]]; then
  echo "${sample_id} not mapped correctly: rerun"
  exit 1
fi

# Sort BAM
samtools sort \
  -@ ${cpu} \
  -l 9 \
  -o "${wd}/data/processed/star/${sample_id}/${new_bam}" \
  -T "${TMPDIR}" \
  "${wd}/data/processed/star/${sample_id}/${bam}"

# Create mapping statistics with samtools
samtools stats -@ ${cpu} \
  "${wd}/data/processed/star/${sample_id}/${new_bam}" > "${wd}/data/processed/star/${sample_id}/${sample_id}_stats.txt"

# Index the bam with samtools
samtools --version
samtools index -@ ${cpu} "${wd}/data/processed/star/${sample_id}/${new_bam}"

# Remove original STAR bam
if [[ -s "${wd}/data/processed/star/${sample_id}/${new_bam}" ]]; then
  rm "${wd}/data/processed/star/${sample_id}/${bam}"
fi

echo "`date` finished indexing ${sample_id}"
