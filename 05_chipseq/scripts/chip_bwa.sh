#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

set -euo pipefail

# Load correct BWA version
module load bwa-mem2/${bwa_version}

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

# Load files
mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_file="${r1_files[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_filename=$(basename ${r1_file})

# Create output dirs
cd "${outdir}"
mkdir -p "bwa/${sample_id}/"

# Check whether script needs to run
if [[ -s "${outdir}/bwa/${sample_id}/${sample_id}.bam" ]]; then
  echo "`date` ${sample_id} already present"
  exit 0
fi

echo "`date` mapping ${sample_id}"

# # Align reads
# bwa-mem2 mem \
#   -t ${threads} \
#   -v 3 \
#   ${bwa_index} \
#   "${outdir}/trimgalore/${sample_id}/${r1_filename}" \
#   "${outdir}/bwa/${sample_id}/${sample_id}.bam"

  # Align reads
bwa-mem2 \
  mem \
  -t ${threads} \
  -v 3 \
  ${bwa_index} \
  "${outdir}/trimgalore/${sample_id}/${r1_filename}" | \
  apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/samtools-1.16.sif samtools \
  sort \
  --threads ${threads} \
  -o "${outdir}/bwa/${sample_id}/${sample_id}.bam" -
