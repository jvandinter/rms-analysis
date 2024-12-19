#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

set -euo pipefail

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

# Load files
mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

# Create output dirs
cd "${outdir}"
mkdir -p "picard/${sample_id}/"

# Check whether script needs to run
if [[ -s "${outdir}/picard/${sample_id}/${sample_id}_dedup.bam" ]]; then
  echo "`date` ${sample_id} already present"
  exit 0
fi

cd $TMPDIR
# Mark Dups
echo "`date` : marking duplicates for $sample_id"
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/picard-2.20.1.sif java \
        -Xmx24G \
        -jar /opt/picard.jar \
        MarkDuplicates \
        I="${outdir}/bwa/${sample_id}/${sample_id}.bam" \
        O="${sample_id}.bam" \
        M="${outdir}/picard/${sample_id}/${sample_id}.MarkDuplicates.metrics.txt"

# Remove duplicates
echo "`date` : removing duplicates for $sample_id"
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/samtools-1.16.sif samtools view -@ ${threads} -b -F 0x400 "${sample_id}.bam" > "${outdir}/picard/${sample_id}/${sample_id}_dedup.bam"
cd "${outdir}/picard/${sample_id}/"
echo "`date` : indexing $sample_id"
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/samtools-1.16.sif samtools index -@ ${threads} "${outdir}/picard/${sample_id}/${sample_id}_dedup.bam"
