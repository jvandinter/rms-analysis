#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 21-12-2023
#
######################################################################

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Check whether script needs to run
if [[ -f "${outdir}/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz" ]]; then
  echo "`date` ${sample_id} FASTQ already present"
  exit 0
fi

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

mkdir -p "${outdir}/bowtie2/${sample_id}/"

echo "create FASTQ for ${sample_id}"
apptainer exec -B /hpc:/hpc ${container_dir}/samtools-1.12.sif samtools fastq \
  -@ ${threads} \
  "${outdir}/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam" > "${outdir}/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz"
