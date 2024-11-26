#!/bin/bash

set -uo pipefail

cpu=$((${SLURM_CPUS_PER_TASK} * 2))

# Load files
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
bam_file="${outdir}/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam"
bam_filename=$(basename ${bam_file})

# Create output dirs
arcas_outdir=${outdir}/arcas-hla/${sample_id}
mkdir -p "${arcas_outdir}/tmp"

# Run arcas-HLA
module load miniconda
# conda install arcas-hla

# 1. Extract chromosome 6 reads and related HLA sequences
arcasHLA \
  extract ${bam_file} \
  --outdir ${arcas_outdir} \
  --temp ${arcas_outdir}/tmp/ \
  --threads ${cpu} \
  --unmapped \
  --verbose

# 2. Run arcas-HLA genotype to predict most likely genotype (no partial alleles)
arcasHLA \
  genotype \
  --outdir ${arcas_outdir} \
  --temp ${arcas_outdir}/tmp/ \
  --threads ${cpu} \
  --verbose \
  ${arcas_outdir}/${bam_filename/.bam/.extracted}.1.fq.gz \
  ${arcas_outdir}/${bam_filename/.bam/.extracted}.2.fq.gz 
