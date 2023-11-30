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

mkdir -p "${outdir}/macs2/"

cd ${TMPDIR}

# RD cell line
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/macs2-2.2.7.1.sif macs2 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs2" \
        --name "RD_PAX3-FOXO1" \
        --treatment "${outdir}/picard/SRR039129/SRR039129_dedup.bam" "${outdir}/picard/SRR039133/SRR039133_dedup.bam" \
        --control "${outdir}/picard/SRR039130/SRR039130_dedup.bam"

# RH cell line
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/macs2-2.2.7.1.sif macs2 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs2" \
        --name "RH4_PAX3-FOXO1" \
        --keep-dup all \
        --treatment "${outdir}/picard/SRR039132/SRR039132_dedup.bam" "${outdir}/picard/SRR039135/SRR039135_dedup.bam" \
        --control "${outdir}/picard/SRR039131/SRR039131_dedup.bam" "${outdir}/picard/SRR039134/SRR039134_dedup.bam"
