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

# Check whether script needs to run
if [[ -f "${outdir}/${pool_id}/${pool_id}_filtered.bam" ]]; then
  echo "`date` ${pool_id}.bam already present"
  exit 0
fi

bams=$(find ${outdir}/star_end2end -maxdepth 2 -name "*.Aligned.sortedByCoord.out.bam" -print)

mkdir -p "${outdir}/${pool_id}/"

echo "create merged BAM"
apptainer exec -B /hpc:/hpc,${TMPDIR}:${TMPDIR} ${container_dir}/samtools-1.12.sif samtools merge -@ ${threads} "${TMPDIR}/${pool_id}.bam" ${bams}
echo "filter merged BAM"
apptainer exec -B /hpc:/hpc,${TMPDIR}:${TMPDIR} ${container_dir}/samtools-1.12.sif samtools view -@ ${threads} -b -q 5 "${TMPDIR}/${pool_id}.bam" | \
  apptainer exec -B /hpc:/hpc ${container_dir}/samtools-1.12.sif samtools sort -@ ${threads} > "${outdir}/${pool_id}/${pool_id}_filtered.bam"
echo "index filtered merged BAM"
apptainer exec -B /hpc:/hpc ${container_dir}/samtools-1.12.sif samtools index -@ ${threads} "${outdir}/${pool_id}/${pool_id}_filtered.bam"
