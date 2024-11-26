#!/bin/bash

#SBATCH --job-name=computeMatrix
#SBATCH --output=log/deeptools/heatmap_%A_%a.out
#SBATCH --array=1-36%36
#SBATCH --cpus-per-task=2
#SBATCH --gres=tmpspace:10G
#SBATCH --mem=20G
#SBATCH --time=02:00:00

set -euo pipefail

output_dir="/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/05_chipseq/analysis/deeptools"

bigwigs=($(find /hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/05_chipseq/analysis/bamcoverage -type f ))
bigwig=${bigwigs[$SLURM_ARRAY_TASK_ID - 1]}

beds=$(find ${output_dir} -type f -name "*.bed" )

bed_test="/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/05_chipseq/analysis/deeptools/RMS_FP_novel_enriched.bed"

gtf="/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/01_rnaseq/analysis/rnaseq_pipeline/customannotation/RMS_full_novel_filtered_corrected.sorted.gtf"

gtfs=($(find "${output_dir}/gtfs" -type f -name "*.gtf" ))

name=$(basename ${bigwig} _normalized.bw)

mkdir -p "${output_dir}/${name}"

cd ${TMPDIR}
for gtf in ${gtfs[@]}; do
  gtf_name=$(basename ${gtf} .gtf)
  echo "computing matrix for $gtf_name"
  apptainer exec -B /hpc:/hpc/,${TMPDIR}:${TMPDIR} /hpc/local/Rocky8/pmc_vanheesch/singularity_images/deeptools_3.5.4--pyhdfd78af_1.sif computeMatrix \
    reference-point \
    --scoreFileName ${bigwig} \
    --regionsFileName ${gtf} \
    --outFileName ${output_dir}/${name}/${name}_${gtf_name}_meta.gz \
    --outFileNameMatrix ${output_dir}/${name}/${name}_${gtf_name}_meta_matrix.txt \
    --referencePoint "TSS" \
    --upstream 2000 \
    --downstream 2000 \
    --smartLabels \
    --metagene \
    -p 4

  apptainer exec -B /hpc:/hpc/,${TMPDIR}:${TMPDIR} /hpc/local/Rocky8/pmc_vanheesch/singularity_images/deeptools_3.5.4--pyhdfd78af_1.sif plotHeatmap \
    --matrixFile ${output_dir}/${name}/${name}_${gtf_name}_meta.gz \
    --outFileName ${output_dir}/${name}/${name}_${gtf_name}_meta_hm.png \
    --colorMap rocket \
    --zMin "-2" \
    --zMax "2" \
    --yMin "-1" \
    --yMax "1" \
    --regionsLabel "${gtf_name}"
done
