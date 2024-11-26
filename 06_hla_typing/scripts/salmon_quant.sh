#!/bin/bash 

set -uo pipefail

cpu=$((${SLURM_CPUS_PER_TASK} * 2))

# Load files
mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
mapfile -t r2_files < ${project_folder}/documentation/r2_files.txt
mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

r1_file="${r1_files[$((SLURM_ARRAY_TASK_ID-1))]}"
r2_file="${r2_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

r1_filename=$(basename ${r1_file})
r2_filename=$(basename ${r2_file})
r1_trimmed="${r1_filename/_R1_/_R1_trimmed_}"
r2_trimmed="${r2_filename/_R2_/_R2_trimmed_}"

# Check if "quant.sf" file in salmon output dir already exists for this sample, if so, skip this sample
if [[ -f "${outdir}/salmon_quant/${sample_id}/quant.sf" ]]; then
  echo "quant.sf file already exists for ${sample_id}, skipping..."
  exit 0
fi

echo "`date` Running salmon quant for ${sample_id}"

mkdir -p "${outdir}/salmon_quant/${sample_id}/"
# Run salmon for transcript counts
apptainer  exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/salmon-1.8.0.sif salmon quant \
  --libType "A" \
  --validateMappings \
  --gcBias \
  --quiet \
  --numGibbsSamples 30 \
  --threads ${cpu} \
  -i "${salmon_index}" \
  -1 "${outdir}/trimgalore/${sample_id}/${r1_trimmed}" \
  -2 "${outdir}/trimgalore/${sample_id}/${r2_trimmed}" \
  --output "${outdir}/salmon_quant/${sample_id}/"

echo "`date` Finished salmon ${sample_id}"
