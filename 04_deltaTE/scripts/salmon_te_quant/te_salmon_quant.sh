#!/bin/bash

# This script quantifies reads mapping to ORFs using Salmon

set -uo pipefail

# Load files
mapfile -t sample_ids < sample_ids.txt

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
r1_filtered="${outdir}/bowtie2/${sample_id}/${sample_id}_filtered.fastq.gz"
index="${outdir}/salmon_index/"

mkdir -p "${outdir}/salmon_quant/${sample_id}/"

echo "`date` Running salmon for ${sample_id}"

# Calculate length and length standard deviation for salmon
fld=($(zcat ${r1_filtered} | awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sqrt(sq/n-m*m));}' - ))

# Run salmon for transcript counts
apptainer exec -B /hpc:/hpc ${container_dir}/salmon-1.8.0.sif salmon quant \
  salmon quant \
    --libType "A" \
    --validateMappings \
    --fldMean ${fld[2]#*=} \
    --fldSD ${fld[3]#*=} \
    --gcBias \
    --quiet \
    --numGibbsSamples 30 \
    --threads ${cpu} \
    -i "${index}" \
    -r "${r1_filtered}" \
    --output "${outdir}/salmon_quant/${sample_id}/"
  
echo "`date` Finished salmon ${sample_id}"
