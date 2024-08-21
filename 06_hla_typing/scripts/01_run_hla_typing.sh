#!/bin/bash

# 01_run_hla_typing.sh

set -uo pipefail

function usage() {
    cat <<EOF
SYNOPSIS
  01_run_hla_typing.sh [-c <config file>] [-h]
DESCRIPTION
  Run the HLA typing pipeline consisting of the following steps:
  1. Find samples from BAM files in the specified folder.
  2. Submit jobs for HLA typing using 'arcas-HLA' and merge outputs.
  Input Requirement: Requires BAM files containing aligned RNA-seq reads stored or symlinked in /projectfolder/data/bams/
OPTIONS
  -c, --config <file>    Configuration file to use
  -h, --help             Display this help message
AUTHOR
  Damon Hofman, MSc
EOF
}

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c|--config)
    CONFIG="$2"
    shift
    shift
    ;;
    -h|--help)
    usage
    exit
    ;;
    "")
    echo "Error: no option provided"
    usage
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    usage
    exit 1
    ;;
esac
done

# Check that configuration file is provided
if [[ -z ${CONFIG+x} ]]; then 
    echo "Error: no configuration file provided"
    usage
    exit 1
fi

# Load configuration variables
source $CONFIG

# Load general functions
source ${scriptdir}/general_functions.sh

# Create a unique prefix for the names for this run_id of the pipeline
# This makes sure that run_ids can be identified
run_id=$(uuidgen | tr '-' ' ' | awk '{print $1}')

# Find samples
echo "$(date '+%Y-%m-%d %H:%M:%S') Finding samples..."
get_samples $project_data_folder $data_folder $paired_end

printf "%s\n" "${r1_files[@]}" > ${project_folder}/documentation/r1_files.txt
printf "%s\n" "${r2_files[@]}" > ${project_folder}/documentation/r2_files.txt
printf "%s\n" "${sample_ids[@]}" > ${project_folder}/documentation/sample_ids.txt

# Create output directories
mkdir -p ${project_folder}/log/${run_id}/{trimgalore,star,samtools,arcas-hla,arcasmerge,salmon}
mkdir -p ${outdir}

##############################################################################

# Step 1a: quality control
trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.trimgalore \
  --output=${project_folder}/log/${run_id}/trimgalore/%A_%a.out \
  --export=ALL \
  "${scriptdir}/trimgalore.sh"
))
if [[ ${#trim_jobid[@]} -eq 0 ]]; then
  fatal "TrimGalore job not submitted successfully, trim_jobid array is empty"
fi
info "TrimGalore jobid: ${trim_jobid[@]}"

# Step 1b: Quantify with SALMON
salmon_quant_jobid=()

salmon_quant_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=6 \
  --time=36:00:00 \
  --gres=tmpspace:50G \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.salmon_quant \
  --output=${project_folder}/log/${run_id}/salmon/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  --export=ALL \
  ${scriptdir}/salmon_quant.sh
))

# Step 2: Align with STAR
star_jobid=()

star_jobid+=($(sbatch --parsable \
  --mem=80G \
  --cpus-per-task=12 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.star_align \
  --output=${project_folder}/log/${run_id}/star/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  --export=ALL \
  "${scriptdir}/star_align.sh"
))
info "STAR alignment jobid: ${star_jobid[@]}"

# Step 3: Create bam index file and generate mapping statistics with samtools
samtools_jobid=()

samtools_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=12 \
  --gres=tmpspace:100G \
  --time=24:00:00 \
  --job-name=${run_id}.samtools \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --output=${project_folder}/log/${run_id}/samtools/%A_%a.out \
  --dependency=aftercorr:${star_jobid} \
  --export=ALL \
  "${scriptdir}/samtools.sh"
))
info "Samtools jobid: ${samtools_jobid[@]}"

# Step 4: Run arcas-HLA
hla_jobid=()

hla_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=2 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.arcashla \
  --output=${project_folder}/log/${run_id}/arcas-hla/%A_%a.out \
  --dependency=aftercorr:${samtools_jobid} \
  --export=ALL \
  "${scriptdir}/run_arcas-hla.sh"
))

# Step 5: Merge arcas-HLA outputs
merge_jobid=()

merge_jobid+=($(sbatch --parsable \
  --cpus-per-task=2 \
  --time=1:00:00 \
  --job-name=${run_id}.arcasmerge \
  --output=${project_folder}/log/${run_id}/arcasmerge/%A_%a.out \
  --dependency=afterok:${hla_jobid} \
  --export=ALL \
  "${scriptdir}/merge_arcas_outputs.sh"
))
