#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 15-11-2023
#
######################################################################

set -euo pipefail

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
source ${scriptdir}/chip_functions.sh

# Create a unique prefix for the names for this run_id of the pipeline
# This makes sure that run_ids can be identified
run_id=$(uuidgen | tr '-' ' ' | awk '{print $1}')

################################################################################
#
# Find fastq samples in directory
#
################################################################################

# Find samples
echo "$(date '+%Y-%m-%d %H:%M:%S') Finding samples..."
get_samples $project_data_folder $data_folder

printf "%s\n" "${r1_files[@]}" > ${project_folder}/documentation/r1_files.txt
printf "%s\n" "${sample_ids[@]}" > ${project_folder}/documentation/sample_ids.txt

# Create output directories
mkdir -p ${project_folder}/log/${run_id}/{trimgalore,bwa,picard,macs} 
echo "`date` using run ID: ${run_id}"
mkdir -p ${outdir}

# make sure there are samples
if [[ ${#samples[@]} -eq 0 ]]; then
  fatal "no samples found in ./raw/ or file containing fastq file locations not present"
fi

info "samples: n = ${#samples[@]}"
for sample in ${samples[@]}; do
  info "    $sample"
done

################################################################################
#
# Run the pipeline
#
################################################################################

echo -e "\n ====== `date` ChIP peak calling pipeline ====== \n"

echo -e "\n`date` Filtering and trimming ..."
echo -e "====================================================================================== \n"

#trim_jobid=()

#trim_jobid+=($(sbatch --parsable \
#  --mem=8G \
#  --cpus-per-task=4 \
#  --time=24:00:00 \
#  --array 1-${#samples[@]}%${simul_array_runs} \
#  --job-name=${run_id}.trimgalore \
#  --output=${project_folder}/log/${run_id}/trimgalore/%A_%a.out \
#  --export=ALL \
#  ${scriptdir}/chip_trimgalore.sh 
#))

echo -e "\n`date` Align reads to genome with bwa mem2 ..."
echo -e "====================================================================================== \n"

echo -e "\n`date` Call peaks with MACS ..."
echo -e "====================================================================================== \n"

macs_jobid=()

macs_jobid+=($(sbatch --parsable \
  --mem=12G \
  --cpus-per-task=1 \
  --time=24:00:00 \
  --gres=tmpspace:30G \
  --job-name=${run_id}.macs \
  --output=${project_folder}/log/${run_id}/macs/macs.out \
  --export=ALL \
  ${scriptdir}/chip_macs.sh
))

info "MACS jobid: ${macs_jobid}"

echo -e "\n ====== `date` Started all jobs! ====== \n"
