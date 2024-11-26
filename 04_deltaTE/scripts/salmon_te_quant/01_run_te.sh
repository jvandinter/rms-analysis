#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 23-06-2023
#
######################################################################

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
source ${scriptdir}/te_functions.sh

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
mkdir -p ${project_folder}/log/${run_id}/{trimgalore,salmon} 
echo "`date` using run ID: ${run_id}"
mkdir -p ${outdir}

# make sure there are samples
if [[ ${#samples[@]} -eq 0 ]]; then
  fatal "no samples found in ./raw/ or file containing fastq file locations not present"
fi

echo "samples: n = ${#samples[@]}"
for sample in ${samples[@]}; do
  echo "    $sample"
done

################################################################################
#
# Run the pipeline
#
################################################################################

echo -e "\n ====== `date` Delta TE Pipeline ====== \n"

echo -e "\n`date` Building index ..."
echo -e "====================================================================================== \n"

# 1. Salmon Index

index_jobid=()

index_jobid+=($(sbatch --parsable \
  --mem=120G \
  --cpus-per-task=2 \
  --time=24:00:00 \
  --job-name=${run_id}.index \
  --output=${project_folder}/log/${run_id}/index.out \
  --export=ALL \
  ${scriptdir}/te_salmon_index.sh
))

if [[ ${#index_jobid[@]} -eq 0 ]]; then
  fatal "TrimGalore job not submitted successfully, trim_jobid array is empty"
fi

echo "Salmon index jobid: ${index_jobid}"

echo -e "\n`date` Filtering and trimming ..."
echo -e "====================================================================================== \n"

# 1. TRIMGALORE. Parallel start of all trimgalore jobs to filter for quality
#                with CUTADAPT and output quality reports with FASTQC

trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=8G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.trimgalore \
  --output=${project_folder}/log/${run_id}/trimgalore/%A_%a.out \
  --dependency=afterok:${index_jobid} \
  --export=ALL \
  ${scriptdir}/te_trimgalore.sh 
))

echo "trimgalore jobid: ${trim_jobid}"

echo -e "\n`date` Removing contaminants ..."
echo -e "====================================================================================== \n"

# 2. BOWTIE2. Use combination of tRNA, rRNA, snRNA, snoRNA, mtDNA fastas to
#             remove those contaminants from RIBO-seq data. Outputs QC stats
#             to a file per contaminant group.

contaminant_jobid+=($(sbatch --parsable \
  --mem=8G \
  --cpus-per-task=12 \
  --time=24:00:00 \
  --gres=tmpspace:100G \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.contaminants \
  --output=${project_folder}/log/${run_id}/bowtie2/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  --export=ALL \
  ${scriptdir}/te_remove_contaminants.sh
))

echo -e "\n`date` Quantification ..."
echo -e "====================================================================================== \n"

# 2. Salmon. Quantifies filtered FASTQ against 

quant_jobid=()

quant_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=6 \
  --time=12:00:00 \
  --gres=tmpspace:50G \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.salmon_quant \
  --output=${project_folder}/log/${run_id}/salmon/%A_%a.out \
  --dependency=aftercorr:${contaminant_jobid} \
  --export=ALL \
  ${scriptdir}/te_salmon_quant.sh
))

echo "quantification jobid: ${quant_jobid}"

echo -e "\n ====== `date` Started all jobs! ====== \n"
