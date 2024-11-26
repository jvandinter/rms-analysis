#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 04-07-2023
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
    "")
    echo "Error: no option provided"
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    exit 1
    ;;
esac
done

# Check that configuration file is provided
if [[ -z ${CONFIG+x} ]]; then 
    echo "Error: no configuration file provided"
    exit 1
fi

# Load configuration variables
source $CONFIG

# Load general functions
source ${scriptdir}/mrp_functions.sh

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

if [[ "${check_samples}" == "true" ]]; then

  get_samples $project_data_folder $data_folder

  printf "%s\n" "${r1_files[@]}" > ${project_folder}/documentation/r1_files.txt
  printf "%s\n" "${sample_ids[@]}" > ${project_folder}/documentation/sample_ids.txt

  # make sure there are samples
  if [[ ${#samples[@]} -eq 0 ]]; then
    fatal "no samples found in ./raw/ or file containing fastq file locations not present"
  fi

  info "samples: n = ${#samples[@]}"
  for sample in ${samples[@]}; do
    info "    $sample"
  done

else

  mapfile -t r1_files < ${project_folder}/documentation/r1_files.txt
  mapfile -t sample_ids < ${project_folder}/documentation/sample_ids.txt

  # make sure there are samples
  if [[ ${#sample_ids[@]} -eq 0 ]]; then
    fatal "no samples found in ./raw/ or file containing fastq file locations not present"
  fi

  info "Samples: n = ${#sample_ids[@]}"
  for i in ${!sample_ids[@]}; do
    info "$((i+1))    ${sample_ids[i]}"
  done
fi 

echo "$(date '+%Y-%m-%d %H:%M:%S') Finding Annotation..."
check_annotation ${reference_annotation} ${reference_gtf} ${reference_annotation_package} ${custom_annotation} ${custom_gtf} ${custom_annotation_package}

echo "`date` using ${annotation_package}"
echo "`date` using ${rannot}"
echo "`date` using ${gtf}"

export annotation_package=${annotation_package}
export rannot=${rannot}
export gtf=${gtf}

# Create output directories
mkdir -p ${project_folder}/log/price_${run_id}/{trimgalore,star_endtoend,bowtie2} 
echo "`date` using run ID: ${run_id}"
mkdir -p ${outdir}

echo -e "\n`date` Run PRICE indexer ..."
echo -e "====================================================================================== \n"

price_annot_jobid=()

price_annot_jobid+=($(sbatch --parsable \
  --mem=8G \
  --cpus-per-task=1 \
  --time=24:00:00 \
  --job-name=${run_id}.price_annot \
  --output=${project_folder}/log/price_${run_id}/price_annot.out \
  --export=ALL \
  ${scriptdir}/mrp_price_annot.sh
))

info "price_annot jobid: ${price_annot_jobid}"

echo -e "\n`date` Filtering and trimming ..."
echo -e "====================================================================================== \n"

trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=8G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.trimgalore \
  --output=${project_folder}/log/price_${run_id}/trimgalore/%A_%a.out \
  --export=ALL \
  ${scriptdir}/mrp_trimgalore.sh 
))

if [[ ${#trim_jobid[@]} -eq 0 ]]; then
  fatal "TrimGalore job not submitted successfully, trim_jobid array is empty"
fi

info "trimgalore jobid: ${trim_jobid}"

echo -e "\n`date` Removing contaminants ..."
echo -e "====================================================================================== \n"

contaminant_jobid=()

contaminant_jobid+=($(sbatch --parsable \
  --mem=8G \
  --cpus-per-task=12 \
  --time=24:00:00 \
  --gres=tmpspace:100G \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.contaminants \
  --output=${project_folder}/log/price_${run_id}/bowtie2/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  --export=ALL \
  ${scriptdir}/mrp_remove_contaminants.sh
))

info "contaminant jobid: ${contaminant_jobid}"

echo -e "\n`date` Align reads to genome with STAR ..."
echo -e "====================================================================================== \n"

# 3. STAR. Align contaminant-depleted read files to supplied genome and
#          transcriptome. If no new custom transcriptome is supplied, 
#          the normal reference transcriptome is used for 
#          guided assembly.

star_jobid=()

star_jobid+=($(sbatch --parsable \
  --mem=60G \
  --cpus-per-task=8 \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run_id}.star_align \
  --output=${project_folder}/log/price_${run_id}/star_endtoend/%A_%a.out \
  --dependency=aftercorr:${contaminant_jobid} \
  --export=ALL \
  ${scriptdir}/mrp_star_endtoend.sh
))

info "alignment jobid: ${star_jobid}"

echo -e "\n`date` Combine & filter bam with SAMTOOLS ..."
echo -e "====================================================================================== \n"

samtools_jobid=()

samtools_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=4 \
  --time=24:00:00 \
  --job-name=${run_id}.samtools \
  --gres=tmpspace:50G \
  --output=${project_folder}/log/price_${run_id}/samtools.out \
  --dependency=afterok:${star_jobid} \
  --export=ALL \
  ${scriptdir}/mrp_samtools_merge.sh
))

info "samtools jobid: ${samtools_jobid}"

echo -e "\n`date` Run PRICE ..."
echo -e "====================================================================================== \n"

price_jobid=()

price_jobid+=($(sbatch --parsable \
  --mem=24G \
  --cpus-per-task=1 \
  --time=24:00:00 \
  --job-name=${run_id}.price \
  --output=${project_folder}/log/price_${run_id}/price.out \
  --dependency=afterok:${samtools_jobid},${price_annot_jobid} \
  --export=ALL \
  ${scriptdir}/mrp_price.sh
))

info "PRICE jobid: ${price_jobid}"

echo -e "\n ====== `date` Started all jobs! ====== \n"
