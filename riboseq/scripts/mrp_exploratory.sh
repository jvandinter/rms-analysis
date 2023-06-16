#!/bin/bash

#SBATCH -t 4:00:00
#SBATCH --job-name=ribo-sebas

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 11-10-2022
#
######################################################################

function info() {
    echo "INFO: $@" >&2
}
function error() {
    echo "ERR:  $@" >&2
}
function fatal() {
    echo "ERR:  $@" >&2
    exit 1
}

# Create a unique prefix for the names for this run of the pipeline. 
# This makes sure that runs can be identified
run=$(uuidgen | tr '-' ' ' | awk '{print $1}')
echo "run ID ${run}"

# Show help message if there was no config file location given on the commandline
if [[ -z $1 ]]; then 

  exit;

fi

# Source all variables from the config file
CONFIG=$1
source ${CONFIG}
source ${scriptdir}/mrp_functions.sh

################################################################################
#
# Find fastq samples in directory
#
################################################################################

get_samples ${wd} ${data_folder}

check_annotation ${reference_annotation} ${reference_gtf} ${custom_annotation} ${custom_gtf} ${custom_annotation_package} ${reference_annotation_package}

echo "`date` using ${annotation_package}"
echo "`date` using ${rannot}"

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

mkdir -p log/${run}/{trimgalore,bowtie2,star_align,riboseqc} data/processed

echo -e "\n ====== `date` Map Riboseq Pipeline ====== \n"

echo -e "\n`date` Filtering and trimming ..."
echo -e "====================================================================================== \n"

# 1. TRIMGALORE. Parallel start of all trimgalore jobs to filter for quality
#                with CUTADAPT and output quality reports with FASTQC

trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.trimgalore \
  --output=log/${run}/trimgalore/%A_%a.out \
  ${scriptdir}/mrp_trimgalore.sh \
  ${scriptdir}/mrp_functions.sh \
  ${CONFIG}
))

info "trimgalore jobid: ${trim_jobid}"

echo -e "\n`date` Removing contaminants ..."
echo -e "====================================================================================== \n"

# 2. BOWTIE2. Use combination of tRNA, rRNA, snRNA, snoRNA, mtDNA fastas to
#             remove those contaminants from RIBO-seq data. Outputs QC stats
#             to a file per contaminant group.

contaminant_jobid=()

contaminant_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.contaminants \
  --output=log/${run}/bowtie2/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  ${scriptdir}/mrp_remove_contaminants.sh \
  ${scriptdir}/mrp_functions.sh \
  ${CONFIG} \
  ${medium_cpu}
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
  --mem=${high_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.star_align \
  --output=log/${run}/star_align/%A_%a.out \
  --dependency=aftercorr:${contaminant_jobid} \
  ${scriptdir}/mrp_star_align.sh \
  ${CONFIG} \
  ${gtf} \
  ${high_cpu}
))

info "alignment jobid: ${star_jobid}"

echo -e "\n`date` Perform QC with RiboseQC ..."
echo -e "====================================================================================== \n"

# 4. RiboseQC. 

riboseqc_jobid=()

riboseqc_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.riboseqc \
  --output=log/${run}/riboseqc/%A_%a.out \
  --dependency=aftercorr:${star_jobid} \
  ${scriptdir}/mrp_riboseqc.sh \
  ${CONFIG} \
  ${rannot} \
  ${annotation_package} \
  ${medium_cpu}
))

info "RiboseQC jobid: ${riboseqc_jobid}"

echo -e "\n`date` Perform QC with RiboseQC ..."
echo -e "====================================================================================== \n"

# 5. MultiQC.

multiqc_jobid=()

multiqc_jobid+=($(sbatch --parsable \
  --mem=${low_mem} \
  --cpus-per-task=${low_cpu} \
  --time=24:00:00 \
  --job-name=${run}.multiqc \
  --output=log/${run}/%A_multiqc.out \
  --dependency=afterok:${riboseqc_jobid} \
  ${scriptdir}/mrp_multiqc.sh \
  ${CONFIG} \
  "riboseq"
))

info "MultiQC jobid: ${multiqc_jobid[@]}"

echo -e "\n ====== `date` Started all jobs! ====== \n"
