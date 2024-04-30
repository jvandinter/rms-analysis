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

echo "$(date '+%Y-%m-%d %H:%M:%S') Finding Annotation..."
check_annotation ${reference_annotation} ${reference_gtf} ${reference_annotation_package} ${custom_annotation} ${custom_gtf} ${custom_annotation_package}

echo "`date` using ${annotation_package}"
echo "`date` using ${rannot}"
echo "`date` using ${gtf}"

export annotation_package=${annotation_package}
export rannot=${rannot}
export gtf=${gtf}

# Create output directories
mkdir -p ${project_folder}/log/${run_id}
echo "`date` using run ID: ${run_id}"
mkdir -p ${outdir}

################################################################################
#
# Run the pipeline
#
################################################################################

mkdir -p ${project_folder}/log/${run_id}

# echo -e "\n`date` Combine & filter bam with SAMTOOLS ..."
# echo -e "====================================================================================== \n"

# samtools_jobid=()

# samtools_jobid+=($(sbatch --parsable \
#   --mem=64G \
#   --cpus-per-task=8 \
#   --time=24:00:00 \
#   --job-name=${run_id}.samtools \
#   --gres=tmpspace:150G \
#   --output=${project_folder}/log/${run_id}/samtools.out \
#   --export=ALL \
#   ${scriptdir}/mrp_pooling.sh
# ))

# info "samtools jobid: ${samtools_jobid}"

# echo -e "\n`date` Perform QC with RiboseQC ..."
# echo -e "====================================================================================== \n"

# riboseqc_jobid=()

# riboseqc_jobid+=($(sbatch --parsable \
#   --mem=80G \
#   --cpus-per-task=1 \
#   --time=24:00:00 \
#   --job-name=${run_id}.riboseqc \
#   --output=${project_folder}/log/${run_id}/riboseqc_merged.out \
#   --dependency=afterok:${samtools_jobid} \
#   --export=ALL \
#   ${scriptdir}/mrp_merged_riboseqc.sh
# ))

# info "RiboseQC jobid: ${riboseqc_jobid}"

# echo -e "\n`date` Creating RiboseQC reports ..."
# echo -e "====================================================================================== \n"

# # 5. RiboseQC report. 

# riboreport_jobid=()

# riboreport_jobid+=($(sbatch --parsable \
#   --mem=80G \
#   --cpus-per-task=1 \
#   --time=24:00:00 \
#   --job-name=${run_id}.riboreport \
#   --output=${project_folder}/log/${run_id}/riboreport.out \
#   --dependency=afterok:${riboseqc_jobid} \
#   --export=ALL \
#   ${scriptdir}/mrp_riboseqc_report.sh
# ))

# info "RiboseQC report jobid: ${riboreport_jobid}"


echo -e "\n`date` Run ORFquant ..."
echo -e "====================================================================================== \n"

orfquant_jobid=()

orfquant_jobid+=($(sbatch --parsable \
  --mem=120G \
  --cpus-per-task=24 \
  --time=144:00:00 \
  --job-name=${run_id}.orfquant \
  --output=${project_folder}/log/${run_id}/orfquant_merged.out \
  --export=ALL \
  ${scriptdir}/mrp_merged_orfquant.sh
))

info "ORFquant job id: ${orfquant_jobid}"
