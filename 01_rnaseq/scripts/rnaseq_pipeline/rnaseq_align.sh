#!/bin/bash

#SBATCH -t 4:00:00
#SBATCH --job-name=rnaseq_aligner

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

set -uo pipefail

function usage() {
    cat <<EOF
SYNOPSIS
  dip_mapping.sh ./path/to/config.config - run part 1 novel isoform detection pipeline
  dip_mapping.sh help - display this help message
DESCRIPTION
  1. Run TRIMGALORE on paired reads for each sample
  a. Trim and filter reads with CUTADAPT
  b. Run FASTQC for each sample
  2. Run STAR on paired reads for each sample
  a. Check correct STAR index
  b. Map reads on human reference genome
  3. Create index and check mapping with SAMTOOLS
  4. Generate raw counts for exploratory analysis with FEATURECOUNTS
  5. Run MULTIQC for standard report output

AUTHOR
  Jip van Dinter, MSc
  Damon Hofman, MSc
EOF
}

# Create a unique prefix for the names for this run of the pipeline. 
# This makes sure that runs can be identified
run=$(uuidgen | tr '-' ' ' | awk '{print $1}')

# Source all variables from the config file
CONFIG=$1

# Show help message if there was no config file location given on the commandline
if [[ -z $1 ]]; then usage; exit; fi

source ${CONFIG}
source ${scriptdir}/dip_functions.sh

################################################################################
#
# Find fastq samples in directory
#
################################################################################

get_samples $wd $data_folder

# make sure there are samples
if [[ ${#samples[@]} -eq 0 ]]; then
  fatal "no samples found in ./data/raw/"
fi

info "samples: n = ${#samples[@]}"
for i in ${!samples[@]}; do
  let "number=i+1"
  info "${number}    ${samples[i]}"
done

################################################################################
#
# Run the pipeline
#
################################################################################

mkdir -p log/${run}/{trimgalore,star,samtools,howarewestrandedhere} data/processed analysis

echo -e "\n ====== `date` Detect Transcript Isoform Pipeline ====== \n"

echo -e "\n`date` Filtering and trimming ..."
echo -e "====================================================================================== \n"

# 1a. TRIMGALORE. Parallel start of all trimgalore jobs to filter for quality
#                with CUTADAPT and output quality reports with FASTQC

trim_jobid=()

trim_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.trimgalore \
  --output=log/${run}/trimgalore/%A_%a.out \
  ${scriptdir}/dip_trimgalore.sh \
  ${scriptdir}/dip_functions.sh \
  ${CONFIG}
))

info "Trimgalore jobid: ${trim_jobid[@]}"

echo -e "\n`date` Checking strandedness ..."
echo -e "====================================================================================== \n"

# 1b. Howarewestrandedhere

strand_jobid=()

strand_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${low_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.howarewestrandedhere \
  --output=log/${run}/howarewestrandedhere/%A_%a.out \
  ${scriptdir}/dip_check_strand.sh \
  ${scriptdir}/dip_functions.sh \
  ${CONFIG}
))

info "Howarewestrandedhere jobid: ${strand_jobid[@]}"

echo -e "\n`date` Mapping with STAR ..."
echo -e "====================================================================================== \n"

# 2. STAR: Parallel start of all STAR mapping jobs after quality filtering with 
#          CUTADAPT and quality control with FASTQC.

star_jobid=()

star_jobid+=($(sbatch --parsable \
  --mem=${high_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.star_align \
  --output=log/${run}/star/%A_%a.out \
  --dependency=aftercorr:${trim_jobid} \
  ${scriptdir}/dip_star_align.sh \
  ${scriptdir}/dip_functions.sh \
  ${CONFIG} \
  ${high_cpu}
))

info "STAR alignment jobid: ${star_jobid[@]}"

echo -e "\n`date` Checking mapping job with samtools ..."
echo -e "====================================================================================== \n"

# 3. SAMTOOLS: Create bam index file and generate mapping statistics.

samtools_jobid=()

samtools_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${high_cpu} \
  --gres=tmpspace:100G \
  --time=24:00:00 \
  --job-name=${run}.samtools \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --output=log/${run}/samtools/%A_%a.out \
  --dependency=aftercorr:${star_jobid} \
  ${scriptdir}/dip_samtools.sh \
  ${scriptdir}/dip_functions.sh \
  ${CONFIG} \
  ${medium_cpu}
))

info "Samtools jobid: ${samtools_jobid[@]}"

echo -e "\n`date` Creating exploratory PCA plot ..."
echo -e "====================================================================================== \n"

# 4. Featurecounts: Count all reads per reference transcript, then load the
#                   the counts in R for DEseq2 analysis and the generation
#                   of a PCA plot.

featurecounts_jobid=()

featurecounts_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --job-name=${run}.featurecounts \
  --output=log/${run}/%A_featurecounts.out \
  --dependency=afterok:${samtools_jobid} \
  ${scriptdir}/dip_featurecounts.sh \
  ${scriptdir}/dip_functions.sh \
  ${CONFIG} \
  ${medium_cpu}
))

info "FeatureCounts jobid: ${featurecounts_jobid[@]}"

echo -e "\n`date` Creating MultiQC file ..."
echo -e "====================================================================================== \n"

# 5. MultiQC: 

multiqc_jobid=()

multiqc_jobid+=($(sbatch --parsable \
  --mem=${low_mem} \
  --cpus-per-task=${low_cpu} \
  --time=24:00:00 \
  --job-name=${run}.multiqc \
  --output=log/${run}/%A_multiqc.out \
  --dependency=afterok:${samtools_jobid} \
  ${scriptdir}/dip_multiqc.sh \
  ${CONFIG} \
  "genes"
))

info "MultiQC jobid: ${multiqc_jobid[@]}"

echo -e "\n ====== `date` Started all jobs! ====== \n"
