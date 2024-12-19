#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

set -uo pipefail

# Load parameters from main script
source $1
source $2

# Load correct modules
module load cutadapt/${cutadapt_version}
module load fastqc/${fastqc_version}
module load trimgalore/${trimgalore_version}

# Get correct files
get_samples $wd $data_folder

full_path_fastq_1="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
full_path_fastq_2=$(echo $full_path_fastq_1 | sed 's/_R1/_R2/')

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

bf1=$(basename ${full_path_fastq_1})
bf2=$(basename ${full_path_fastq_2})

echo "`date` running trimgalore for ${sample_id}"

# Create output dirs
cd "${wd}/data/processed"
mkdir -p "trimgalore/${sample_id}/"

# Check whether script needs to run
if [[ -f "${wd}/data/processed/trimgalore/${sample_id}/${bf1}" ]]; then
  echo "`date` ${sample_id} file already present"
  exit 0
fi

# Run trimgalore on both reads
cutadapt --version
fastqc --version

cd "$wd/data/processed/trimgalore/${sample_id}/"

trim_galore "${full_path_fastq_1}" "${full_path_fastq_2}" \
  --cores 2 \
  --paired \
  --gzip \
  --fastqc \
  --fastqc_args "--outdir ${wd}/data/processed/trimgalore/${sample_id}/" \
  --output_dir "$wd/data/processed/trimgalore/${sample_id}"

# Change names of validated trimgalore output to basename
mv "${wd}/data/processed/trimgalore/${sample_id}/${bf1%.*.*}_val_1.fq.gz" "${wd}/data/processed/trimgalore/${sample_id}/${bf1}"
mv "${wd}/data/processed/trimgalore/${sample_id}/${bf2%.*.*}_val_2.fq.gz" "${wd}/data/processed/trimgalore/${sample_id}/${bf2}"

# Calculate trimcounts per paired fastq
tot_reads=$(zcat "${full_path_fastq_1}" | echo $((`wc -l`/4)))

trimmed_reads=$(zcat "${wd}/data/processed/trimgalore/${sample_id}/${bf1}" | echo $((`wc -l`/4)))
trimmed_percentage=`awk -vn=248 "BEGIN{print(${trimmed_reads}/${tot_reads}*100)}"`

# Add read trimming info to run QC file
printf '%s\t%s\t%s\t%s\n' "${sample_id}" "Trimmed" $trimmed_reads $trimmed_percentage >> "$wd/data/processed/trim_stats.txt"

echo "`date` finished ${sample_id}"