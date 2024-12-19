#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
#
# Date: 24-06-2022
#
######################################################################

set -uo pipefail

# Load parameters from main script
source $1
source $2
cpu=$3

# Load correct modules
module load STAR/${star_version}

# Get correct files
get_samples $wd $data_folder

fastq="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"
new_bam="${sample_id}.Aligned.sortedByCoord.out.bam"

read_1=$(basename ${fastq})
read_2=$(echo $read_1 | sed 's/_R1/_R2/')

# Create output dirs
cd "${wd}/data/processed"
mkdir -p "star/${sample_id}/"

echo "`date` running STAR for ${sample_id}"

# Check whether script needs to run
if [[ -s "${wd}/data/processed/star/${sample_id}/${sample_id}.Aligned.out.bam" ]]; then
  echo "`date` ${sample_id} BAM already present"
  exit 0
elif [[ -s "${wd}/data/processed/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam.bai" ]]; then
  echo "`date` ${sample_id} BAM index already present"
  exit 0
fi

# Check first 10k reads for read length for star index
read_length=$(gunzip -c "${fastq}" | \
  head -n 10000 | \
  awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count - 1}')

# Make sure read length picks are correct
if [[ "${read_length}" =~ "." ]]; then
  read_length=${read_length%.*}
fi

if (( ${read_length} >= 70 && ${read_length} <= 85 )); then
  used_index=74
elif (( ${read_length} >= 86 && ${read_length} <= 111 )); then
  used_index=99
elif (( ${read_length} >= 112 && ${read_length} <= 137 )); then
  used_index=124
elif (( ${read_length} >= 138 && ${read_length} <= 163 )); then
  used_index=149
else
  used_index=99
fi

# Use STAR for mapping the reads
STAR --version
STAR --genomeDir "${star_index_basedir}/${used_index}nt" \
  --sjdbGTFfile ${reference_gtf} \
  --readFilesIn "${wd}/data/processed/trimgalore/${sample_id}/${read_1}" "${wd}/data/processed/trimgalore/${sample_id}/${read_2}" \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --runThreadN ${cpu} \
  --runDirPerm All_RWX \
  --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}   \
  --outFilterType BySJout \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM NM MD jM jI MC ch \
  --outSAMstrandField intronMotif \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "${wd}/data/processed/star/${sample_id}/${sample_id}." \
  --outFilterMismatchNmax 6 \
  --outTmpKeep None \
  --alignSJoverhangMin 10 \
  --outFilterMultimapNmax 10 \
  --outFilterScoreMinOverLread 0.75

  echo "`date` finished mapping ${sample_id}"
