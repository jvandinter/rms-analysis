#!/bin/bash

#SBATCH --cpus-per-task 12 
#SBATCH --job-name merge_BAM
#SBATCH --mem=24G 
#SBATCH --time=24:00:00

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 04-06-2021
#
######################################################################

# Load parameters from main script
cpu=12
threads=$((cpu * 2))
wd=`pwd`

# Load software modules
module load samtools/1.12

pool_ids=("RNA_tumoroid_eRMS" "RNA_tumoroid_aRMS" "RNA_patient_eRMS" "RNA_patient_aRMS")
check_ids=("tum_erms" "tum_arms" "tis_erms" "tis_arms")

for i in ${!pool_ids[@]}; do

  # Create output dirs
  mkdir -p "data/processed/${pool_ids[i]}"

  # Find bam files to merge
  bams=$(find ${wd}/data/processed/star/${check_ids[i]} -maxdepth 2 -name "*.Aligned.sortedByCoord.out.bam" -print)

  # Merge bams
  echo -e "\n `date` Merging bam files ..."
  samtools merge -@ ${threads} "${wd}/data/processed/${pool_ids[i]}/${pool_ids[i]}.bam" $bams

done
echo -e "\n `date` Merging bam files ... complete! "
