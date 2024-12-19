#!/bin/bash

#SBATCH --cpus-per-task 4
#SBATCH --job-name create_salmon_index
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --gres=tmpspace:20G
#SBATCH --array=1-10%10

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

set -uo pipefail

# Load parameters from main script
cpu=8
gffread_version=0.12.6
salmon_version=1.8.0
wd=`pwd`
resource_dir="/hpc/pmc_vanheesch/shared_resources"
species="Homo_sapiens"
genome_version="GRCh38"
annot="102"
scriptdir="${wd}/scripts"
reference_genome="/${resource_dir}/GENOMES/${species}.${genome_version}/${annot}/${species}.${genome_version}.dna.primary_assembly.fa"

gtf_basenames=(NBL_complete RMS_full ATRT_full EPN_full EWS_full AML_full B-ALL_full T-ALL_full WT_full MBL_full)

merged_gtf_basename=${gtf_basenames[$((SLURM_ARRAY_TASK_ID-1))]}

# Load correct modules
module load gffread/${gffread_version}
module load salmon/${salmon_version}

# Create output dirs
mkdir -p "${wd}/data/processed/salmon_index/"

# Create fasta from merged annotated GTF
gffread \
  -w "${wd}/data/processed/salmon_index/${merged_gtf_basename}_transcripts.fa" \
  -g ${reference_genome} \
  "${wd}/data/processed/customannotation/${merged_gtf_basename}_novel_filtered_corrected.gtf"

# Create "gentrome" for salmon index
gentrome="${wd}/data/processed/salmon_index/gentrome_${merged_gtf_basename}.fa"
cat "${wd}/data/processed/salmon_index/${merged_gtf_basename}_transcripts.fa" ${reference_genome} > ${gentrome}

# Grep seq names from fa and create decoys for salmon
decoy="${wd}/data/processed/salmon_index/decoys_${merged_gtf_basename}.txt"
grep "^>" ${reference_genome} | cut -d " " -f 1 > ${decoy}
sed -i.bak -e 's/>//g' ${decoy}

# Create salmon index 
salmon index \
  --transcripts ${gentrome} \
  --decoys ${decoy} \
  --index "${wd}/data/processed/salmon_index/${merged_gtf_basename}" \
  --threads ${cpu}

echo "`date` Finished creating salmon index"
