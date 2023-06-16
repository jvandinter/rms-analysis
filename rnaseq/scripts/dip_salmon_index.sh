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
cpu=$3

# Load correct modules
module load gffread/${gffread_version}
module load salmon/${salmon_version}

# Create output dirs
mkdir -p "${wd}/data/processed/salmon_index/"

# Create fasta from merged annotated GTF
gffread \
  -w "${wd}/data/processed/salmon_index/${merged_gtf_basename}_transcripts.fa" \
  -g ${reference_genome} \
  "${wd}/data/processed/customannotation/${merged_gtf_basename}_novel_filtered.gtf"

# Create "gentrome" for salmon index
gentrome="${wd}/data/processed/salmon_index/gentrome.fa"
cat "${wd}/data/processed/salmon_index/${merged_gtf_basename}_transcripts.fa" ${reference_genome} > ${gentrome}

# Grep seq names from fa and create decoys for salmon
decoy="${wd}/data/processed/salmon_index/decoys.txt"
grep "^>" ${reference_genome} | cut -d " " -f 1 > ${decoy}
sed -i.bak -e 's/>//g' ${decoy}

# Create salmon index 
salmon index \
  --transcripts ${gentrome} \
  --decoys ${decoy} \
  --index "${wd}/data/processed/salmon_index/${merged_gtf_basename}" \
  --threads ${cpu}

echo "`date` Finished creating salmon index"