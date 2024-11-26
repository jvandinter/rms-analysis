#!/bin/bash

#SBATCH --cpus-per-task 4
#SBATCH --job-name create_salmon_index
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --gres=tmpspace:50G

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

# This script generates a salmon index for the NBL_complete dataset, based on the NBL_complete filtered corrected GTF file 
# generated with the rnaseq_customannotation scripts

# Default parameters
cpu=8
resource_dir="/hpc/pmc_vanheesch/shared_resources"
species="Homo_sapiens"
genome_version="GRCh38"
annot="102"
transcriptome_basename="${species}.${genome_version}.${annot}"
reference_genome="${resource_dir}/GENOMES/${species}.${genome_version}/${annot}/${species}.${genome_version}.dna.primary_assembly.fa"
reference_transcriptome="${resource_dir}/GENOMES/${species}.${genome_version}/${annot}/annotation/${species}.${genome_version}.${annot}.gtf"
index_dir="${resource_dir}/GENOMES/${species}.${genome_version}/${annot}/salmon/1.8.0"
wd="${index_dir}"

# Create fasta from merged annotated GTF
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 /hpc/local/Rocky8/pmc_vanheesch/singularity_images/gffread-0.12.6.sif gffread \
  -w "${wd}/${transcriptome_basename}_transcripts.fa" \
  -g ${reference_genome} \
  "${reference_transcriptome}"

# Create "gentrome" for salmon index
gentrome="${wd}/gentrome_${transcriptome_basename}.fa"
cat "${wd}/${transcriptome_basename}_transcripts.fa" ${reference_genome} > ${gentrome}

# Grep seq names from fa and create decoys for salmon
decoy="${wd}/decoys_${transcriptome_basename}.txt"
grep "^>" ${reference_genome} | cut -d " " -f 1 > ${decoy}
sed -i.bak -e 's/>//g' ${decoy}

# Create salmon index
apptainer  exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 /hpc/local/Rocky8/pmc_vanheesch/singularity_images/salmon-1.8.0.sif salmon index \
  --transcripts ${gentrome} \
  --decoys ${decoy} \
  --index "${index_dir}/" \
  --threads ${cpu}

echo "`date` Finished creating salmon index"
