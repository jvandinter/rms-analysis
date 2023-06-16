#!/bin/bash

#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --gres=tmpspace:50G

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

# Load parameters from main script
wd="/hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts"
r_version=4.1.2
scriptdir=${wd}/scripts
twobit="/hpc/pmc_vanheesch/shared_resources/GENOMES/Homo_sapiens.GRCh38/102/Homo_sapiens.GRCh38.dna.primary_assembly.2bit"
merged_gtf_basename="EWS_full"

echo "`date` Generating custom annotation ..."

# Load correct modules
module load R/${r_version}

# Prepare RiboseQC and ORFquant annotation files based on 
# merged annotated stringtie GTF
Rscript "${scriptdir}/dip_prepare_custom_annotation.R" \
  ${twobit} \
  "${wd}/data/processed/customannotation/${merged_gtf_basename}_novel_filtered_corrected.gtf" \
  "${wd}/data/processed/customannotation/${merged_gtf_basename}/" \
  ${merged_gtf_basename}

echo "`date` Generating custom annotation finished"
