#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

# Load parameters from main script
source $1

# Check whether script needs to run
if [[ -f "${wd}/data/processed/customannotation/${merged_gtf_basename}/${merged_gtf_basename}_novel_filtered.gtf_Rannot" ]]; then
  echo "Custom annotation package already generated!"
  exit 0
fi

echo "`date` Generating custom annotation ..."

# Load correct modules
module load R/${r_version}

# Prepare RiboseQC and ORFquant annotation files based on 
# merged annotated stringtie GTF
Rscript "${scriptdir}/dip_prepare_custom_annotation.R" \
  ${twobit} \
  "${wd}/data/processed/customannotation/${merged_gtf_basename}_novel_filtered.gtf" \
  "${wd}/data/processed/customannotation/${merged_gtf_basename}/" \
  ${merged_gtf_basename}

echo "`date` Generating custom annotation finished"