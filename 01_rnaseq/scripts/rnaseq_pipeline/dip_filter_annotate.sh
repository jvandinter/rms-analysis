#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:   24-06-2022
#
######################################################################

# Load parameters from main script
source $1

# Load correct modules
module load R/${r_version}

# Check whether script needs to run
if [[ -f "${wd}/data/processed/customannotation/${merged_gtf_basename}_novel_filtered.gtf" ]]; then
  echo "Merged filtered annotated GTF already present"
  exit 0
fi

# Create output dirs
mkdir -p "${wd}/data/processed/customannotation"

echo "`date` Filter and annotate novel GTF"

# Process and filter novel GTF
Rscript "${scriptdir}/dip_filter_annotate.R" \
  "${wd}" \
  "${reference_gtf}" \
  "${wd}/data/processed/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_gffcompare.annotated.gtf" \
  "${refseq_gtf}" \
  "${wd}/data/processed/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_matching.tracking" \
  "${min_occurence}" \
  "${min_tpm}" \
  "${wd}/data/processed/customannotation/${merged_gtf_basename}_novel_filtered.gtf"

echo "`date` Finished GTF filtering"