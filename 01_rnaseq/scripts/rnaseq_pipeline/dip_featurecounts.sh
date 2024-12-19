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
module load subread/${subread_version}
module load R/${r_version}

# Create output dir
mkdir -p "${wd}/data/processed/featurecounts"

echo "`date` running FeatureCounts for ${merged_gtf_basename} samples"

# Check whether script needs to run
if [[ -f "${wd}/data/processed/${merged_gtf_basename}.pdf" ]]; then
  echo "`date` PCA plot already present"
  exit 0
elif [[ -f "${wd}/data/processed/featurecounts/${merged_gtf_basename}.counts" ]]; then
  echo "`date` Count file already present"
  echo "`date` generating PCA for ${merged_gtf_basename} samples"

  # Load featurecounts table in R for DESeq2 analysis
  Rscript "${scriptdir}/dip_create_exploratory_PCA.R" \
    "${wd}" \
    "${wd}/data/processed/featurecounts/${merged_gtf_basename}.counts" \
    "${metadata}" \
    "${merged_gtf_basename}"

  echo "`date` Finished ${merged_gtf_basename} samples"
fi

# Get correct files
get_samples $wd $data_folder

# Check strandedness for downstream tools
for i in ${!sample_ids[@]}; do
  bams[i]="${wd}/data/processed/star/${sample_ids[i]}/${sample_ids[i]}.Aligned.sortedByCoord.out.bam"
  strandedness[i]=$( awk -v sid="${sample_ids[i]}" '$1 ~ sid {print $2}' ${wd}/data/processed/check_strandedness/strandedness_all.txt)
    if [[ ${strandedness[i]} == "RF/fr-firststrand" ]]; then
      strandtype[i]=2
    elif [[ ${strandedness[i]} == "FR/fr-secondstrand" ]]; then
      strandtype[i]=1
    else
      strandtype[i]=0
    fi
done

strandtype=$(IFS=,; echo "${strandtype[*]}")

echo "`date` Running featureCounts"

featureCounts -v

featureCounts \
  -s ${strandtype} \
  -T ${cpu} \
  -p \
  -t "CDS" \
  -g "gene_id" \
  -J \
  -G ${reference_genome} \
  -a "${reference_gtf}" \
  -o "${wd}/data/processed/featurecounts/${merged_gtf_basename}.counts" \
  ${bams[@]}

echo "`date` generating PCA for ${merged_gtf_basename} samples"

# Load featurecounts table in R for DESeq2 analysis
Rscript "${scriptdir}/dip_create_qc_pca.R" \
  "${wd}" \
  "${wd}/data/processed/featurecounts/${merged_gtf_basename}.counts" \
  "${merged_gtf_basename}"

echo "`date` Finished ${merged_gtf_basename} samples"