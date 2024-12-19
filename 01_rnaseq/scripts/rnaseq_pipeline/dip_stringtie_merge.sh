#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:  24-06-2022
#
######################################################################

# Load parameters from main script
source $1
source $2
cpu=$3
gtfmergefile="${wd}/data/processed/stringtie/gtfmergefile.txt"

# Load correct modules
module load stringtie/${stringtie_version}
module load gffcompare/${gffcompare_version}

# Check whether script needs to run
if [[ -f "${wd}/data/processed/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf" ]]; then
  echo "Merged GTF already present"
  exit 0
fi

# Create output dirs
mkdir -p "${wd}/data/processed/gffcompare/${merged_gtf_basename}/"
mkdir -p "${wd}/data/processed/stringtie/${merged_gtf_basename}/"
cd "${wd}/data/processed/stringtie/${merged_gtf_basename}/"

#! we should probably experiment with minimum isoform fraction and minimum tx length

echo "`date` running Stringtie --merge"

# Run stringtie merge
stringtie --version
stringtie ${gtfmergefile} \
  --merge \
  -G ${reference_gtf} \
  -f ${iso_frac} \
  -m 50 \
  -T 5 \
  -o "${wd}/data/processed/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf" \
  -p ${cpu}

echo "`date` finished Stringtie Merge"

echo "`date` running GFFcompare for annotation"

# Create output dirs
cd "${wd}/data/processed/gffcompare/${merged_gtf_basename}/"

# Run GFFcompare to annotate novel assembly .GTF
gffcompare \
  -V \
  -r ${reference_gtf} \
  -s ${masked_fasta} \
  -o "${merged_gtf_basename}_gffcompare" \
  "${wd}/data/processed/stringtie/${merged_gtf_basename}/${merged_gtf_basename}.gtf"

echo "`date` running GFFcompare for transcript occurence"

# Run GFFcompare to generate .tracking file

# Get correct sample
get_samples $wd $fastq_file $data_folder

for i in ${!sample_ids[@]}; do
  gtf_list[i]="${wd}/data/processed/stringtie/${sample_ids[i]}/${sample_ids[i]}.gtf"
done

gffcompare \
  -V \
  -r "${wd}/data/processed/gffcompare/${merged_gtf_basename}/${merged_gtf_basename}_gffcompare.annotated.gtf" \
  -o "${merged_gtf_basename}_matching" \
  ${gtf_list[@]}

echo "`date` finished running GFFcompare"