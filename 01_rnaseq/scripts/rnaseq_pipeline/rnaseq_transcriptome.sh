#!/bin/bash

#SBATCH -t 4:00:00
#SBATCH --job-name=rnaseq_transcript_detection

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
######################################################################

function usage() {
    cat <<EOF
SYNOPSIS
  dip_transcripts.sh ./path/to/config.config - run novel isoform detection pipeline
  dip_transcripts.sh help - display this help message
DESCRIPTION
  1. Assemble transcriptomes for each sample with STRINGTIE, using BAM files as input
  2. Run GFFCOMPARE to compare custom GTF to reference GTF 
  3. Merge GTF files per sample group with STRINGTIE --merge
  4. Annotate and filter merged GTF with GFFCOMPARE and custom R scripts
  5. Create custom annotation with custom merged GTF for ribo-seq analysis
  6. Generate Salmon index
  7. Quantify reads with salmon quant
  8. Run MultiQC on new data

AUTHOR
  Jip van Dinter, MSc
  Damon Hofman, MSc
EOF
}

# Create a unique prefix for the names for this run of the pipeline. 
# This makes sure that runs can be identified
run=$(uuidgen | tr '-' ' ' | awk '{print $1}')

# Source all variables from the config file
CONFIG=$1

# Show help message if there was no config file location given on the commandline
if [[ -z $1 ]]; then usage; exit; fi

source ${CONFIG}
source ${scriptdir}/dip_functions.sh

################################################################################
#
# Find fastq samples in directory
#
################################################################################

mkdir -p log/${run}/{gffcompare,stringtie,salmon_quant}
get_samples $wd $fastq_file $data_folder

# make sure there are samples
if [[ ${#samples[@]} -eq 0 ]]; then
  fatal "no samples found in ./data/raw/"
fi

info "samples: n = ${#samples[@]}"
for i in ${!samples[@]}; do
  let "number=i+1"
  info "${number}    ${samples[i]}"
done

################################################################################
#
# Run the pipeline
#
################################################################################

echo -e "\n ====== `date` Detect Transcript Isoform Pipeline ====== \n"

echo -e "\n`date` Assembling transcriptome with Stringtie ..."
echo -e "====================================================================================== \n"

# 1. STRINGTIE: Assembles transcriptome of each sample with Stringtie
#               to locate potential novel transcript isoforms and
#               novel gene loci

stringtie_jobid=()

stringtie_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.stringtie \
  --output=log/${run}/stringtie/%A_%a.out \
  ${scriptdir}/dip_stringtie.sh \
  ${scriptdir}/dip_functions.sh \
  ${CONFIG} \
  ${medium_cpu}
))

info "Stringtie jobid: ${stringtie_jobid[@]}"

echo -e "\n`date` Checking precision and sensitivity with gffcompare ..."
echo -e "====================================================================================== \n"

#2. GFFCOMPARE: Compares the newly created GTFs with the reference GTF

gffcompare_jobid=()

gffcompare_jobid+=($(sbatch --parsable \
  --mem=${low_mem} \
  --cpus-per-task=${low_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.gffcompare \
  --output=log/${run}/gffcompare/%A_%a \
   --dependency=aftercorr:${stringtie_jobid} \
  ${scriptdir}/dip_gffcompare.sh \
  ${scriptdir}/dip_functions.sh \
  ${CONFIG}
))

info "Gffcompare jobid: ${gffcompare_jobid[@]}"

echo -e "\n`date` Merging transcriptomes with Stringtie ..."
echo -e "====================================================================================== \n"

# 3. STRINGTIE merge: Merge all single-sample GTFs into a single GTF
#    by keeping overlapping transcripts, using the reference GTF as 
#    a guide.

stringtie_merge_jobid=()

stringtie_merge_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${medium_cpu} \
  --time=24:00:00 \
  --job-name=${run}.stringtie_merge \
  --output=log/${run}/%A_stringtie_merge.out \
  --dependency=afterany:${stringtie_jobid} \
  ${scriptdir}/dip_stringtie_merge.sh \
  ${CONFIG} \
  ${scriptdir}/dip_functions.sh \
  ${medium_cpu}
))

info "Stringtie merge jobid: ${stringtie_merge_jobid[@]}"

echo -e "\n`date` Annotating and filtering novel assembly GTF ..."
echo -e "====================================================================================== \n"

# 4. GFFCOMPARE + R: Annotate the transcripts in the merged GTF
#                    with their associated gffcompare class codes
#                    and filter the transcripts based on strand,
#                    number of exons, and the transcript sample 
#                    occurence.

filter_annotate_jobid=()

filter_annotate_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${low_cpu} \
  --time=24:00:00 \
  --job-name=${run}.filter_annotate \
  --output=log/${run}/%A_filter_annotate.out \
  --dependency=afterany:${stringtie_merge_jobid} \
  ${scriptdir}/dip_filter_annotate.sh \
  ${CONFIG}
))

info "Filter and annotation jobid: ${filter_annotate_jobid[@]}"

echo -e "\n`date` Create custom annotation for downstream RIBO-seq ..."
echo -e "====================================================================================== \n"

# 5. R: Creates a custom annotation file using the merged, filtered,
#       annotated GTF for further RIBO-seq pipeline analysis

custom_annotation_jobid=()

if [[ ${create_annotation} =~ "TRUE" ]]; then

  if [[ $(find ${wd}/ -name '*.Rannot' | wc -l) -eq 0 ]]; then

    custom_annotation_jobid+=($(sbatch --parsable \
      --mem=${low_mem} \
      --cpus-per-task=${low_cpu} \
      --gres=tmpspace:50G \
      --time=24:00:00 \
      --job-name=${run}.custom_annotation \
      --output=log/${run}/%A_custom_annotation.out \
      --dependency=afterok:${filter_annotate_jobid} \
      ${scriptdir}/dip_custom_annotation.sh \
      ${CONFIG}
    ))

    info "Custom annotation jobid: ${custom_annotation_jobid[@]}"

  else
  echo "Annotation file already present"

  fi

else
  echo "Creation of annotation not specified"

fi

# 6. SALMON INDEX: Use the custom GTF file to generate a 
#                   salmon index for use with the Salmon quant tool

salmon_index_jobid=()

salmon_index_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --job-name=${run}.salmon_index \
  --output=log/${run}/%A_salmon_index.out \
  --dependency=afterok:${filter_annotate_jobid} \
  ${scriptdir}/dip_salmon_index.sh \
  ${CONFIG} \
  ${scriptdir}/dip_functions.sh \
  ${high_cpu}
))

info "Salmon index jobid: ${salmon_index_jobid[@]}"

echo -e "\n`date` Calculating counts with salmon ..."
echo -e "====================================================================================== \n"

# 7. SALMON: Use the trimmed fastq files from trimgalore
#             to estimate transcript abundance, now using
#             the merged, filtered, annotated GTF to
#             also include novel transcript isoforms.

salmon_jobid=()

salmon_jobid+=($(sbatch --parsable \
  --mem=${medium_mem} \
  --cpus-per-task=${high_cpu} \
  --time=24:00:00 \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.salmon_quant \
  --output=log/${run}/salmon_quant/%A_%a.out \
  --dependency=afterok:${salmon_index_jobid} \
  ${scriptdir}/dip_salmon_quant.sh \
  ${CONFIG} \
  ${scriptdir}/dip_functions.sh \
  ${high_cpu}
))

info "Salmon jobid: ${salmon_jobid[@]}"


# 5. MultiQC: 

multiqc_jobid=()

multiqc_jobid+=($(sbatch --parsable \
  --mem=${low_mem} \
  --cpus-per-task=${low_cpu} \
  --time=24:00:00 \
  --job-name=${run}.multiqc \
  --output=log/${run}/%A_multiqc.out \
  --dependency=afterok:${samtools_jobid} \
  ${scriptdir}/dip_multiqc.sh \
  ${CONFIG} \
  "transcripts"
))

info "MultiQC jobid: ${multiqc_jobid[@]}"

echo -e "\n ====== `date` Started all jobs! ====== \n"