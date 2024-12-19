#!/bin/bash 


indices=(/hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/NBL_complete /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/RMS_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/EPN_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/ATRT_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/EWS_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/B-ALL_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/T-ALL_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/AML_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/MBL_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/WT_full)
data_folder="/hpc/pmc_vanheesch/data/external/rnaseq"
cpu=12
wd=`pwd`

# Load correct modules
module load salmon/1.8.0
module load cutadapt/3.4
module load fastqc/0.11.9
module load trimgalore/0.6.6

# Get correct files
fastq_names=($(ls ${wd}/data/raw_single))
for i in ${!fastq_names[@]}; do

    fastq_files[i]=$(find "${data_folder}" -name ${fastq_names[i]} 2> >(grep -v 'Permission denied' >&2))
    full_name[i]=${fastq_files[i]%.fastq.gz}
    sample_ids[i]=$(basename ${full_name[i]})

done

full_path_fastq_1="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "`date` running ${sample_id} for ${index_name}"

bf1=$(basename ${full_path_fastq_1})

# Create output dirs


echo "`date` running trimgalore for ${sample_id}"

# Run trimgalore on both reads

cd "${TMPDIR}"

trim_galore "${full_path_fastq_1}" \
  --cores 2 \
  --gzip \
  --output_dir "${TMPDIR}"

ls ${TMPDIR}

# Change names of validated trimgalore output to basename
mv "${TMPDIR}/${bf1%.*.*}_trimmed.fq.gz" "${TMPDIR}/${bf1}"

loc="${TMPDIR}/${bf1}"

# Calculate length and length standard deviation for salmon
fld=($(zcat ${loc} | awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sqrt(sq/n-m*m));}' - ))

echo "`date` Running salmon quant for ${sample_id}"

for used_index in ${indices[@]}; do
  index_name=$(basename ${used_index} | cut -f 1 -d "_")
  echo "`date` Running salmon quant with ${used_index} for ${sample_id}"
  mkdir -p "${wd}/data/processed/salmon_quant/${index_name}/${sample_id}/"
  # Run salmon for transcript counts
  salmon quant \
    --libType "A" \
    --validateMappings \
    --fldMean ${fld[2]#*=} \
    --fldSD ${fld[3]#*=} \
    --gcBias \
    --quiet \
    --numGibbsSamples 30 \
    --threads ${cpu} \
    -i "${used_index}" \
    -r "${TMPDIR}/${bf1}" \
    --output "${wd}/data/processed/salmon_quant/${index_name}/${sample_id}/"
done

echo "`date` Finished salmon ${sample_id}"
