#!/bin/bash 


indices=(/hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/NBL_complete /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/RMS_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/EPN_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/ATRT_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/EWS_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/B-ALL_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/T-ALL_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/AML_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/MBL_full /hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/data/processed/salmon_index/WT_full)
data_folder="/hpc/pmc_vanheesch/data"
cpu=12
wd=`pwd`

# Load correct modules
module load salmon/1.8.0
module load cutadapt/3.4
module load fastqc/0.11.9
module load trimgalore/0.6.6

# Get correct files
fastq_names=($(find "${wd}/data/raw" -maxdepth 1 -name "*_R1_001*" -exec basename {} \; 2> >(grep -v 'Permission denied' >&2) | sort -u ))

for i in ${!fastq_names[@]}; do

    fastq_files[i]=$(find "${data_folder}" -name ${fastq_names[i]} 2> >(grep -v 'Permission denied' >&2))
    sample_ids[i]=$(basename ${fastq_files[i]} | cut -f 1 -d "_")

done

full_path_fastq_1="${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]}"
full_path_fastq_2=$(echo $full_path_fastq_1 | sed 's/_R1_001/_R2_001/')

sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

if [[ -f "${wd}/data/processed/salmon_quant/${index_name}/${sample_id}/quant.sf" ]]; then
    echo "`date` Salmon count files for ${index_name}-${sample_id} already present"
    exit 0
fi

echo "`date` running ${sample_id} for ${index_name}"

bf1=$(basename ${full_path_fastq_1})
bf2=$(basename ${full_path_fastq_2})

# Create output dirs


echo "`date` running trimgalore for ${sample_id}"

# Run trimgalore on both reads

cd "${TMPDIR}"

trim_galore "${full_path_fastq_1}" "${full_path_fastq_2}" \
  --cores 2 \
  --paired \
  --gzip \
  --output_dir "${TMPDIR}"

# Change names of validated trimgalore output to basename
mv "${TMPDIR}/${bf1%.*.*}_val_1.fq.gz" "${TMPDIR}/${bf1}"
mv "${TMPDIR}/${bf2%.*.*}_val_2.fq.gz" "${TMPDIR}/${bf2}"

echo "`date` Running salmon quant for ${sample_id}"

for used_index in ${indices[@]}; do
  index_name=$(basename ${used_index} | cut -f 1 -d "_")
  echo "`date` Running salmon quant with ${used_index} for ${sample_id}"
  mkdir -p "${wd}/data/processed/salmon_quant/${index_name}/${sample_id}/"
  # Run salmon for transcript counts
  salmon quant \
    --libType "A" \
    --validateMappings \
    --gcBias \
    --quiet \
    --numGibbsSamples 30 \
    --threads ${cpu} \
    -i "${used_index}" \
    -1 "${TMPDIR}/${bf1}" \
    -2 "${TMPDIR}/${bf2}" \
    --output "${wd}/data/processed/salmon_quant/${index_name}/${sample_id}/"
done

echo "`date` Finished salmon ${sample_id}"
