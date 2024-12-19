#!/bin/bash 

#SBATCH --cpus-per-task 8
#SBATCH --job-name chipseq
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --gres=tmpspace:10G
#SBATCH --array=7

# Load modules
module load deepTools/3.5.1
module load bwa/0.7.17
module load trimgalore/0.6.6
module load cutadapt/3.4
module load fastqc/0.11.9
module load picardtools/2.21.1
module load samtools/1.12

wd=`pwd`
bwa_ref="/hpc/pmc_vanheesch/shared_resources/GENOMES/Homo_sapiens.GRCh38/102/bwa/Homo_sapiens.GRCh38"
data_folder="/hpc/pmc_vanheesch/data/external/dnaseq"
threads=$(( ${SLURM_CPUS_PER_TASK} * 2 ))
echo $threads

fastq_names=($(find "${wd}/data/raw" -maxdepth 1 -exec basename {} \; 2> >(grep -v 'Permission denied' >&2) | sort -u ))

for i in ${!fastq_names[@]}; do

  fastq_files[i]=$(find "${data_folder}" -name ${fastq_names[i]} 2> >(grep -v 'Permission denied' >&2))

done

fastq_name=${fastq_names[$((SLURM_ARRAY_TASK_ID-1))]]}
sample_id=${fastq_name%".fastq.gz"}

mkdir -p ${wd}/data/processed/{picard,bwa-mem,trimgalore}/${sample_id}

# # trimgalore
# cutadapt --version
# fastqc --version

# trim_galore ${fastq_files[$((SLURM_ARRAY_TASK_ID-1))]]} \
#   --cores 4 \
#   --gzip \
#   --length 25 \
#   --trim-n \
#   --fastqc \
#   --fastqc_args "--outdir ${wd}/data/processed/trimgalore/${sample_id}/" \
#   --output_dir "$wd/data/processed/trimgalore/${sample_id}/"

# bwa-mem
bwa mem \
  -t ${threads} \
  -M \
  ${bwa_ref} \
  ${wd}/data/processed/trimgalore/${sample_id}/${sample_id}_trimmed.fq.gz | \
  samtools sort -@ 2 -o ${wd}/data/processed/bwa-mem/${sample_id}/${sample_id}.bam -

# Picard MarkDuplicates
java -Xmx24G -Djava.io.tmpdir=$TMPDIR -jar $PICARD MarkDuplicates \
      I=${wd}/data/processed/bwa-mem/${sample_id}/${sample_id}.bam \
      O=${wd}/data/processed/picard/${sample_id}/${sample_id}.bam \
      M=${wd}/data/processed/picard/${sample_id}/${sample_id}_marked_dup_metrics.txt \
      VALIDATION_STRINGENCY=SILENT

samtools index -@ ${threads} ${wd}/data/processed/picard/${sample_id}/${sample_id}.bam

# Remove duplicates
samtools view -@ ${threads} -b -F 0x400 ${wd}/data/processed/picard/${sample_id}/${sample_id}.bam > ${wd}/data/processed/picard/${sample_id}/${sample_id}_dedup.bam
samtools index -@ ${threads} ${wd}/data/processed/picard/${sample_id}/${sample_id}_dedup.bam

# Subset on chr2
samtools view -@ ${threads} -bSh ${wd}/data/processed/picard/${sample_id}/${sample_id}_dedup.bam 2 > ${wd}/data/processed/picard/${sample_id}/${sample_id}_chr2_dedup.bam
