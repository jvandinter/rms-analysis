#!/bin/bash
#SBATCH --job-name=STAR-fusion
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --array=1-92%23

# Note: this script doesn't work with symlinks to the fastq files, so you need to directly link to the fastq files.

#fastq_dir="/hpc/pmc_vanheesch/data/maxima/rnaseq/tot_rna_depleted/RMS"
fastq_dir="/hpc/pmc_vanheesch/data/external/rnaseq/20220527_StJudeCloud/RMS/"
genome_lib="/hpc/pmc_gen/references/RNA-Seq/starfusion/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir"
wd=`pwd`
# Get sample IDs
SAMPLE_FILES=($(find "${fastq_dir}" -name '*R1.fastq.gz' -exec basename {} \; | sort | uniq))

# Get the sample name for this array index
sample=${SAMPLE_FILES[$SLURM_ARRAY_TASK_ID-1]}
sample_2=${sample/R1/R2}
#sample_id=$(echo ${sample} | rev | cut -d '_' -f 3- | rev)
sample_id=$(echo ${sample} | cut -d '_' -f 1)
outdir="${wd}/analysis/starfusion/${sample_id}"

mkdir -p ${outdir}

# Run STAR-fusion
singularity exec -e -B /hpc:/hpc \
    "/hpc/local/Rocky8/pmc_vanheesch/singularity_images/starfusion-1.8.0.sif" \
    STAR-Fusion \
    --left_fq "${fastq_dir}/${sample}" \
    --right_fq "${fastq_dir}/${sample_2}" \
    --genome_lib_dir ${genome_lib} \
    --output_dir ${outdir} \
    --CPU 32 \
    --FusionInspector validate
