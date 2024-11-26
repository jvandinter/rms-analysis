#!/bin/bash

#SBATCH --job-name=bamCompareArray
#SBATCH --output=log/bamcompare/bamCompare_%A_%a.out
#SBATCH --array=1-36%36
#SBATCH --cpus-per-task=4
#SBATCH --gres=tmpspace:10G
#SBATCH --mem=20G
#SBATCH --time=02:00:00

set -euo pipefail

output_dir="/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/05_chipseq/analysis/bamcoverage"
cd ${TMPDIR}

# Create an array of antibody BAM files
ANTIBODY_BAMS=($(find /hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/05_chipseq/analysis/bwa -type f -name '*ChIP*.bam' ! -name '*-control.bam' ))

# Get the antibody BAM file based on the SLURM array task ID
ANTIBODY_BAM=${ANTIBODY_BAMS[$SLURM_ARRAY_TASK_ID - 1]}

# Extract the cellline name from the filename
CELLLINE=$(echo "$ANTIBODY_BAM" | cut -d'-' -f2)

# Using find to search for matching files
matching_files=($(find /hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/05_chipseq/analysis/bwa -type f -name "*${CELLLINE}*control.bam"))

# Checking the number of matches
num_matches=${#matching_files[@]}

# Checking if there is exactly one match
if [ $num_matches -eq 1 ]; then
    # Storing the single file in a variable
    INPUT_CONTROL=${matching_files[0]}
else
    INPUT_CONTROL=${matching_files}
fi

# Output file name
OUTPUT_BW=$(basename "$ANTIBODY_BAM" .bam)_normalized.bw

# Execute bamCompare
apptainer exec -B /hpc:/hpc/,${TMPDIR}:${TMPDIR} /hpc/local/Rocky8/pmc_vanheesch/singularity_images/deeptools_3.5.4--pyhdfd78af_1.sif bamCompare \
  -b1 $ANTIBODY_BAM \
  -b2 $INPUT_CONTROL \
  -o ${output_dir}/$OUTPUT_BW \
  --operation log2 \
  --binSize 20 \
  -p 8
