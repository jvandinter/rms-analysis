#!/bin/bash 
#SBATCH --cpus-per-task 1
#SBATCH --job-name="custom_annot"
#SBATCH --mem=10G 
#SBATCH --gres=tmpspace:50G
#SBATCH --time=24:00:00

srun -c 1 --mem 10G --gres=tmpspace:50G --time=24:00:00 --pty bash

scriptdir="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/scripts/container"

outdir="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/analysis/rnaseq_pipeline"

merged_gtf_basename="RMS_container"

container_dir="/hpc/local/Rocky8/pmc_vanheesch/singularity_images"

twobit="/hpc/pmc_vanheesch/shared_resources/GENOMES/Homo_sapiens.GRCh38/102/Homo_sapiens.GRCh38.dna.primary_assembly.2bit"

package_install_loc="/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_4.1.2_libs"

custom_gtf="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/analysis/rnaseq_pipeline/customannotation/RMS_full_novel_filtered_corrected.gtf"

echo "`date` Generating custom annotation ..."

# Prepare RiboseQC and ORFquant annotation files based on 
# merged annotated stringtie GTF
apptainer exec -B "/hpc:/hpc",${TMPDIR}:${TMPDIR} ${container_dir}/orfquant-4.1.2.sif \
  Rscript "${scriptdir}/prepare_custom_annotation.R" \
  ${twobit} \
 ${custom_gtf} \
  "${outdir}/customannotation/${merged_gtf_basename}/" \
  ${merged_gtf_basename} \
  ${package_install_loc}

echo "`date` Generating custom annotation finished"
