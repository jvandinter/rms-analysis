#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --job-name deseq2_results
#SBATCH --gres=tmpspace:20G
#SBATCH --mem=120G
#SBATCH --time=24:00:00

version="4.3.0"

# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
singularity_dir="/hpc/local/Rocky8/pmc_vanheesch/singularity_images"
wd="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq"
workdir=${TMPDIR}

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment

cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=${SLURM_JOB_CPUS_PER_NODE}
export R_LIBS_USER=/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_${version}_libs
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export APPTAINER_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server,/hpc/pmc_vanheesch,/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_${version}_libs:/usr/local/lib/R/site-library"

apptainer exec --cleanenv ${singularity_dir}/rstudio_${version}_bioconductor.sif \
  Rscript /hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/scripts/quantification/create_deseq2_results.R \
