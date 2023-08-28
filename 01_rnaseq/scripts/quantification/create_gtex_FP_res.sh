#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --job-name gtex_results
#SBATCH --gres=tmpspace:20G
#SBATCH --mem=144G
#SBATCH --time=144:00:00


# Run this script with 'sbatch /hpc/local/Rocky8/pmc_vanheesch/small_tools/rstudio_server.sh'
#
# This script launches an RStudio Server Apptainer image, which can then be accessed from a local
# workstation by following the instructions that are printed to the user's home directory.
# This script was adapted from the Rocker project (https://www.rocker-project.org/use/singularity).
#
# The apptainer image is stored in '/hpc/local/Rocky8/pmc_vanheesch/singularity_images/'; this
# script can be easily adjusted to work with images of newer R versions.
# To download a new apptainer image, simply start an interactive job with 'srun' (apptainer
# is only available on compute nodes), cd to the apptainer images directory, and pull the image
# with 'apptainer pull docker://bioconductor/bioconductor_docker:devel', replacing the version number as needed.
#
# Be consistent with naming the apptainer image use 'rstudio_${version}_bioconductor.sif'
# Create a directory for the version of Rstudio you want to use in
# '/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/' and use 'Rstudio_${version}_libs'
# and point the version="...." below to the version of R you used for naming of the image and directory.
#
# Damon Hofman
# 22-06-2022
#
# conversion to Rocky8 and using seperate directories for packages in different R versions
# 06-06-2023 by Amalia Nabuurs

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

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0

apptainer exec --cleanenv ${singularity_dir}/rstudio_${version}_bioconductor.sif \
  Rscript /hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/scripts/quantification/create_gtex_FP_res.R
