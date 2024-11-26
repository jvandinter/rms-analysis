#!/bin/bash
#SBATCH --job-name calc_CDS
#SBATCH --cpus-per-task 6
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --gres=tmpspace:10G

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -s|--script)
    SCRIPT="$2"
    shift
    shift
    ;;
    "")
    echo "Error: no option provided"
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    exit 1
    ;;
esac
done

version="4.3.0"

# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
singularity_dir="/hpc/local/Rocky8/pmc_vanheesch/singularity_images"
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
  Rscript ${SCRIPT}
