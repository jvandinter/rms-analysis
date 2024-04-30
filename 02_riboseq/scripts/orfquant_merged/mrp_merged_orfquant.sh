#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 04-06-2021
#
######################################################################

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

# Check whether script needs to run
if [[ -f "${outdir}/ORFquant/${pool_id}/${pool_id}.html" ]]; then
  echo "File already present"
  exit 0
fi

echo "`date` calling ORFs for ${pool_id}"

# Create output dirs
cd "${outdir}/"
mkdir -p "ORFquant/${pool_id}/"

apptainer exec -B /hpc:/hpc,${TMPDIR}:${TMPDIR} --env LC_CTYPE="en_US.UTF-8" ${container_dir}/orfquant-4.1.2b.sif Rscript "${scriptdir}/mrp_orfquant.R" \
  ${outdir} \
  "${outdir}/merged_p_sites/for_ORFquant_merged" \
  "${pool_id}" \
  "${annot_name}" \
  "${rannot}" \
  "${threads}" \
  "${pandoc_dir}" \
  "${resource_dir}" \
  "${annotation_package}"  \
  ${package_install_loc}

  echo "`date` finished ${pool_id}"
