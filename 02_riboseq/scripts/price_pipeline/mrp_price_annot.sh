#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 21-12-2023
#
######################################################################

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

module load gedi/${gedi_version}

# Check whether script needs to run
if [[ -s "${outdir}/${pool_id}_index/${pool_id}.oml" ]]; then
  echo "`date` PRICE index already present"
  exit 0
fi

gedi -e IndexGenome -s "${reference_genome}" -a "$gtf" \
    -n "${pool_id}" \
    -f "${outdir}/${pool_id}_index" \
	-o "${outdir}/${pool_id}_index/${pool_id}.oml" \
    -nobowtie \
    -nostar \
    -nokallisto
