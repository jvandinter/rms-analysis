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

mkdir -p "${outdir}/price/${pool_id}"

module load gedi/1.0.6
module load java/1.8.0

gedi Nashorn -e \
  'load("'${outdir}/price/${pool_id}.orfs.cit'").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getGeneId()+"__"+o.data.getTranscript()+"__"+o.data.getType()+"__"+o.data.getOrfid()+"__"+o.data.getStartCodon())))).print()' > "${outdir}/price/${pool_id}.orfs.bed"
