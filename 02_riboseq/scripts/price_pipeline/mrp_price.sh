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

gedi -e Price \
    -reads "${outdir}/${pool_id}/${pool_id}_filtered.bam" \
    -genomic "${outdir}/${pool_id}_index/${pool_id}.oml" \
    -prefix "${outdir}/price/${pool_id}" \
    -plot

#gedi Nashorn -e 'load("'<price_orfs.cit>'").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getGeneId()+"__"+o.data.getTranscript()+"__"+o.data.getType()+"__"+o.data.getOrfid()+"__"+o.data.getStartCodon())))).print()' > "<price_orfs.cit.bed>"
