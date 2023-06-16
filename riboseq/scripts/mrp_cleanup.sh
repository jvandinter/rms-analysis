#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 21-06-2021
#
######################################################################

# Load parameters from main script
source $1

# Remove unnecessary files
cd "${wd}/processed"

# Remove trimmed fastq files
for f in */trimgalore/*.fastq.gz; do
    rm $f
done

# Remove bowtie2 files
for f in */bowtie2/*; do
    rm $f
done