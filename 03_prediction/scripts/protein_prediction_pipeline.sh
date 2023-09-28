#!/bin/bash

######################################################################
# 
# Protein prediction project pipeline
# Pipeline that predicts several characteristics of a microprotein and
# puts these results in a .tsv file which can be examined. Currenly the 
# pipeline excists of 4 tools; OmegaFold, BLASTp, DeepTMHMM and SignalP. 
#
# Resources 
# Resources need to be alted depending on the size of the input fasta
# the BLASTp, deepTMHMM, SignalP 6.0 and merging of files don't need a lot of resources
# so this can stay the same. For OmegaFold on a GPU all 8 microproteins 
# can be ran in under a minute but need 5Gb of memory -> adjust accordingly.
#
# List of output files per tool:
# OmegaFold:    .pdb file per microprotein sequence
# BLASTp:       all hits found with BLASTp
# DeepTMHMM:    3 files containing info about whether there is a trans 
#               membrane domain present. 
# SignalP 6.0:  .txt files per protein for a plot. For every amino acid
#               a score is predicted whether there is a signal present.
#               And a overview file of all the proteins and signals.
# IUPred 3.0:   -
# Characteristics:  -
# NetMHCpan:    -
# ELM search:   -
# 
# Author:   Amalia Nabuurs (A.J.M.Nabuurs-3@prinsesmaximacentrum.nl)
# Date:     26-06-2023
#
######################################################################

set -euo pipefail

function usage() {
    cat <<EOF
SYNOPSIS
  protein_prediction_pipeline.sh [-c <config file>] [-h]
DESCRIPTION
  Run the protein prediction pipeline consisting of the following steps:
  1. Predict 3D structure with OmegaFold
  2. Check for homology with BLASTp
  3. Look for transmembrane domains with DeepTMHMM
OPTIONS
  -c, --config <file>    Configuration file to use
  -h, --help             Display this help message
AUTHOR
  Amalia Nabuurs
EOF
}

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c|--config)
    CONFIG="$2"
    shift
    shift
    ;;
    -h|--help)
    usage
    exit
    ;;
    "")
    echo "Error: no option provided"
    usage
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    usage
    exit 1
    ;;
esac
done

# Check that configuration file is provided
if [[ -z ${CONFIG+x} ]]; then 
    echo "Error: no configuration file provided"
    usage
    exit 1
fi

# Load configuration variables
source $CONFIG

# Create a unique prefix for the names for this run of the pipeline
# This makes sure that runs can be identified
export run=$(uuidgen | tr '-' ' ' | awk '{print $1}')

######################################################################

# Step 1: running of several tools
# Run OmegaFold

mkdir -p ${wd}/log/omegafold
omegafold_jobid=$(sbatch --parsable \
    -p gpu \
    -c 2 \
    --gpus-per-node=1 \
    --mem=80G \
    --time=10:00:00 \
    --job-name=${run}.omegafold \
    --output=${wd}/log/omegafold/%A.out \
    --export=ALL \
    "${scriptdir}/omegafold/omegafold_model2.sh")
echo "OmegaFold jobid: ${omegafold_jobid}"

# Run BLASTp

#mkdir -p ${wd}/log/BLASTp
#BLASTp_jobid=$(sbatch --parsable \
#    --mem=30G \
#    --time=6:00:00 \
#    --job-name=${run}.BLASTp \
#    --output=${wd}/log/BLASTp/%A.out \
#    --export=ALL \
#    "${scriptdir}/BLASTp/BLASTp_remote.sh")
#echo "BLASTp jobid: ${BLASTp_jobid}"

# Run DeepTMHMM

mkdir -p ${wd}/log/deepTMHMM
deepTMHMM_jobid=$(sbatch --parsable \
    --mem=10G \
    --time=4:00:00 \
    --job-name=${run}.deepTMHMM \
    --output=${wd}/log/deepTMHMM/%A.out \
    --export=ALL \
    "${scriptdir}/deepTMHMM/deepTMHMM.sh")
echo "DeepTMHMM jobid: ${deepTMHMM_jobid}"

# Run SignalP 6.0

mkdir -p ${wd}/log/signalP6
signalp6_jobid=$(sbatch --parsable \
    --mem=10G \
    --time=2:00:00 \
    --job-name=${run}.signalP6 \
    --output=${wd}/log/signalP6/%A.out \
    --export=ALL \
    "${scriptdir}/signalp6_fast/signalp6_fast.sh")
echo "signalP6 jobid: ${signalp6_jobid}"

# Run IUPred3

mkdir -p ${wd}/log/IUPred3
iupred3_jobid=$(sbatch --parsable \
    --mem=10G \
    --time=3:00:00 \
    --job-name=${run}.IUPred3 \
    --gres=tmpspace:10G \
    --output=${wd}/log/IUPred3/%A.out \
    --export=ALL \
    "${scriptdir}/iupred3/iupred3.sh")
echo "IUPred3 jobid: ${iupred3_jobid}"

# Run characteristics R script

mkdir -p ${wd}/log/characteristics
characteristics_jobid=$(sbatch --parsable \
    --mem=30G \
    --time=2:00:00 \
    --job-name=${run}.characteristics \
    --gres=tmpspace:30G \
    --output=${wd}/log/characteristics/%A.out \
    --export=ALL \
    "${scriptdir}/characteristics/calculate_characteristics.sh")
echo "characteristics jobid: ${characteristics_jobid}"

# Run netMHCpan
mkdir -p ${wd}/log/netMHCpan
netMHCpan_jobid=$(sbatch --parsable \
    --mem=10G \
    --time=18:00:00 \
    --job-name=${run}.netMHCpan \
    --gres=tmpspace:10G \
    --output=${wd}/log/netMHCpan/%A.out \
    --export=ALL \
    "${scriptdir}/netMHCpan/netMHCpan_script.sh")
echo "netMHCpan jobid: ${netMHCpan_jobid}"

# Run ELM search after IUPred3 is finished
mkdir -p ${wd}/log/ELM_search
ELM_search_jobid=$(sbatch --dependency=afterok:${iupred3_jobid} \
    --parsable \
    --mem=10G \
    --time=2:00:00 \
    --job-name=${run}.ELM_search \
    --gres=tmpspace:10G \
    --output=${wd}/log/ELM_search/%A.out \
    --export=ALL \
    "${scriptdir}/ELM_search/ELM_SLiM_search.sh")
echo "ELM search jobid: ${ELM_search_jobid}"

# Step 2: combine the results

sbatch --dependency=afterok:${omegafold_jobid},${deepTMHMM_jobid},${signalp6_jobid},${iupred3_jobid},${characteristics_jobid},${netMHCpan_jobid} \
    --mem=10G \
    --time=2:00:00 \
    --job-name=${run}.merge_files \
    --output=${wd}/log/%A.out \
    --export=ALL \
    "${scriptdir}/merge_files.sh"

# ,${BLASTp_jobid}
