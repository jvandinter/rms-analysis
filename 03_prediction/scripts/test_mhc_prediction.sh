#!/bin/bash

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

# Run netMHCpan
input_netMHCpan_jobid=$(sbatch --parsable \
    --mem=10G \
    --time=24:00:00 \
    --job-name=${run}.input_netMHCpan \
    --gres=tmpspace:10G \
    --output=${wd}/log/netMHCpan/%A_input.out \
    --export=ALL \
    "${scriptdir}/netMHCpan/app_prepare_input.sh")
echo "netMHCpan jobid: ${input_netMHCpan_jobid}"

mkdir -p ${wd}/log/netMHCpan
netMHCpan_jobid=$(sbatch --parsable \
    --mem=10G \
    --time=24:00:00 \
    --job-name=${run}.netMHCpan \
    --gres=tmpspace:10G \
    --output=${wd}/log/netMHCpan/%A.out \
    --dependency=afterok:${input_netMHCpan_jobid} \
    --export=ALL \
    "${scriptdir}/netMHCpan/netMHCpan_script.sh")
echo "netMHCpan jobid: ${netMHCpan_jobid}"
