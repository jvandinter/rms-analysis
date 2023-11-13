#!/bin/bash

set -euo pipefail

# Load configuration variables
export wd="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/netmhcpan_single"
export project_data_folder="${wd}/data"
export outdir="${wd}/analysis"
export scriptdir="${wd}/scripts"
export input_fasta="${project_data_folder}/rms_neoantigen_strict.fasta"
export apptainer_dir="/hpc/local/Rocky8/pmc_vanheesch/singularity_images"
export protdir="/hpc/pmc_vanheesch/shared_resources/proteome"

# Create a unique prefix for the names for this run of the pipeline
# This makes sure that runs can be identified
export run="RMS_neoantigen_strict"

######################################################################

# Run netMHCpan

mkdir -p ${wd}/log/${run}/netMHCpan
# input_netMHCpan_jobid=$(sbatch --parsable \
#     --mem=10G \
#     --time=18:00:00 \
#     --job-name=${run}.input_netMHCpan \
#     --gres=tmpspace:10G \
#     --output=${wd}/log/${run}/netMHCpan/%A_input.out \
#     --export=ALL \
#     "${scriptdir}/app_prepare_input.sh")
# echo "netMHCpan input jobid: ${input_netMHCpan_jobid}"

# process_fasta_jobid=$(sbatch --parsable \
#     --mem=10G \
#     --time=18:00:00 \
#     --job-name=${run}.process.fasta \
#     --gres=tmpspace:10G \
#     --output=${wd}/log/${run}/netMHCpan/%A_process.out \
#     --dependency=afterok:${input_netMHCpan_jobid} \
#     --export=ALL \
#     "${scriptdir}/create_temp_fasta.sh")
# echo "tempfasta parser jobid: ${process_fasta_jobid}"

# netMHCpan_jobid=$(sbatch --parsable \
#     --mem=10G \
#     --time=72:00:00 \
#     --job-name=${run}.netMHCpan \
#     --gres=tmpspace:10G \
#     --output=${wd}/log/${run}/netMHCpan/%A_%a.out \
#     --dependency=afterok:${process_fasta_jobid} \
#     --export=ALL \
#     --array=1-60%15 \
#     "${scriptdir}/netMHCpan_script.sh")
# echo "netMHCpan jobid: ${netMHCpan_jobid}"

netMHCpan_jobid=$(sbatch --parsable \
    --mem=10G \
    --time=72:00:00 \
    --job-name=${run}.netMHCpan \
    --gres=tmpspace:10G \
    --output=${wd}/log/${run}/netMHCpan/%A_%a.out \
    --export=ALL \
    --array=1-60%15 \
    "${scriptdir}/netMHCpan_script.sh")
echo "netMHCpan jobid: ${netMHCpan_jobid}"
