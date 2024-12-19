#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --job-name="netMHCpan"
#SBATCH --mem=8G 
#SBATCH --time=24:00:00 
#SBATCH --gres=tmpspace:50G 
#SBATCH --array=1-60%5

set -euox pipefail

outdir=/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis/03_prediction
apptainer_dir="/hpc/local/Rocky8/pmc_vanheesch/singularity_images"
processed_fasta="${outdir}/data/RMS_enriched.fasta"
netmhcpan_array=(8-HLA-A01:01 9-HLA-A01:01 10-HLA-A01:01 11-HLA-A01:01 12-HLA-A01:01 8-HLA-A02:01 9-HLA-A02:01 10-HLA-A02:01 11-HLA-A02:01 12-HLA-A02:01 8-HLA-A03:01 9-HLA-A03:01 10-HLA-A03:01 11-HLA-A03:01 12-HLA-A03:01 8-HLA-A24:02 9-HLA-A24:02 10-HLA-A24:02 11-HLA-A24:02 12-HLA-A24:02 8-HLA-A26:01 9-HLA-A26:01 10-HLA-A26:01 11-HLA-A26:01 12-HLA-A26:01 8-HLA-B07:02 9-HLA-B07:02 10-HLA-B07:02 11-HLA-B07:02 12-HLA-B07:02 8-HLA-B08:01 9-HLA-B08:01 10-HLA-B08:01 11-HLA-B08:01 12-HLA-B08:01 8-HLA-B27:05 9-HLA-B27:05 10-HLA-B27:05 11-HLA-B27:05 12-HLA-B27:05 8-HLA-B39:01 9-HLA-B39:01 10-HLA-B39:01 11-HLA-B39:01 12-HLA-B39:01 8-HLA-B40:01 9-HLA-B40:01 10-HLA-B40:01 11-HLA-B40:01 12-HLA-B40:01 8-HLA-B58:01 9-HLA-B58:01 10-HLA-B58:01 11-HLA-B58:01 12-HLA-B58:01 8-HLA-B15:01 9-HLA-B15:01 10-HLA-B15:01 11-HLA-B15:01 12-HLA-B15:01)
# Create a temp name file and fasta file because NetMHCPan cant handle long input names in the fasta
item="${netmhcpan_array[$((SLURM_ARRAY_TASK_ID-1))]}"
length="${item%%-*}"  # Everything before the first '-'
hla_type="${item#*-}" # Everything after the first '-'

output_file="${outdir}/analysis/netMHCpan/RMS"

# Run NetMHCpan
apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch,$TMPDIR:$TMPDIR "${apptainer_dir}/netmhcpan-4.1b.sif" \
    /app/package/netMHCpan-4.1/netMHCpan \
        -BA \
        -l ${length} \
        -a ${hla_type} \
        -xls \
        -xlsfile "${output_file}_${hla_type}_${length}_temp.tsv" \
        -f "$processed_fasta"
