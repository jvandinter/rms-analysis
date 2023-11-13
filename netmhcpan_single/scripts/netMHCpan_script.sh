#!/bin/bash

# $4 = output .xls file -> 1 big table all HLA types after eachother

outdir="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/netmhcpan_single/analysis"
container_dir="/hpc/local/Rocky8/pmc_vanheesch/singularity_images/"
processed_fasta="${outdir}/netMHCpan/RMS_neoantigen_strict/peptides.fasta"
length_alleles=(8_HLA-A01:01 8_HLA-A02:01 8_HLA-A03:01 8_HLA-A24:02 8_HLA-A26:01 8_HLA-B07:02 8_HLA-B08:01 8_HLA-B27:05 8_HLA-B39:01 8_HLA-B40:01 8_HLA-B58:01 8_HLA-B15:01 9_HLA-A01:01 9_HLA-A02:01 9_HLA-A03:01 9_HLA-A24:02 9_HLA-A26:01 9_HLA-B07:02 9_HLA-B08:01 9_HLA-B27:05 9_HLA-B39:01 9_HLA-B40:01 9_HLA-B58:01 9_HLA-B15:01 10_HLA-A01:01 10_HLA-A02:01 10_HLA-A03:01 10_HLA-A24:02 10_HLA-A26:01 10_HLA-B07:02 10_HLA-B08:01 10_HLA-B27:05 10_HLA-B39:01 10_HLA-B40:01 10_HLA-B58:01 10_HLA-B15:01 11_HLA-A01:01 11_HLA-A02:01 11_HLA-A03:01 11_HLA-A24:02 11_HLA-A26:01 11_HLA-B07:02 11_HLA-B08:01 11_HLA-B27:05 11_HLA-B39:01 11_HLA-B40:01 11_HLA-B58:01 11_HLA-B15:01 12_HLA-A01:01 12_HLA-A02:01 12_HLA-A03:01 12_HLA-A24:02 12_HLA-A26:01 12_HLA-B07:02 12_HLA-B08:01 12_HLA-B27:05 12_HLA-B39:01 12_HLA-B40:01 12_HLA-B58:01 12_HLA-B15:01)

picked_length_allele="${length_alleles[$((SLURM_ARRAY_TASK_ID-1))]}"
length=${picked_length_allele%%_*}
allele=${picked_length_allele#*_}
output_file_prefix="${outdir}/netMHCpan/RMS_neoantigen_strict/netMHCpan_${picked_length_allele}"

temp_fasta="${outdir}/netMHCpan/RMS_neoantigen_strict/temp_${length}.fasta"
name_table="${outdir}/netMHCpan/RMS_neoantigen_strict/temp_names.tsv"

# Run NetMHCpan
apptainer exec -B /hpc/:/hpc/,$TMPDIR:$TMPDIR "${container_dir}/netmhcpan-4.1b.sif" \
    /app/package/netMHCpan-4.1/netMHCpan \
        -l $length \
        -BA \
        -a ${allele} \
        -xls -xlsfile "${output_file_prefix}_temp.tsv" \
        -f "$temp_fasta"

# Change the names back to the original long names

# Step 2: Use awk to process the temp_file and replace the spaces with tabs
awk -F"\t" 'FNR==NR { map[$2]=$1; next } FNR==NR { print; next } { if ($3 in map) $3=map[$3]} 1' "$name_table" "${output_file_prefix}_temp.tsv" | awk 'BEGIN { OFS="\t" } { gsub(/ /, "\t"); print }' > "${output_file_prefix}.tsv"

# print fragment number, peptide, ID and count of binders.
# copy columns of interest with awk (this changes when the amount of HLA changes)
# awk '{print $1,"\t",$2,"\t",$3,"\t",$NF}' $1 > $2

echo -e "Pos\tPeptide\tID\t${alleles[$((SLURM_ARRAY_TASK_ID-1))]}" > "${output_file_prefix}_summary.tsv"
tail -n +3 "${output_file_prefix}.tsv" | awk 'BEGIN { OFS="\t" } {print $1, $2, $3}' >> "${output_file_prefix}_summary.tsv"

# # Run the R script to make a overview table with ORF_id column and a column with the weak binding peptides and a column with the strong binding peptides
# # Rscript SB_WB_netMHCpan_overview.R "${output_file_prefix}_summary.tsv" "${output_file_prefix}_SB_WB.tsv"
# apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch,$TMPDIR:$TMPDIR,/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_4.3.0_libs:/usr/local/lib/R/site-library "${apptainer_dir}/rstudio_4.3.0_bioconductor.sif" \
#     Rscript "${scriptdir}/netMHCpan/SB_WB_netMHCpan_overview.R" \
#         "${output_file_prefix}_summary.tsv" \
#         "${output_file_prefix}_SB_WB.tsv"
