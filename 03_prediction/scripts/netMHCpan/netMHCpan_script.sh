#!/bin/bash

# $4 = output .xls file -> 1 big table all HLA types after eachother

mkdir "${outdir}/netMHCpan_${run}" || { echo "Error: Failed to create the output directory."; exit 1; }

temp_fasta="${project_data_folder}/temp_${run}.fasta"
name_table="${project_data_folder}/temp_${run}_names.tsv"
output_file="${outdir}/netMHCpan_${run}/netMHCpan"

# Check if output files exist already
if [[ -f "$temp_fasta" || -f "$name_table" ]]; then
    read -p "The output files already exist. Do you want to overwrite them? (y/n): " overwrite_choice

    if [[ $overwrite_choice != "y" && $overwrite_choice != "Y" ]]; then
        echo "Operation aborted. Please choose different names for output files."
        exit 1
    fi
fi

# Create a temp name file and fasta file because NetMHCPan cant handle long input names in the fasta
counter=1

while IFS= read -r line; do
    if [[ $line =~ ^\>(.*)$ ]]; then
        name="${BASH_REMATCH[1]}"
        new_name=$(printf "%015d" $counter)
        echo -e "$name\t$new_name" >> "$name_table"
        echo ">$new_name" >> "$temp_fasta"
        ((counter++))
    else
        echo "$line" >> "$temp_fasta"
    fi
done < "${input_fasta}"

# Run NetMHCpan
apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch,$TMPDIR:$TMPDIR "${apptainer_dir}/netmhcpan-4.1b.sif" \
    /app/package/netMHCpan-4.1/netMHCpan \
        -BA \
        -a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01,HLA-B15:01 \
        -xls -xlsfile "${output_file}_temp.tsv" \
        -f "$temp_fasta" 

# Change the names back to the original long names

# Step 2: Use awk to process the temp_file and replace the spaces with tabs
awk -F"\t" 'FNR==NR { map[$2]=$1; next } FNR==NR { print; next } { if ($3 in map) $3=map[$3]} 1' "$name_table" "${output_file}_temp.tsv" | awk 'BEGIN { OFS="\t" } { gsub(/ /, "\t"); print }' > "${output_file}.tsv"

# print fragment number, peptide, ID and count of binders.
# copy columns of interest with awk (this changes when the amount of HLA changes)
# awk '{print $1,"\t",$2,"\t",$3,"\t",$NF}' $1 > $2

echo -e "Pos\tPeptide\tID\tHLA-A01:01\tHLA-A02:01\tHLA-A03:01\tHLA-A24:02\tHLA-A26:01\tHLA-B07:02\tHLA-B08:01\tHLA-B27:05\tHLA-B39:01\tHLA-B40:01\tHLA-B58:01\tHLA-B15:01" > "${output_file}_summary.tsv"
tail -n +3 "${output_file}.tsv" | awk 'BEGIN { OFS="\t" } {print $1, $2, $3, $7, $13, $19, $25, $31, $37, $43, $49, $55, $61, $67, $73}' >> "${output_file}_summary.tsv"

# Cleanup: Remove the temporary file
rm "$temp_file"
rm "${output_file}_temp.tsv"
rm "$name_table"
rm "$temp_fasta"

# Run the R script to make a overview table with ORF_id column and a column with the weak binding peptides and a column with the strong binding peptides
# Rscript SB_WB_netMHCpan_overview.R "${output_file}_summary.tsv" "${output_file}_SB_WB.tsv"
apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch,$TMPDIR:$TMPDIR,/hpc/local/Rocky8/pmc_vanheesch/Rstudio_Server_Libs/Rstudio_4.3.0_libs:/usr/local/lib/R/site-library "${apptainer_dir}/rstudio_4.3.0_bioconductor.sif" \
    Rscript "${scriptdir}/netMHCpan/SB_WB_netMHCpan_overview.R" \
        "${output_file}_summary.tsv" \
        "${output_file}_SB_WB.tsv"