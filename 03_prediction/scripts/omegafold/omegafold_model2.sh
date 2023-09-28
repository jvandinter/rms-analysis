#!/bin/bash

# Script that executes OmegaFold and predicts 3D structure with it

mkdir "${outdir}/omegafold_${run}" || { echo "Error: Failed to create the output directory."; exit 1; }

apptainer exec --nv -B /hpc:/hpc "${apptainer_dir}/omegafold-1.0.0.sif" \
    omegafold \
        --weights_file=${model_omegafold} \
        --model 2 \
        "${input_fasta}" \
        "${outdir}/omegafold_${run}/pdb_files" \
        || { echo "Error: Failed to execute OmegaFold."; exit 1; }

# create file with paths to the created pdb files
ls -1 ${outdir}/omegafold_${run}/pdb_files/*.pdb > "${outdir}/omegafold_${run}/overview_pdb_files.txt"

# A simple bash loop to calculate the mean pLDDT score of a predicted OmegaFold structure
# Initialize variables
output_file="${outdir}/omegafold_${run}/parsed_omegafold_${run}.tsv"

# Add the header to the output file
echo -e "protein_id\tmean_pLDDT_OmegaFold\tvalue_per_residue_OmegaFold" > "${output_file}"

# Read each line from the file specified in $2/overview_pdb_files.txt
while IFS= read -r line
do
  # Calculate the average of the 11th column using awk
  average=$(awk '{ total += $11; count++ } END { print total/count }' $line)

  # Extract all the values in the 11th column separated by commas
  values=$(awk '{ printf (NR>1 ? "," : "") $11 }' "$line")
  
  # Extract the protein_id from the line
  protein_id=$(basename $line | sed 's/\.pdb$//')
  
  # Append the protein_id and average to the output file
  echo -e "${protein_id}\t${average}\t${values}" >> "${output_file}" 
done < "${outdir}/omegafold_${run}/overview_pdb_files.txt"