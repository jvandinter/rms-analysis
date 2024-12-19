#!/bin/bash

mkdir -p "${outdir}/iupred3_${run}/per_protein"

# Read the FASTA file line by line
while IFS= read -r line
do
  # Check if the line starts with '>'
  if [[ $line == ">"* ]]; then
    # Extract the protein sequence ID
    sequence_id="${line:1}"
    echo "Submitting sequence: $sequence_id"
  else
    # Create a temporary file with the protein sequence
    tmp_file=$(mktemp)
    echo "$line" > "$tmp_file"

    apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch,$TMPDIR:$TMPDIR ${apptainer_dir}/iupred-3.0.sif\
        python3 /app/iupred3.py "$tmp_file" short > "${outdir}/iupred3_${run}/per_protein/${sequence_id}.txt"
    
    rm "$tmp_file"
  fi
done < "${input_fasta}"


# create file with paths to the created pdb files
ls -1 ${outdir}/iupred3_${run}/per_protein/*.txt > "${outdir}/iupred3_${run}/overview_iupred3_files.tsv"

# A simple bash loop to calculate the mean pLDDT score of a predicted OmegaFold structure
# Initialize variables
output_file="${outdir}/iupred3_${run}/parsed_iupred3_${run}.tsv"

# Add the header to the output file
echo -e "protein_id\tmean_iupred3_score\tvalue_per_residue_iupred3" > "${output_file}"

# Read each line from the file specified in $2/overview_pdb_files.txt
while IFS= read -r line
do
  # Calculate the average of the 11th column using awk
  average=$(awk 'NR > 12 { total += $3; count++ } END { print total/count }' $line)

  # Extract all the values in the 3rd column separated by commas
  values=$(awk 'NR > 12 { printf (NR>13 ? "," : "") $3 }' "$line")
  
  # Extract the protein_id from the line
  protein_id=$(basename $line | sed 's/\.txt$//')
  
  # Append the protein_id and average to the output file
  echo -e "${protein_id}\t${average}\t${values}" >> "${output_file}" 
done < "${outdir}/iupred3_${run}/overview_iupred3_files.tsv"