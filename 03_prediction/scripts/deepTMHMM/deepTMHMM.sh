#!/bin/bash  

# script that executes deepTMHMM and parses results to tsv format
# it is NOT runned locally, we can change this with machin = 'local' in the python script if we want to

mkdir "${outdir}/deepTMHMM_${run}" || { echo "Error: Failed to create the output directory."; exit 1; }

apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch "${apptainer_dir}/deepTMHMM-1.0.24.sif" \
    python3 /app/deeptmhmm.py \
        -i "${input_fasta}" \
        -o "${outdir}/deepTMHMM_${run}" \
        || { echo "Error: Failed to execute deepTMHMM."; exit 1; }

# parse the results from deepTMHMM to format we need/want

apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch "${apptainer_dir}/parser_deepTMHMM-1.0.sif" \
    python3 /app/results_parser.py \
        -i "${outdir}/deepTMHMM_${run}/predicted_topologies.3line" \
        -o "${outdir}/deepTMHMM_${run}/parsed_deepTMHMM_${run}.tsv" \
        || { echo "Error: Failed to parse deepTMHMM results."; exit 1; }