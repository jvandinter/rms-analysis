#!/bin/bash

apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch ${apptainer_dir}/signalp-6.0-fast.sif \
    signalp6 \
        --fastafile ${input_fasta} \
        --output_dir "${outdir}/signalp_${run}" \
        --format txt \
        --organism eukarya \
        --mode fast

echo -e "protein_id\tsignalP_result" > "${outdir}/signalp_${run}/parsed_signalp_${run}.tsv"
sed 's/#//' "${outdir}/signalp_${run}/prediction_results.txt" | tail -n +3 | awk '{print $1"\t"$2}' >> "${outdir}/signalp_${run}/parsed_signalp_${run}.tsv"