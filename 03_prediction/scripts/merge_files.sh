#!/bin/sh
#SBATCH --time=00:10:00           
#SBATCH --mem=1G     
#SBATCH --output=/hpc/pmc_vanheesch/projects/Amalia/20230503_AN_peptide_prediction_project/20230623_AN_pipeline/log/slurm_%j.out

# script to merge several tsv files with overlapping column to a table

# Run the command with the input file arguments
apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch /hpc/local/Rocky8/pmc_vanheesch/singularity_images/merge_tsvfiles_python.sif \
    python3 /app/dataframe_merge.py \
        -i "${outdir}/characteristics_${run}/characteristics_proteins.tsv" "${outdir}/omegafold_${run}/parsed_omegafold_${run}.tsv" "${outdir}/deepTMHMM_${run}/parsed_deepTMHMM_${run}.tsv" "${outdir}/signalp_${run}/parsed_signalp_${run}.tsv" "${outdir}/iupred3_${run}/parsed_iupred3_${run}.tsv" "${outdir}/netMHCpan_${run}/netMHCpan_SB_WB.tsv" "${outdir}/ELM_search_${run}/ELM_SLiMs_overview.tsv" \
        -o "${outdir}/run_${run}_table.tsv"

# "${outdir}/BLASTp_${run}/parsed_BLASTp_${run}.tsv"