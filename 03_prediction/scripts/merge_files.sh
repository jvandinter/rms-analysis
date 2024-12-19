#!/bin/sh
#SBATCH --time=00:10:00           
#SBATCH --mem=1G     
#SBATCH --output=/hpc/pmc_vanheesch/projects/Amalia/20230503_AN_peptide_prediction_project/20230623_AN_pipeline/log/slurm_%j.out

# script to merge several tsv files with overlapping column to a table

# Run the command with the input file arguments
apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch /hpc/local/Rocky8/pmc_vanheesch/singularity_images/merge_tsvfiles_python.sif \
    python3 /app/dataframe_merge.py \
        -i "${outdir}/characteristics/characteristics_proteins.tsv" "${outdir}/omegafold/parsed_omegafold.tsv" "${outdir}/deepTMHMM/parsed_deepTMHMM.tsv" "${outdir}/signalp/parsed_signalp.tsv" "${outdir}/iupred3/parsed_iupred3.tsv" "${outdir}/netMHCpan/netMHCpan_SB_WB.tsv" "${outdir}/ELM_search/ELM_SLiMs_overview.tsv" \
        -o "${outdir}/run_table.tsv"

# "${outdir}/BLASTp/parsed_BLASTp.tsv"
