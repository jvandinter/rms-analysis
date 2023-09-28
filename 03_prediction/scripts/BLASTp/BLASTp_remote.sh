#!/bin/bash

# node needed for executing apptainer, but due to using the -remote flag it will be submitted on the 
# blast servers which means there are minimal resouces needed to run this job
# apptainer executes command to perform remote blastp search on fasta file containing protein sequences
# it outputs a tsv file containing only the proteins WITH hits

mkdir "${outdir}/BLASTp_${run}" || { echo "Error: Failed to create the output directory."; exit 1; }

apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch "${apptainer_dir}/ncbi_blast-2.14.0.sif" \
    blastp \
        -db nr \
        -query "${input_fasta}" \
        -out "${outdir}/BLASTp_${run}/BLASTp_${run}.out" \
        -remote \
        -evalue 0.05 \
        -max_target_seqs 100 \
        -seg yes \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
        || { echo "Error: Failed to execute blastp search."; exit 1; }

# put header line on top of the output file

sed '1i qseqid  sseqid  pident  length  mismatch    gapopen qstart  qend    sstart  send    evalue  bitscore    qseq    sseq' "${outdir}/BLASTp_${run}/BLASTp_${run}.out" \
    || { echo "Error: Failed to add header line to the output file."; exit 1; }

# script to iterate over results to find proteins that do have and that do not have results

apptainer exec -B /hpc/pmc_vanheesch:/hpc/pmc_vanheesch "${apptainer_dir}/blast_parser_python-2.0.sif" \
    python3 /app/blastout_to_list.py \
        -b "${outdir}/BLASTp_${run}/BLASTp_${run}.out" \
        -f "${input_fasta}" \
        -o "${outdir}/BLASTp_${run}/parsed_BLASTp_${run}.tsv" \
        || { echo "Error: Failed to parse BLASTp results."; exit 1; }