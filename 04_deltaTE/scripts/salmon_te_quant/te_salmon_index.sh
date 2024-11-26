#!/bin/bash

set -euo pipefail

mkdir -p "${outdir}/salmon_index"

# Set locations

merged_cds_gtf_output="${outdir}/RMS_CDS.gtf"
#merged_cds_gtf_output="/hpc/pmc_vanheesch/shared_resources/GENOMES/Homo_sapiens.GRCh38/102/annotation/Homo_sapiens.GRCh38.102.gtf"
fixed_cds_gtf_out="${outdir}/salmon_index/tx_cds.gtf"
fixed_cds_fasta_out="${outdir}/salmon_index/tx_cds.fa"
gentrome_out="${outdir}/salmon_index/gentrome.fa"
decoy_out="${outdir}/salmon_index/decoys.txt"

echo "`date` fixing GTF"

# Fix GTF file to give all CDS regions an exon and parent transcript
apptainer exec -B /hpc:/hpc ${container_dir}/gffread-0.12.6.sif gffread \
  -g ${reference_genome} \
  -o "${fixed_cds_gtf_out}" \
  -F \
  --keep-exon-attrs \
  --force-exons \
  -T \
  ${merged_cds_gtf_output}

echo "`date` Converting to FASTA"

# Extract CDS sequences into fasta file
apptainer exec -B /hpc:/hpc ${container_dir}/gffread-0.12.6.sif gffread \
  -g "${reference_genome}" \
  -x "${fixed_cds_fasta_out}" \
  "${fixed_cds_gtf_out}"

echo "`date` creating gentrome"

# Make gentrome
cat "${fixed_cds_fasta_out}" "${reference_genome}" > ${gentrome_out}

echo "`date` creating salmon decoy file"

# Make salmon decoy
grep "^>" "${reference_genome}" | cut -d " " -f 1 > ${decoy_out}
sed -i.bak -e 's/>//g' ${decoy_out}

echo "`date` creating salmon index"

# Create salmon index
apptainer exec -B /hpc:/hpc ${container_dir}/salmon-1.8.0.sif salmon index \
  --transcripts ${gentrome_out} \
  --decoys ${decoy_out} \
  --index "${outdir}/salmon_index/" \
  --kmerLen 13 \
  --keepDuplicates \
  --threads 4

  echo "`date` Done!"
