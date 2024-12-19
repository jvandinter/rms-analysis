#!/bin/bash
#SBATCH --job-name create_count_matrix
#SBATCH --mem=120G
#SBATCH --time=12:00:00

module load R/4.2.1

/usr/bin/env Rscript -<<EOF
library(ggplot2)
library(DESeq2)
library(magrittr)
library(readr)

basedir="/hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/"
workdir=paste(basedir,"20221020_JD_quant_tumor_cohorts", sep = "/")
savedir=paste(workdir, "analysis", sep = "/")
nbl_gtf=paste(basedir,"/20221020_JD_quant_tumor_cohorts/data/processed/customannotation/NBL_complete_novel_filtered_corrected.gtf", sep = "/")
nbl_txdb=paste(basedir,"/20221020_JD_quant_tumor_cohorts/data/processed/customannotation/NBL_complete/NBL_complete_novel_filtered_corrected.gtf_TxDb", sep = "/")
gtf <- rtracklayer::import.gff(nbl_gtf)
txdb <- AnnotationDbi::loadDb(nbl_txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

count_files <- list.files(paste(workdir,"data/processed/salmon_quant/NBL",sep  = "/"), recursive = T, pattern = "quant.sf", full.names = T)
names(count_files) <- basename(gsub("/quant.sf","",count_files))

txi <- tximport::tximport(count_files, type = "salmon", tx2gene = tx2gene, dropInfReps = T, countsFromAbundance = "scaledTPM")

saveRDS(txi, file = paste(savedir, "NBL_counts_full.RDS", sep = "/"))
print("Count Matrix successfully created!")
EOF
