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

tumor_type="RMS"
basedir="/hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes/"
workdir=paste(basedir,"20221020_JD_quant_tumor_cohorts", sep = "/")
savedir=paste(workdir, "analysis", sep = "/")
tumor_txdb=paste0(basedir,"/20221020_JD_quant_tumor_cohorts/data/processed/customannotation/",
                 tumor_type,"_full/",
                 tumor_type,"_full_novel_filtered_corrected.gtf_TxDb")
txdb <- AnnotationDbi::loadDb(tumor_txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

count_files <- list.files(paste0(workdir,"/data/processed/salmon_quant/",tumor_type), recursive = T, pattern = "quant.sf", full.names = T)
names(count_files) <- basename(gsub("/quant.sf","",count_files))

txi <- tximport::tximport(count_files, type = "salmon", tx2gene = tx2gene, dropInfReps = T, countsFromAbundance = "scaledTPM")

saveRDS(txi, file = paste(savedir, paste0(tumor_type,"_counts_full.RDS"), sep = "/"))
print(paste0("Count Matrix successfully created for ",tumor_type,"!"))
EOF
