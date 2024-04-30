library(ggplot2)
library(DESeq2)
library(magrittr)
library(readr)

tumor_type="RMS"
basedir="/hpc/pmc_vanheesch/projects/jvandinter"
workdir=paste(basedir,"rms_analysis","01_rnaseq","analysis", sep = "/")
savedir=paste(workdir,"quantification_all_cohorts", sep = "/")
tumor_txdb=paste0(workdir,"/rnaseq_pipeline/customannotation/",
                 tumor_type,"_container/",
                 tumor_type,"_full_novel_filtered_corrected.gtf_TxDb")
txdb <- AnnotationDbi::loadDb(tumor_txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
meta_cohort <- read.delim(paste(basedir,"rms_analysis","01_rnaseq","documentation", "tidy_rms_quant_meta.txt", sep="/"), sep = " ") %>%
  dplyr::distinct()

print(paste(Sys.time(),"Finding salmon quant files ..."))
count_files <- list.files(paste0(basedir,"/quant_to_end_all_quants/analysis/quantification/",tumor_type), recursive = T, pattern = "quant.sf", full.names = T)
names(count_files) <- basename(gsub("/quant.sf","",count_files))
count_files <- count_files[which(names(count_files) %in% meta_cohort$sample_id)]

print(paste(Sys.time(),"Loading salmon quant files ..."))

txi <- tximport::tximport(count_files, type = "salmon", tx2gene = tx2gene, dropInfReps = T, countsFromAbundance = "scaledTPM")

print(paste(Sys.time(),"Loading salmon quant files. Done!"))

saveRDS(txi, file = paste(savedir, paste0(tumor_type,"_counts_all_cohorts.RDS"), sep = "/"))
print(paste(Sys.time(),paste0("Count Matrix successfully created for ",tumor_type,"!")))

  # Combine sample names per cohort
  meta_cohort$cohort <- factor(meta_cohort$cohort, levels = unique(meta_cohort$cohort))
  rownames(meta_cohort) <- meta_cohort$sample_id
  samples <- meta_cohort$sample_id
  print(nrow(meta_cohort))
  
  txi_subset <-
    list(
      abundance = txi$abundance[, which(colnames(txi$abundance) %in% c(samples))],
      counts = txi$counts[, which(colnames(txi$counts) %in% c(samples))],
      length = txi$length[, which(colnames(txi$length) %in% c(samples))],
      countsFromAbundance = "scaledTPM"
    )
  
  all(rownames(meta_cohort) == colnames(txi_subset$counts))
  all(colnames(txi_subset$counts) == colnames(txi_subset$length))
  
  meta_cohort <- meta_cohort[colnames(txi_subset$counts),]
  all(rownames(meta_cohort) == colnames(txi_subset$counts))

  print(paste(Sys.time(),"Creating DEseq2 count matrix ..."))
  
  dds <- DESeq2::DESeqDataSetFromTximport(txi = txi_subset,
                                               colData = meta_cohort,
                                               design = ~ 1)

  saveRDS(dds,
          file = paste(savedir,"RMS_dds_incl_OS_tumor_cohorts.RDS", sep = "/"))

  print(paste(Sys.time(),paste0("Count Matrix successfully created for ",tumor_type,"!")))
