suppressPackageStartupMessages({
  library(magrittr)
  library(DESeq2)
})

basedir = "/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq"

rms_fp_cohort = read.delim(paste(basedir, "documentation", "RMS_fp.txt", sep = "/"), sep = ";")
rms_fn_cohort = read.delim(paste(basedir, "documentation", "RMS_fn.txt", sep = "/"), sep = ";")
r2_cohort = read.delim(paste(basedir, "documentation", "RMS_R2_atlas.txt", sep = "/"), sep = ";")
gtex_cohort = read.delim(paste(basedir, "documentation", "RMS_GTEx.txt", sep = "/"), sep = ";")
evo_cohort = read.delim(paste(basedir, "documentation", "RMS_EVO.txt", sep = "/"), sep = ";")
meta_cohort = rbind(rms_fp_cohort,rms_fn_cohort,r2_cohort,gtex_cohort,evo_cohort)

txi <-
  readRDS(paste(basedir, "analysis","quantification","deseq2", "RMS_counts_full.RDS", sep = "/"))

  rds_loc = paste(basedir, "analysis", "quantification", "deseq2",sep = "/")

  # Combine sample names per cohort
  meta_cohort$condition <- factor(meta_cohort$condition, levels = unique(meta_cohort$condition))
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
  
  dds <- DESeq2::DESeqDataSetFromTximport(txi = txi_subset,
                                               colData = meta_cohort,
                                               design = ~ 1)

  saveRDS(dds,
          file = paste(rds_loc,"RMS_counts_combined.RDS", sep = "/"))
