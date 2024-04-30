suppressPackageStartupMessages({
  library(magrittr)
  library(DESeq2)
})

basedir = "/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq"

rms_df <- read.delim(paste(basedir, "documentation", "RMS_df.txt", sep = "/"), sep = ";")
meta_cohort = read.delim(paste(basedir, "documentation", "RMS_R2_atlas.txt", sep = "/"), sep = ";")

cohort = "R2_atlas"

txi <-
  readRDS(paste(basedir, "analysis","quantification","deseq2", "RMS_counts_full.RDS", sep = "/"))

for (j in 1:nrow(rms_df)) {

  analysis_name = rms_df$analysis_type[j]
  savedir = paste(basedir, "results","quantification", analysis_name, sep = "/")
  rds_loc = paste(basedir, "analysis", "quantification", "deseq2",sep = "/")
  meta_rms = read.delim(paste(basedir, "documentation", rms_df$metadata_name[j], sep = "/"), sep = ";")
  
  print(paste(Sys.time(),"Analysing;",analysis_name,"for cohort;",cohort))
  
  # Combine sample names per cohort
  meta <- rbind(meta_rms, meta_cohort)
  meta$condition <- factor(meta$condition, levels = c(unique(meta_rms$condition),unique(meta_cohort$condition)))
  samples <- meta$sample_id
  print(nrow(meta))
  
  txi_subset <-
    list(
      abundance = txi$abundance[, which(colnames(txi$abundance) %in% c(samples))],
      counts = txi$counts[, which(colnames(txi$counts) %in% c(samples))],
      length = txi$length[, which(colnames(txi$length) %in% c(samples))],
      countsFromAbundance = "scaledTPM"
    )
  
  all(rownames(meta) == colnames(txi_subset$counts))
  all(colnames(txi_subset$counts) == colnames(txi_subset$length))
  
  meta <- meta[colnames(txi_subset$counts),]
  all(rownames(meta) == colnames(txi_subset$counts))
  
  gtex_dds <- DESeq2::DESeqDataSetFromTximport(txi = txi_subset,
                                               colData = meta,
                                               design = ~ condition)
  
  gtex_dds <- DESeq2::DESeq(gtex_dds)
  
  saveRDS(gtex_dds,
          file = paste(rds_loc,paste0(cohort,"_",analysis_name,"_condition.RDS"), sep = "/"))
  
  resultsNames(gtex_dds)
  
  for (i in 2:length(levels(meta$condition))) {
    print(paste("writing results for", levels(meta$condition)[i]))
    res = results(gtex_dds, contrast = c("condition", levels(meta$condition)[1], 
                                         levels(meta$condition)[i]))
    write.table(res, file = paste(savedir, paste0(cohort,"_",analysis_name,"_", 
                                                  gsub(" ", "_", levels(meta$condition)[i]), ".txt"),
    sep = "/"))
  }
}

for (j in 1:nrow(rms_df)) {
  
  print(rms_df$analysis_type[j])
  print(paste(basedir, "results","quantification", analysis_name, sep = "/"))
  print(paste(basedir, "analysis", "quantification", "deseq2",sep = "/"))
  print(paste(basedir, "documentation", rms_df$metadata_name[j], sep = "/"), sep = ";")
}
