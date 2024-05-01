suppressPackageStartupMessages({
  library(magrittr)
  library(DESeq2)
  library(dplyr)
})

# Set global variables ----------------------------------------------------

basedir = "/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis"
rnadir = paste(basedir,"01_rnaseq" , sep = "/")
meta_cohort = read.delim(paste(rnadir, "documentation", "tidy_rms_quant_meta.txt", sep = "/"), sep = " ") %>%
  dplyr::filter(cohort %in% c("AML","ATRT","B-ALL","EPN","EWS","MBL","NBL","OS","T-ALL","WT")) %>%
  dplyr::distinct()
meta_rms = read.delim(paste(rnadir, "documentation", "tidy_rms_quant_meta.txt", sep = "/"), sep = " ") %>%
  dplyr::filter(cohort == "RMS")

txi <-
  readRDS(paste(rnadir, "analysis","quantification_all_cohorts", "RMS_counts_all_cohorts.RDS", sep = "/"))

savedir = paste(rnadir, "results","quantification", "pediatric", sep = "/")
rds_loc = paste(rnadir, "analysis", "quantification", "deseq2","pediatric",sep = "/")
  
print(paste(Sys.time(),"Analysing pediatric cohorts"))
  
# Combine sample names per cohort
meta <- rbind(meta_rms, meta_cohort)
meta$cohort <- factor(meta$cohort, levels = c(unique(meta_rms$cohort),unique(meta_cohort$cohort)))
samples <- meta$sample_id
rownames(meta) <- samples
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
                                             design = ~ cohort)

gtex_dds <- DESeq2::DESeq(gtex_dds)

saveRDS(gtex_dds,
        file = paste(rds_loc,"pediatrics_cohort.RDS", sep = "/"))

resultsNames(gtex_dds)

for (i in 2:length(levels(meta$cohort))) {
  print(paste("writing results for", levels(meta$cohort)[i]))
  res = results(gtex_dds, contrast = c("cohort", levels(meta$cohort)[1], 
                                       levels(meta$cohort)[i]))
  write.table(res, file = paste(savedir, paste0(levels(meta$cohort)[1],"_vs_",
                                                levels(meta$cohort)[i], ".txt"),
                                sep = "/"))
}

print(paste(Sys.time(),"Finished DE for pediatric cohorts"))
