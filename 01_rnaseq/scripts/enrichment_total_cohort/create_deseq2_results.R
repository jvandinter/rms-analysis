suppressPackageStartupMessages({
  library(magrittr)
  library(DESeq2)
  library(dplyr)
})

# Set global variables ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
tumor_type <- args[1]
group <- args[2]

base_dir = "/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis"
rna_dir = paste(base_dir,"01_rnaseq" , sep = "/")
meta <- read.delim(paste(rna_dir,"documentation/all_cohorts","rms_meta_all_combined.csv", sep = "/"), 
                   sep = ",")

txi <-
  readRDS(paste(rna_dir, "analysis","quantification", "RMS_counts_all_cohorts.RDS", sep = "/"))

save_dir = paste(rna_dir, "results","quantification", tumor_type, sep = "/")
rds_loc = paste(rna_dir, "analysis", "quantification", tumor_type, sep = "/")
  
print(paste(Sys.time(),"Analysing pediatric cohorts for",tumor_type))
  
# Combine sample names per cohort
if(tumor_type == "RMS") {
  meta <- meta %>%
    dplyr::mutate(type = dplyr::case_when(cohort == tumor_type ~ "RMS",
                                          TRUE ~ type))
  meta_rms <- meta %>%
    dplyr::filter(type == tumor_type)
} else {
  meta_rms <- meta %>%
    dplyr::filter(type == tumor_type)
}

meta_compare <- meta %>% 
  dplyr::filter(!(cohort == "RMS"))
comparisons <- unique(meta_compare$type)

if (group %in% c(1,4,7)) {
  meta_subset <- meta_compare %>%
    dplyr::filter(type %in% comparisons[1:19])
} else if (group %in% c(2,5,8)) {
  meta_subset <- meta_compare %>%
    dplyr::filter(type %in% comparisons[20:38])
} else if (group %in% c(3,6,9)) {
  meta_subset <- meta_compare %>%
    dplyr::filter(type %in% comparisons[39:57])
}

meta <- rbind(meta_subset,meta_rms)

meta$type <- factor(meta$type, levels = c(tumor_type,unique(meta_subset$type)))
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
                                             design = ~ type)

gtex_dds <- DESeq2::DESeq(gtex_dds)

saveRDS(gtex_dds,
        file = paste0(rds_loc,"/rms_dds_",tumor_type,"_",group,".RDS"))

resultsNames(gtex_dds)

for (i in 2:length(levels(meta$type))) {
  print(paste("writing results for", levels(meta$type)[i]))
  res = results(gtex_dds, contrast = c("type", levels(meta$type)[1], 
                                       levels(meta$type)[i]))
  write.table(res, file = paste(save_dir, paste0(levels(meta$type)[1],"_vs_",
                                                levels(meta$type)[i], ".txt"),
                                sep = "/"))
}

print(paste(Sys.time(),"Finished DE for pediatric cohorts"))
