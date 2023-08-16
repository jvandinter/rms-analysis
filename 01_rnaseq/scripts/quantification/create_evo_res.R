suppressPackageStartupMessages({
    library(magrittr)
    library(DESeq2)
})

tumor_type="RMS"
pmc_tumor_code=c("Alveolar rhabdomyosarcoma","Embryonal rhabdomyosarcoma, NOS","Spindle cell rhabdomyosarcoma","Rhabdomyosarcoma, NOS")
sj_tumor_code=c("ARMS","ERMS","SCRMS","RMS")
basedir="/hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes"
workdir=paste(basedir,"20221020_JD_quant_tumor_cohorts", sep = "/")
savedir=paste(workdir, "analysis",tumor_type, sep = "/")

meta_evo=read.delim(paste(savedir,"metadata","metadata_evo.txt", sep = "/"), sep = ";")
txi <- readRDS(paste(savedir,"processed", "RMS_counts_full.RDS", sep = "/"))

# Combine sample names per cohort
samples_evo <- meta_evo$sample_id

meta_evo$tissue <- factor(meta_evo$tissue, levels = c(tumor_type,"forebrain","hindbrain","heart","kidney","liver","ovary","testis"))
meta_evo$batch <- factor(meta_evo$batch, levels = c("PMC","SJ","EVO"))
meta_evo$stage <- factor(meta_evo$stage, levels = c("adult","child","fetal",tumor_type))
meta_evo$fetal <- factor(meta_evo$fetal, levels = c("postbirth","fetal",tumor_type))
meta_evo$fetal_tissue <- factor(meta_evo$fetal_tissue, levels = unique(meta_evo$fetal_tissue))

# RMS + evo
txi_evo <- list(abundance = txi$abundance[,which(colnames(txi$abundance) %in% c(samples_evo))],
                counts = txi$counts[,which(colnames(txi$counts) %in% c(samples_evo))],
                length = txi$length[,which(colnames(txi$length) %in% c(samples_evo))],
                countsFromAbundance = "scaledTPM")

# Final check
rownames(meta_evo) <- meta_evo$sample_id

all(rownames(meta_evo) == colnames(txi_evo$counts))
all(colnames(txi_evo$counts) == colnames(txi_evo$length))

meta_evo <- meta_evo[colnames(txi_evo$counts),]
all(rownames(meta_evo) == colnames(txi_evo$counts))

evo_dds <- DESeq2::DESeqDataSetFromTximport(txi = txi_evo,
                                            colData = meta_evo,
                                            design = ~ fetal_tissue)

# Keep samples with more than 10 counts in more than 3 samples
keep <- rowSums(counts(evo_dds) >= 10) >= 10
evo_dds <- evo_dds[keep,]

evo_dds <- DESeq2::DESeq(evo_dds)

saveRDS(evo_dds, file = paste(savedir,"results","EVO_DESEQ_tissue.RDS", sep = "/"))

for (i in 2:(length(unique(meta_evo$fetal_tissue)))) {
  print(paste("writing results for",unique(meta_evo$fetal_tissue)[i]))
  res = results(evo_dds, contrast = c("fetal_tissue",unique(meta_evo$fetal_tissue)[1],unique(meta_evo$fetal_tissue)[i]))
  write.table(res, file = paste(savedir, "results",paste0("EVO_res_",gsub(" ","_",unique(meta_evo$fetal_tissue)[i]),".txt"), sep= "/"))
}
