library(magrittr)
library(DESeq2)

tumor_type="RMS"
pmc_tumor_code=c("Alveolar rhabdomyosarcoma","Embryonal rhabdomyosarcoma, NOS","Spindle cell rhabdomyosarcoma","Rhabdomyosarcoma, NOS")
sj_tumor_code=c("ARMS","ERMS","SCRMS","RMS")
basedir="/hpc/pmc_vanheesch/projects/Jip/custom_transcriptomes"
workdir="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq"
savedir="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/results/"

meta_gtex=read.delim(paste(workdir,"metadata","metadata_gtex_FP.txt", sep = "/"), sep = ";")
txi <- readRDS(paste(savedir,"processed", "RMS_counts_full.RDS", sep = "/"))

# Combine sample names per cohort
samples_gtex <- meta_gtex$sample_id

# RMS + GTEX
txi_gtex <- list(abundance = txi$abundance[,which(colnames(txi$abundance) %in% c(samples_gtex))],
                counts = txi$counts[,which(colnames(txi$counts) %in% c(samples_gtex))],
                length = txi$counts[,which(colnames(txi$length) %in% c(samples_gtex))],
                countsFromAbundance = "scaledTPM")

meta_gtex$type <- ifelse(meta_gtex$type %in% c("Cervix Uteri","Uterus","Vagina","Ovary","Fallopian Tube"),"female_reproductive",
                         ifelse(meta_gtex$type %in% c("Esophagus","Colon","Small Intestine","Stomach"),"digestive_system",
                                ifelse(meta_gtex$batch %in% c("SJ","PMC"),tumor_type,
                                       ifelse(meta_gtex$type %in% c("Prostate","Testis"),"male_reproductive",meta_gtex$type))))
meta_gtex$type <- factor(meta_gtex$type,levels = unique(meta_gtex$type))
meta_gtex$sex <- factor(meta_gtex$sex, levels = c("male", "female","not available"))
meta_gtex$batch <- factor(meta_gtex$batch, levels = c("PMC","SJ","GTEX"))
meta_gtex$condition <- factor(meta_gtex$condition, levels = c("healthy",tumor_type))

# Final check
rownames(meta_gtex) <- meta_gtex$sample_id

all(rownames(meta_gtex) == colnames(txi_gtex$counts))
all(colnames(txi_gtex$counts) == colnames(txi_gtex$length))

meta_gtex <- meta_gtex[colnames(txi_gtex$counts),]
all(rownames(meta_gtex) == colnames(txi_gtex$counts))

gtex_dds <- DESeq2::DESeqDataSetFromTximport(txi = txi_gtex,
                                            colData = meta_gtex,
                                            design = ~ type)

# Keep samples with more than 10 counts in more than 10 samples
keep <- rowSums(counts(gtex_dds) >= 10) >= 10
gtex_dds <- gtex_dds[keep,]

gtex_dds <- DESeq2::DESeq(gtex_dds)

saveRDS(gtex_FP_dds, file = paste(savedir,"results","GTEX_DESEQ_type.RDS", sep = "/"))

resultsNames(gtex_dds)

for (i in 2:length(levels(meta_gtex$type))) {
  print(paste("writing results for",levels(meta_gtex$type)[i]))
  res = results(gtex_dds, contrast = c("type",levels(meta_gtex$type)[1],levels(meta_gtex$type)[i]))
  write.table(res, file = paste(savedir, "results",paste0("GTEX_res_FP_",gsub(" ","_",levels(meta_gtex$type)[i]),".txt"), sep= "/"))
}

EOF
