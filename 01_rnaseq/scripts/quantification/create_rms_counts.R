basedir = "/hpc/pmc_vanheesch/projects/Jip/rms_analyis/01_rnaseq"
tumor_type = 'RMS'
rds_loc = paste(basedir,'analysis','quantification','deseq2',sep="/")

rms_fp_cohort = read.delim(paste(basedir, "documentation", "RMS_fp.txt", sep = "/"), sep = ";")
rms_fn_cohort = read.delim(paste(basedir, "documentation", "RMS_fn.txt", sep = "/"), sep = ";")
rms_fp_org = read.delim(paste(basedir, "documentation", "RMS_fp_org.txt", sep = "/"), sep = ";")
rms_fn_org = read.delim(paste(basedir, "documentation", "RMS_fn_org.txt", sep = "/"), sep = ";")

rms_meta <- rbind(rms_fp_cohort,rms_fn_cohort,rms_fp_org,rms_fn_org)
rownames(rms_meta) <- rms_meta$sample_id

count_files <- list.files(paste(basedir,"analysis","quantification","salmon_quant",
                                sep  = "/"),
                          recursive = T,
                          pattern = "quant.sf",
                          full.names = T)
names(count_files) <- basename(gsub("/quant.sf","",count_files))

txdb = paste(basedir,"analysis","rnaseq_pipeline","customannotation",
             paste0(tumor_type,"_container"),
             paste0(tumor_type, "_full_novel_filtered_corrected.gtf_TxDb"), sep = "/")

txdb <- AnnotationDbi::loadDb(txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

count_files_rms <- count_files[which(names(count_files) %in% rms_meta$sample_id)]
rms_meta <- rms_meta[which(rms_meta$sample_id %in% names(count_files_rms)),]

txi_rms <- tximport::tximport(count_files_rms, type = "salmon", 
                              tx2gene = tx2gene, 
                              dropInfReps = T, 
                              countsFromAbundance = "scaledTPM")

all(rownames(rms_meta) == colnames(txi_rms$counts))
all(colnames(txi_rms$counts) == colnames(txi_rms$length))

rms_meta <- rms_meta[colnames(txi_rms$counts),]
all(rownames(rms_meta) == colnames(txi_rms$counts))

dds_rms <- DESeq2::DESeqDataSetFromTximport(txi = txi_rms,
                                            colData = rms_meta,
                                            design = ~ condition)

saveRDS(dds_rms,
        file = paste(rds_loc,"RMS_counts_org_patient.RDS", sep = "/"))

