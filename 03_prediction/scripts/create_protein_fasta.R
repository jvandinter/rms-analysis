library(magrittr)
library(dplyr)
library(tidyr)

basedir = "/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis"
rna_dir = paste(basedir,"01_rnaseq",sep = "/")
ribo_dir = paste(basedir,"02_riboseq", sep = "/")
save_dir = paste(basedir,"poster_plots/figures",sep="/")
metadata_dir = paste(rna_dir,"documentation",sep="/")
txome_gtf = paste(rna_dir,
                  "analysis/rnaseq_pipeline/customannotation",
                  "RMS_full_novel_filtered_corrected.sorted.gtf", sep = "/")
tumor_type = "RMS"
rds_loc = paste(rna_dir,"analysis","quantification", sep = "/")


fp_table <- read.delim(paste(rna_dir,
                             "results/quantification/",
                             "RMS-FP_enrichment.csv", sep="/"), 
                       sep = ",")
fp_specific <- fp_table %>% 
  dplyr::filter(specific == T &
                  gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
  dplyr::pull(gene_id)

fn_table <- read.delim(paste(rna_dir, 
                             "results/quantification/",
                             "RMS-FN_enrichment.csv", sep="/"), 
                       sep = ",")
fn_specific <- fn_table %>% 
  dplyr::filter(specific == T &
                  gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
  dplyr::pull(gene_id)

rms_table <- read.delim(paste(rna_dir, "results/quantification/",
                              "RMS_enrichment.csv", sep="/"), 
                        sep = ",")
rms_specific <- rms_table %>% 
  dplyr::filter(specific == T &
                  gene_biotype %in% c("lncRNA","protein_coding","stringtie")) %>%
  dplyr::pull(gene_id)

specific_genes <- c(fp_specific, 
                    rms_specific,
                    fn_specific) %>% unique()


# Translated ORFs
intorfs <- read.delim(paste(ribo_dir, 
                            "results/orf_reannotation",
                            "jorge_psites_intorfs.csv",
                            sep = "/"),sep = ",")

orf_table <- read.delim(file = paste(ribo_dir,"results/orf_reannotation",
                                     "RMS_harmonised_ORF_table.csv", sep = "/"),
                        sep = ",")

orf_categories_df <- orf_table %>%
  dplyr::filter(translated == T &
                gene_id %in% specific_genes) %>%
  dplyr::mutate(orf_category_new = dplyr::case_when(orf_category_new == "nested_ORF" &
                                                      !(orf_id %in% intorfs$orf_id) ~ "ORF_annotated",
                                                    grepl("truncation|extension",orf_category_new) ~ "ORF_annotated",
                                                    TRUE ~ orf_category_new),
                Protein = ifelse(grepl("\\*$",Protein),
                                 Protein,
                                 paste0(Protein,"*")))

write.table(orf_categories_df, file = paste(basedir, "03_prediction/documentation",
                                            "orf_table.csv",
                                            sep = "/"),
            sep = ",", quote = F, row.names = F)

# Create peptides fasta

df_to_fasta <- function(df, file) {
  fasta <- apply(df, 1, function(row) {
    paste0(">", row["orf_id"], "\n", row["Protein"])
  })
  writeLines(fasta, file)
}

# Convert dataframe to FASTA file
df_to_fasta(orf_categories_df, paste(basedir,"03_prediction/documentation","RMS_enriched.fasta", sep = "/"))


