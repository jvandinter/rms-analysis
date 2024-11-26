library(tidyverse)
library(rtracklayer)
library(biomaRt)
library(Biostrings)
library(parallel)

no_cores <- 12

basedir = "/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis"
ribo_dir = paste(basedir,
                 "02_riboseq", sep = "/")

mat_cutoff_summary <- read.delim(paste(ribo_dir,
                                       "analysis/orf_quantification/ORF_sample_sharing_summary.csv",
                                       sep = "/"), 
                                 sep = ",",
                                 row.names = 1)

orfquant_results_location <- paste(ribo_dir,
                                   "analysis/ORFquant/RMS_merged_psites",
                                   "RMS_merged_psites_final_ORFquant_results",sep="/")

# Load ORF definitions
rms_orfquant_orfs <- get(load(orfquant_results_location))
orfs_tx <- rms_orfquant_orfs$ORFs_tx
orfs_tx_df <-  data.frame(mcols(orfs_tx)[, c("ORF_id_tr", "Protein", "gene_id", 
                                             "gene_biotype", "gene_name", 
                                             "transcript_id", "transcript_biotype", 
                                             "ORF_category_Tx", "ORF_category_Gen", 
                                             "P_sites_raw", "P_sites_raw_uniq")]) %>%
  distinct()

# Generate table genomic ORF locations
orfs_locations <- data.frame(rms_orfquant_orfs$ORFs_gen) %>%
  mutate(ORF_id_tr = names(rms_orfquant_orfs$ORFs_gen)) %>%
  group_by(ORF_id_tr) %>%
  mutate(ORF_ranges = paste0(seqnames, ":", min(start), "-", max(end))) %>%
  dplyr::select(c("ORF_id_tr", "ORF_ranges")) %>%
  distinct()

# Add uniprot IDs
require('biomaRt')
mart = useEnsembl("ENSEMBL_MART_ENSEMBL")  # Use same ensembl version as used for generating custom annotation
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'uniprot_gn_id'),
  uniqueRows=FALSE) %>%
  group_by(ensembl_gene_id) %>% 
  mutate(uniprot_gn_ids = paste0(uniprot_gn_id, collapse = ";")) %>%
  ungroup() %>%
  dplyr::select(c(ensembl_gene_id, uniprot_gn_ids)) %>%
  distinct()

# Generate ORF table
orfs_extended <- mat_cutoff_summary %>%
  tibble::rownames_to_column(var = "orf_id") %>%
  left_join(orfs_tx_df, by = c("orf_id" = "ORF_id_tr")) %>%
  left_join(orfs_locations, by = c("orf_id" = "ORF_id_tr")) %>%
  left_join(annotLookup, by = c("gene_id" = "ensembl_gene_id")) # Add uniprot column

# Load UniProt Reference Proteome fasta files
uniprot_fasta <- rtracklayer::import("/hpc/pmc_vanheesch/projects/Damon/Neuroblastoma_neoantigens/UP000005640_9606.fasta", type="AA")
uniprot_fasta_additional <- rtracklayer::import("/hpc/pmc_vanheesch/projects/Damon/Neuroblastoma_neoantigens/UP000005640_9606_additional.fasta", type="AA")

uniprot_fasta_full <- c(uniprot_fasta, uniprot_fasta_additional)
# uniprot_fasta_full <- c(uniprot_fasta)

uniprot_fasta_names <- uniprot_fasta_full 
names(uniprot_fasta_names) <- sapply(names(uniprot_fasta_full), function(x) {
  strsplit(x, "\\|")[[1]][2]
})

prot_ids <- orfs_extended$uniprot_gn_ids
prot_ids_split <- unique(unlist(strsplit(x = prot_ids, ";")))

process_sequence <- function(i) {
  prot_ids <- orfs_extended$uniprot_gn_ids[i]
  if (is.na(prot_ids)) return(NA)
  
  orf_sequence <- AAString(orfs_extended$Protein[i])
  prot_ids_split <- strsplit(prot_ids, ";")[[1]]
  
  prot_ids_split <- prot_ids_split[prot_ids_split %in% names(uniprot_fasta_names)]
  fasta_entries <- uniprot_fasta_names[prot_ids_split]
  fasta_entries <- fasta_entries[!sapply(fasta_entries, is.null)]
  
  if (length(fasta_entries) == 0) return(NA)
  
  # prot_seqs <- AAStringSet(unlist(fasta_entries))
  alignment <- pairwiseAlignment(pattern = fasta_entries, subject = orf_sequence, type = "local")
  similarity_score <- max(nmatch(alignment) / width(fasta_entries)) * 100
  
  return(similarity_score)
}

orfs_extended$similarity_score <- parallel::mclapply(1:nrow(orfs_extended), process_sequence, mc.cores = no_cores)
orfs_extended$similarity_score <- as.numeric(orfs_extended$similarity_score)

write.csv(orfs_extended, paste(ribo_dir,
                               "results/orf_reannotation",
                               "RMS_merged_ORFquant_orfs_simScore.csv",
                               sep = "/"),
          quote = F,
          sep = ",")
