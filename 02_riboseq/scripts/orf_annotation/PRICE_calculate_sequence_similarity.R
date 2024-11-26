library(tidyverse)
library(rtracklayer)
library(biomaRt)
library(Biostrings)
library(parallel)

no_cores <- 12

basedir = "/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis"
ribo_dir = paste(basedir,
                 "02_riboseq", sep = "/")

orfs_extended <- read.delim(paste(ribo_dir,"results/orf_reannotation","price_orf_table.csv", sep = "/"), sep = ",")

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
                               "RMS_merged_PRICE_orfs_simScore.csv",
                               sep = "/"),
          quote = F,
          sep = ",")
