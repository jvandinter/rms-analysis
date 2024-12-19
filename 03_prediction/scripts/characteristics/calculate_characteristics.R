library(Peptides)
library(dplyr)
library(Biostrings)

# define the arguments needed to run the script
args <- commandArgs(trailingOnly = TRUE)

# check whether the arguments are inputted
if (length(args) < 2) {
  stop("Usage: Rscript calculate_characteristics.R input_file output_file")
}

# assign to variable
input_file <- args[1]
output_file <- args[2]

# fasta to table
fasta_sequences <- readAAStringSet(input_file)

# Extract protein IDs and sequences
protein_ids <- names(fasta_sequences)
protein_sequences <- as.character(fasta_sequences)

# Create a data frame with protein IDs and sequences
protein_table <- data.frame(protein_id = protein_ids, protein_sequence = protein_sequences, stringsAsFactors = FALSE)

# calculate hydrophobicity, weight, charge etc
protein_characteristics <- protein_table %>% 
  mutate(length = nchar(protein_sequence),
         hydrophobicity = hydrophobicity(protein_sequence, "Miyazawa"),
         moleculairweight = mw(protein_sequence),
         mass_over_charge = mz(protein_sequence),
         instability =instaIndex(protein_sequence),
         pi_isoelectricpoint = pI(protein_sequence), #iso electric point
         charge = charge(protein_sequence)) 

# Write results to analysis directory on HPC
write.table(protein_characteristics, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)