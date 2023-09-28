library(tidyverse)
library(magrittr)
library(Biostrings)
library(stringr)
library(dplyr)

# load arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript SB_WB_netMHCpan_overview.R input_fasta_file input_directory elm_classes_input output_file")
}

input_fasta_file <- args[1]
input_directory <- args[2]
elm_classes_input <- args[3]
output_file <- args[4]

# load in fasta file and regex overview

elm_classes <- read.table(elm_classes_input, sep = '\t', header = TRUE)

fasta_sequences <- readAAStringSet(input_fasta_file)

# regex search in whole protien
# Create an empty data frame to store the results
regex_whole_protein <- data.frame(Protein_Name = character(),
                                  Regex_Name_WP = character(),
                                  Matched_Sequence_WP = character(),
                                  stringsAsFactors = FALSE)

# Iterate over each sequence in the FASTA file
for (i in 1:length(fasta_sequences)) {
  # Get the current sequence and its name
  current_sequence_WP <- fasta_sequences[[i]]
  current_sequence_name_WP <- names(fasta_sequences)[i]
  
  # Iterate over regex patterns
  for (j in 1:nrow(elm_classes)) {
    # Get the regex pattern and its name from the table
    regex_pattern <- elm_classes$Regex[j]
    regex_name <- elm_classes$ELMIdentifier[j]
    
    # Convert the sequence to a character string
    current_sequence_string <- as.character(current_sequence_WP)
    
    # Use gregexpr to find all matches
    match_positions <- gregexpr(regex_pattern, current_sequence_string, ignore.case = TRUE)
    
    # Store the matching sequences in the results data frame if there are matches
    matching_sequences <- unlist(lapply(match_positions, function(x) substring(current_sequence_string, x)))
    
    if (length(matching_sequences) > 0 && !is.element(current_sequence_string, matching_sequences)) {
      regex_whole_protein <- rbind(regex_whole_protein, data.frame(Protein_Name = current_sequence_name_WP,
                                                                   Regex_Name_WP = regex_name,
                                                                   Matched_Sequence_WP = matching_sequences,
                                                                   stringsAsFactors = FALSE))
    }
  }
}

# Group the results by Protein_Name
summary_table_whole_protein <- regex_whole_protein %>%
  group_by(Protein_Name) %>%
  summarize(Regex_Names_WP = paste(unique(Regex_Name_WP), collapse = ","),
            Sequences_WP = paste(unique(Matched_Sequence_WP), collapse = ","),
            Regex_Count_WP = n_distinct(Regex_Name_WP))

# regex search in disorderd regions
# first process the disorder files
# Directory containing the disorder data files
disorder_data_directory <- input_directory
disorder_data_files <- list.files(disorder_data_directory, pattern = "\\.txt$", full.names = TRUE)

# Function to read a score file and process the data
process_disorder_data_file <- function(file_path) {
  data <- read.table(file_path)
  
  # Assuming the residue column is the 1st column (index 1) and the score column is the 2nd column (index 2)
  colnames(data) <- c("group", "residue", "score")
  
  # Convert the score column to numeric if it's not already numeric
  if (!is.numeric(data$score)) {
    data$score <- as.numeric(data$score)
  }
  
  # Create a new column "grp" to identify consecutive residues with score > 0.5
  data$grp <- cumsum(c(TRUE, diff(data$score > 0.5) != 0))
  
  # Filter the data and perform grouping
  filtered_data <- data %>%
    filter(score > 0.5) %>%
    group_by(grp) %>%
    summarize(peptide = paste(residue, collapse = ""), score = max(score))
  
  # Add the "file" column and remove .txt extension from file name
  filtered_data$file <- gsub("\\.txt$", "", basename(file_path))
  
  return(filtered_data)
}

# Process each file separately and store results in a list
disorder_results_list <- lapply(disorder_data_files, process_disorder_data_file)

# Combine results from all files
disordered_parts <- do.call(rbind, disorder_results_list)

# Create an empty data frame to store the results
regex_disordered_parts <- data.frame(Protein_Name = character(),
                                     Regex_Name_DP = character(),
                                     Sequences_DP = character(),
                                     stringsAsFactors = FALSE)

# Iterate over each peptide sequence in the disordered_parts data frame
for (i in 1:nrow(disordered_parts)) {
  # Get the current sequence and its name
  current_sequence <- disordered_parts$peptide[i]
  current_sequence_name <- disordered_parts$file[i]  # Assuming the file column contains the name
  
  # Iterate over regex sequences
  for (j in 1:nrow(elm_classes)) {
    # Get the regex pattern and its name from the table
    regex_pattern <- elm_classes$Regex[j]
    regex_name <- elm_classes$ELMIdentifier[j]
    
    # Use gregexpr to find all matches
    match_positions <- gregexpr(regex_pattern, current_sequence, ignore.case = TRUE)
    
    # Store the matching sequences in the regex_disordered_parts data frame if there are matches
    matching_sequences <- unlist(lapply(match_positions, function(x) substring(current_sequence, x)))
    
    if (length(matching_sequences) > 0 && !is.element(current_sequence, matching_sequences)) {
      regex_disordered_parts <- rbind(regex_disordered_parts, data.frame(Protein_Name = current_sequence_name,
                                                                         Regex_Name_DP = regex_name,
                                                                         Sequences_DP = matching_sequences,
                                                                         stringsAsFactors = FALSE))
    }
  }
}

# Group the results by Protein_Name
summary_table_disordered_parts <- regex_disordered_parts %>%
  group_by(Protein_Name) %>%
  summarize(Regex_Names_DP = paste(unique(Regex_Name_DP), collapse = ","),
            Sequences_DP = paste(unique(Sequences_DP), collapse = ","),
            Regex_Count_DP = n_distinct(Regex_Name_DP))

# Get the unique Protein_Name values from the FASTA file
unique_proteins <- unique(names(fasta_sequences))

# Create a data frame with all unique proteins
all_proteins <- data.frame(Protein_Name = unique_proteins, stringsAsFactors = FALSE)

# Left join all_proteins with final_result to include matching proteins
summary_table_disordered_parts_complete <- left_join(all_proteins, summary_table_disordered_parts, by = "Protein_Name")

# combine all results
overview_SLiMs <- merge(summary_table_whole_protein, summary_table_disordered_parts_complete, by = "Protein_Name", all = TRUE)

# Rename ID column of protein_id
names(overview_SLiMs)[names(overview_SLiMs) == "Protein_Name"] <- "protein_id"

# write results to hpc
write.table(overview_SLiMs, file = output_file, sep = '\t', quote = FALSE, row.names = FALSE)
