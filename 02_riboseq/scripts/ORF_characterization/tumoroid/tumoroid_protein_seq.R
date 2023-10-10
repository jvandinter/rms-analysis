library(magrittr)
library(Biostrings)
library(dplyr)

fasta_locs = grep("ORG",list.dirs("/hpc/pmc_vanheesch/projects/Jip/rms_analysis/02_riboseq/analysis/ORFquant", recursive = F), value = T)
savedir = "/hpc/pmc_vanheesch/projects/Jip/rms_analysis/02_riboseq/results"
interesting_orfs = read.delim(paste(savedir,"tumoroid","RMS_tumoroid_all_ORFs.txt",sep="/"), sep =";")
savefile = paste(savedir,"tumoroid","tumoroid_sequence_table.txt", sep = "/")

# code to extract all proteins from the file(s)
processed_ids <- character(0)
matching_seq <- character(0)
sequence_table <- data.frame(orf_id = character(0), sequence = character(0), width = integer(0), names = character(0), stringsAsFactors = FALSE)
filtered_interesting_orfs <- interesting_orfs %>%
  replace(is.na(x = .), "c") %>%
  dplyr::filter(!(class_code == "k"))

for (i in 1:nrow(filtered_interesting_orfs)) {
  id <- filtered_interesting_orfs$orf_id[i]
  source_file_number <- filtered_interesting_orfs$source[i]
  if (!(id %in% processed_ids)) {
    fasta_file <- list.files(grep(source_file_number,fasta_locs, value = T), pattern = "fasta", full.names = T)
    fasta <- Biostrings::readAAStringSet(fasta_file)
    matching_fasta <- fasta[grepl(id, names(fasta))]
    # Extract sequences, widths, and names
    matching_sequences <- as.character(matching_fasta)
    matching_widths <- width(matching_fasta)
    matching_names <- names(matching_fasta)
    # Add the extracted data to the sequence_table data frame
    if (length(matching_sequences) > 0) {
      new_rows <- data.frame(orf_id = id, sequence = matching_sequences, width = matching_widths, names = matching_names, stringsAsFactors = FALSE)
      sequence_table <- rbind(sequence_table, new_rows)
    }
    processed_ids <- c(processed_ids, id)
  }
}

write.table(x = sequence_table, file = savefile,quote=F,sep="\t", row.names = F)
