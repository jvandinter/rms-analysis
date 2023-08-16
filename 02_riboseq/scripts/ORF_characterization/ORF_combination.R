message(paste(Sys.time(),"Loading required libraries ..."))
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
  library(RColorBrewer)
  library(parallel)
})

wd <- "/hpc/pmc_vanheesch/projects/Jip/rms_analysis/02_riboseq"
savedir <- paste(wd, "results", sep = "/")

gtf_file <- "/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/analysis/rnaseq_pipeline/customannotation/RMS_full_novel_filtered_corrected.gtf"
gtf_rannot_file <- paste("/hpc/pmc_vanheesch/projects/Jip/rms_analysis/01_rnaseq/analysis/rnaseq_pipeline/customannotation/","RMS_full","RMS_full_novel_filtered_corrected.gtf_Rannot", sep = "/")

thresholds <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

orf_methods <- c("ORFquant", "PRICE", "Ribo-TISH", "Ribotricer")

orfquant_files = list.files(paste(wd,"analysis","ORFquant", sep = "/"),
                            pattern = "*final_ORFquant_results",
                            recursive = T,
                            full.names = T)

prepare_orfs <- function(file_list, orf_method, ncores) {
  # Return a grangeslist with metadata columns representing the ORFs.
  #
  # Parameters:
  #   file_list - list, the input files from an ORF caller.
  #   orf_method - string, name of the ORF caller.
  #   ncores - int, number of cores to use.
  #
  # Returns:
  #   orfs - grangeslist, detected ORFs including the stop codon.
  if (orf_method == "PRICE") {
    orfs <- prepare_price(file_list, ncores)
    orfs <- lapply(orfs, function(x) {
      orf <- correct_price_geneid(x)
      m_cols <- mcols(orf)
      orf <- orf %>% as.data.frame() %>% mutate(group_name = name) %>% 
        makeGRangesListFromDataFrame(split.field = "group_name") 
      mcols(orf) <- m_cols
      return(orf)
    })
  } else if (orf_method == "ORFquant") {
    orfs <- prepare_orfquant(file_list)
    orfs <- lapply(orfs, function(orf) {
      orf <- add_stop_codon(orf, gtf_exon)
      return(orf)
    })
  } else if (orf_method == "Ribo-TISH") {
    orfs <- get_tish_orfs(file_list, gtf_exon, ncores)
  } else if (orf_method == "Ribotricer") {
    tricer_ref <- get_all_tricer_orfs(tricer_allorf_file)
    orfs <- prepare_ribotricer(file_list, tricer_ref, ncores)
    orfs <- mclapply(orfs, function(orf) {
      orf <- add_stop_codon(orf, gtf_exon)
      return(orf)
    }, mc.cores = ncores)
  }
  orfs <- mclapply(orfs, function(orf) {
    orf_cats <- orf %>% as.data.frame() %>% 
        mutate(transcript_id = as.character(str_match(group_name, 
                                                      "ENST\\d+|MSTRG[.]\\d+[.]\\d+"))) %>% 
        group_by(group_name, transcript_id, strand) %>% 
        summarise(orf_sta = min(start), orf_end = max(end), 
                  .groups = "keep") %>% annotate_categories(gtf_ref) %>% 
      ungroup() %>% select(group_name, orf_category)
    mcols(orf) <- mcols(orf) %>% as.data.frame() %>% 
      left_join(orf_cats, c("orf_id" = "group_name"))
    return(orf)
  }, mc.cores = ncores)
  return(orfs)
}

annotate_categories <- function(orf_seq, gtf_cds) {
  # Determine ORF categories based on transcript and GTF annotation
  #
  # Parameters:
  #   orf_seq - dataframe, ORFs with transcript ids and genomic start and end 
  #     coordinates
  #   gtf_cds - dataframe, gtf based reference with transcript ids and genomic 
  #     start and end coordinates
  #   ncores - int, number of cores to use.
  #
  # Returns:
  #   orf_cat_compare - dataframe, input dataframe with added orf_category 
  #     column
  orf_cat_compare <- orf_seq %>% left_join(gtf_cds, "transcript_id") %>% 
    mutate(orf_category = case_when(
      gtf_sta == orf_sta & gtf_end == orf_end ~ "ORF_annotated",
      gtf_sta > orf_sta & gtf_end == orf_end & strand == "+" ~ "N_extension",
      gtf_sta > orf_sta & gtf_end == orf_end & strand == "-" ~ "C_extension",
      gtf_sta < orf_sta & gtf_end == orf_end & strand == "+" ~ "N_truncation",
      gtf_sta < orf_sta & gtf_end == orf_end & strand == "-" ~ "C_truncation",
      gtf_end < orf_sta & gtf_end < orf_end & strand == "+" ~ "dORF",
      gtf_end < orf_sta & gtf_end < orf_end & strand == "-" ~ "uORF",
      gtf_sta < orf_sta & gtf_end < orf_end & strand == "+" ~ "overl_dORF",
      gtf_sta < orf_sta & gtf_end < orf_end & strand == "-" ~ "overl_uORF",
      gtf_sta > orf_sta & gtf_sta > orf_end & strand == "+" ~ "uORF",
      gtf_sta > orf_sta & gtf_sta > orf_end & strand == "-" ~ "dORF",
      gtf_sta > orf_sta & gtf_end > orf_end & strand == "+" ~ "overl_uORF",
      gtf_sta > orf_sta & gtf_end > orf_end & strand == "-" ~ "overl_dORF",
      gtf_sta < orf_sta & gtf_end > orf_end ~ "nested_ORF",
      gtf_sta > orf_sta & gtf_end < orf_end ~ "NC_extension",
      gtf_sta == orf_sta & gtf_end < orf_end & strand == "+" ~ "C_extension",
      gtf_sta == orf_sta & gtf_end < orf_end & strand == "-" ~ "N_extension",
      gtf_sta == orf_sta & gtf_end > orf_end & strand == "+" ~ "C_truncation",
      gtf_sta == orf_sta & gtf_end > orf_end & strand == "-" ~ "N_truncation",
      TRUE ~ "novel"
    )
  )
}

add_stop_codon <- function(orf_granges, gtf_exon) {
  # Extend all input ORFs to include the stop codon
  #
  # Parameters:
  #   orf_granges - GRangesList, ORFs to which a stop codon needs to be added
  #   gtf_exon - dataframe, gtf based reference annotation
  #
  # Returns:
  #   orf_shift - GRangesList, same as input, but with stopcodon included
  m_cols <- mcols(orf_granges)
  orf_df <- orf_granges %>% as.data.frame() %>% 
    mutate(transcript_id = as.character(str_match(group_name, "ENST\\d+|MSTRG[.]\\d+[.]\\d+"))) %>% 
    dplyr::rename(orf_start = start, orf_end = end, orf_strand = strand)
  orf_plus <- orf_df %>% filter(orf_strand == "+") %>% 
    inner_join(gtf_exon, "transcript_id", relationship = "many-to-many") %>% 
    filter(start <= orf_start, end >= orf_end) %>% 
    group_by(group_name) %>% 
    mutate(new_end = ifelse(orf_end == max(orf_end), orf_end + 3, orf_end),
           dif = ifelse(new_end > end, new_end - end, NA))
  orf_plus_extra <- orf_plus %>% filter(!is.na(dif)) %>% 
    left_join(gtf_exon, c("transcript_id", "strand"), 
              relationship = "many-to-many") %>% 
    filter(start.y > orf_start) %>% 
    group_by(group_name, dif, seqnames, strand) %>% 
    summarise(orf_start = min(start.y), .groups = "keep") %>% 
    mutate(new_end = orf_start + dif - 1) %>% ungroup() %>% select(-dif)
  orf_plus_shift <- orf_plus %>% 
    select(group_name, seqnames, strand, orf_start, new_end) %>% 
    rbind(orf_plus_extra) %>% dplyr::rename(start = orf_start, end = new_end)
  
  orf_minus <- orf_df %>% filter(orf_strand == "-") %>% 
    inner_join(gtf_exon, "transcript_id", relationship = "many-to-many") %>% 
    filter(start <= orf_start, end >= orf_end) %>% 
    group_by(group_name) %>% 
    mutate(new_start = ifelse(orf_start == min(orf_start), orf_start - 3, 
                              orf_start),
           dif = ifelse(new_start < start, start - new_start, NA))
  orf_minus_extra <- orf_minus %>% filter(!is.na(dif)) %>% 
    left_join(gtf_exon, c("transcript_id", "strand"), 
              relationship = "many-to-many") %>% 
    filter(end.y < orf_end) %>% 
    group_by(group_name, dif, seqnames, strand) %>% 
    summarise(orf_end = max(end.y), .groups = "keep") %>% 
    mutate(new_start = orf_end - dif + 1) %>% ungroup() %>% select(-dif)
  orf_minus_shift <- orf_minus %>% 
    select(group_name, seqnames, strand, new_start, orf_end) %>% 
    rbind(orf_minus_extra) %>% dplyr::rename(start = new_start, end = orf_end)
  orf_shift <- rbind(orf_plus_shift, orf_minus_shift) %>% distinct() %>% 
    makeGRangesListFromDataFrame(split.field = "group_name") 
  orf_shift <- orf_shift[m_cols$orf_id]
  mcols(orf_shift) <- m_cols
  return(orf_shift)
}

get_correct_geneid <- function(t_id, g_id) {
  # Retrieves correct gene id if there are multiple gene ids
  #
  # Parameters:
  #   t_id - string, the transcript id
  #   g_id - string, gene id determined by PRICE
  #   gid_tid_table - data frame, maps transcript ids to gene ids
  #
  # Returns:
  #   correct_gene_id - string, g_id if it consisted of a single gene id, the 
  #                     gene id found based on the transcript id if it consisted
  #                     of multiple gene ids
  split_id <- str_split(g_id, "_", simplify = T)
  if (length(split_id) == 1) {
    return(g_id)
  }
  correct_gene_id <- tid_gid$gene_id[which(tid_gid$transcript_id == t_id)]
  return(correct_gene_id)
}

prepare_orfquant <- function(files) {
  # Read ORFquant files into GRangesLists
  #
  # Parameters:
  #   files - vector, contains the file paths to the ORFquant RData files
  #
  # Returns:
  #   orfs_list - list, contains GRangesLists with the ORFquant orfs
  orfs_list <- vector(mode = "list", length = length(files))
  for (i in 1:length(files)) {
    load(files[i])
    orfs <- split(ORFquant_results$ORFs_gen, names(ORFquant_results$ORFs_gen))
    gene_ids <- data.frame(transcript_id = 
                             ORFquant_results$ORFs_tx$transcript_id, 
                           gene_id = ORFquant_results$ORFs_tx$gene_id, 
                           orf_type = ORFquant_results$ORFs_tx$ORF_category_Tx,
                           orf_id = names(ORFquant_results$ORFs_tx),
                           p_sites = ORFquant_results$ORFs_tx$P_sites_raw,
                           p_value = ORFquant_results$ORFs_tx$pval)
    orf_ids <- data.frame(orf_id = names(orfs))
    mcols(orfs) <- left_join(orf_ids, gene_ids, "orf_id")
    orfs_list[[i]] <- orfs
  }
  return(orfs_list)
}

get_combined_orfs <- function(orfs_list, files) {
  # Combine a list of GRangesLists into a single GRangeList
  #
  # Parameters:
  #   orfs_list - list, contains GRangesLists with the ORFs
  #   files - vector, contains the file paths to the input files used to create
  #                   orfs_list
  #
  # Returns:
  #   orfs_list_combined - GRangesLists, contains GRangesLists with the ORFquant 
  #                        orfs. Has an extra metacolumn source with the index
  #                        indicating from which sample the ORF was
  for (i in 1:length(orfs_list)) {
    orf_mcols <- mcols(orfs_list[[i]])
    orf_mcols$source <- rep(i, NROW(orfs_list[[i]]))
    mcols(orfs_list[[i]]) <- orf_mcols
  }
  orfs_list_unique <- lapply(orfs_list, remove_duplicates)
  orfs_list_combined <- do.call("c", orfs_list_unique)
  duplicate_hits <- findOverlaps(orfs_list_combined, orfs_list_combined, type = "equal")
  duplicate_count <- table(table(duplicate_hits@from))
  duplicate_count <- duplicate_count / as.integer(names(duplicate_count))
  sample_counts <- data.frame(sample=sapply(files, basename, USE.NAMES = FALSE),
                              all_orfs=sapply(orfs_list, NROW), 
                              unique_orfs=sapply(orfs_list_unique, NROW))
  print(sample_counts)
  print(paste0("All ORFs total: ", sum(sample_counts$all_orfs)))
  print(paste0("Unique ORFs total: ", sum(sample_counts$unique_orfs)))
  print("Counts of ORFs occuring in multiple samples:")
  print(duplicate_count)
  print(paste0("Sum of unique ORFs: ", sum(duplicate_count)))
  return(orfs_list_combined)
}

remove_duplicates <- function(g_range_list) {
  # Removes duplicates from the given GRangesList based on exact overlap
  #
  # Parameters:
  #   g_range_list - GRangesList
  #
  # Returns:
  #   g_range_unique - GRangesList, same as input, but without duplicates
  duplicate_hits <- findOverlaps(g_range_list, g_range_list, type = "equal")
  if (length(duplicate_hits) == length(g_range_list)) {
    return(g_range_list)
  }
  g_range_unique <- g_range_list[-duplicate_hits@to[
    duplicated(duplicate_hits@from)]]
  return(g_range_unique)
}

create_overlap_matrix <- function(com_orfs, thd, overlap_list, orf_lengths, 
                                  orf_sources, sample_num) {
  # Compute overlap between a set of ORFs for a give required overlap fraction
  #
  # Parameters:
  #   com_orfs - GRangesLists, the ORFs, including a source meta column
  #   thd - double, the required overlap fraction
  #   overlap_list - list, maps the indexes of overlapping ORFs
  #   orf_lengths - vector, the length of each ORF in com_orfs
  #   orf_sources - vector, the source index of each ORF in com_orfs
  #   sample_num - int, the total number of samples
  #
  # Returns:
  #   orfs_list_combined - GRangesLists, contains GRangesLists with the ORFquant 
  #                        orfs
  overlap_matrix <- matrix(0, ncol = sample_num + 1, nrow = NROW(com_orfs))
  for (i in 1:NROW(com_orfs)) {
    source_idx <- orf_sources[i]
    pos_parents <- overlap_list[[i]]
    overlap_matrix[i,source_idx] <- overlap_matrix[i,source_idx] + 1
    if (length(pos_parents) == 1) {
      next
    }
    child_width <- orf_lengths[i]
    len_parents <- sapply(pos_parents, function(x) child_width/orf_lengths[x])
    parent <- pos_parents[min(which(len_parents == 
                                      min(len_parents[len_parents >= thd])))]
    if (parent == i) {
      next
    }
    overlap_matrix[i,sample_num + 1] <- parent
    while (parent > 0) {
      overlap_matrix[parent, source_idx] = 
        overlap_matrix[parent, source_idx] + 1
      parent <- overlap_matrix[parent, sample_num + 1]
    }
  }
  parent_overlap_counts <- overlap_matrix[overlap_matrix[,sample_num + 1] == 
                                            0,1:sample_num]
  duplicate_counts <- table(rowSums(parent_overlap_counts > 0))
  print(paste0("Counts of ORFs occuring in multiple when requiring ",
               "minimum overlap of ",thd, ":"))
  print(duplicate_counts)
  print(paste0("Sum of unique ORFs: ", sum(duplicate_counts)))
  return(overlap_matrix)
}

get_partial_overlap <- function(com_orfs, thds) {
  # Computes full or partial overlap between a set of ORFs
  #
  # Parameters:
  #   com_orfs - GRangesLists, the ORFs, including a source meta column
  #   thds - vector, the required overlap fractions for which to compute the 
  #          overlap
  #
  # Returns:
  #   all_overlap - matrix, contains per threshold-source combination a column
  #                 named threshold_source and per threshold a column named
  #                 threshold_parent. Number of rows is equal to the length of
  #                 com_orfs. Per ORF and threshold it is indicated whether it
  #                 has an overlapping parent by an index in the parent column.
  #                 The threshold_source columns indicate which samples contain
  #                 an overlapping ORF.
  within_overlap <- findOverlaps(com_orfs, com_orfs, type = "within")
  minus_end_overlap <- findOverlaps(com_orfs[strand(com_orfs) == "-"], 
                                    com_orfs, type = "start")
  plus_end_overlap <- findOverlaps(com_orfs[strand(com_orfs) == "+"], 
                                   com_orfs, type = "end")
  within_overlap_list <- split(within_overlap@to, within_overlap@from)
  minus_end_overlap_list <- split(to(minus_end_overlap), 
                                  from(minus_end_overlap))
  plus_end_overlap_list <- split(to(plus_end_overlap), from(plus_end_overlap))
  unsorted_end_overlap_list <- c(minus_end_overlap_list, plus_end_overlap_list)
  end_overlap_list <- vector(mode = "list", 
                             length = length(unsorted_end_overlap_list))
  for (i in 1:length(unsorted_end_overlap_list)) {
    row_data <- unsorted_end_overlap_list[[i]]
    for (j in 1:length(row_data)) {
      end_overlap_list[[row_data[j]]] <- row_data
    }
  }
  overlap_list <- mapply(intersect, within_overlap_list, end_overlap_list)
  sample_num <- max(mcols(com_orfs)$source)
  all_overlap <- matrix(0, ncol = (sample_num+1) * length(thds), 
                        nrow = NROW(com_orfs))
  col_names_a <- rep(thds, each = sample_num+1)
  col_names_b <- rep.int(c(1:sample_num, "parent"), length(thds))
  colnames(all_overlap) <- paste(col_names_a, col_names_b, sep="_")
  orf_lengths <- sum(width(com_orfs))
  orf_sources <- mcols(com_orfs)$source
  for (i in 1:length(thds)) {
    overlap_matrix <- create_overlap_matrix(com_orfs, thds[i], overlap_list, 
                                            orf_lengths, orf_sources, sample_num)
    all_overlap[,((sample_num+1)*(i-1)+1):((sample_num+1)*i)] <- overlap_matrix
  }
  return(all_overlap)
}

get_overlap_counts <- function(m_overlap, thds, pred_source) {
  # Get count of shared ORFs between 1, 3 & 6 samples for each threshold 
  #
  # Parameters:
  #   m_overlap - matrix, created by get_partial_overlap function
  #   thds - vector, the required overlap fractions for which to compute the 
  #          overlap
  #   pred_source - string, name of the orfcaller
  #
  # Returns:
  #   graph_data - dataframe, contains for each combination of overlap fraction 
  #     and whether an ORF needs to be shared between 1, 3 & 6 samples the total
  #     number of unique orfs
  sample_num <- ncol(m_overlap)/length(thds) - 1
  graph_data <- matrix(0, ncol = 6, nrow = length(thds))
  colnames(graph_data) <- c("x", "All", "Half","Common","Viable", "Source")
  graph_data <- as.data.frame(graph_data)
  graph_data[,"x"] <- thds
  for (i in 1:length(thds)) {
    overlap_cols <- sapply(1:sample_num, function(x) paste(thds[i], x, sep="_"))
    parent_col <- paste(thds[i], "parent", sep="_")
    overlap_counts <- m_overlap[m_overlap[,parent_col] == 0,overlap_cols]
    duplicate_counts <- table(rowSums(overlap_counts > 0))
    total_orf <- sum(duplicate_counts)
    common_orf <- sum(duplicate_counts[sample_num])
    half_orf <- sum(duplicate_counts[round(sample_num/2):sample_num])
    viable_orf <- sum(duplicate_counts[3:sample_num])
    graph_data[i,"All"] <- total_orf
    graph_data[i, "Common"] <- common_orf
    graph_data[i, "Half"] <- half_orf
    graph_data[i, "Viable"] <- viable_orf
    
    
  }
  graph_data <- graph_data %>% pivot_longer(cols = c("All", "Half","Common","Viable")) %>% 
    mutate(Source = pred_source)
  return(graph_data)
}

filter_orfs <- function(m_overlap, m_orfs) {
  # Determine for each unique ORFs between how many samples it is shared
  # Parameters:
  #   m_overlap - matrix, created by get_partial_overlap function
  #   m_orfs - GRangesList, all ORFs from an ORF caller
  #
  # Returns:
  #   orfs_filtered - dataframe, rownames are ORF ids, the column indicates
  #     between how many samples the ORF was shared
  orfs_filtered <- m_overlap %>% 
    as.data.frame(row.names = names(m_orfs)) %>% 
    filter(`1_parent` == 0) %>% 
    dplyr::select(starts_with("1_") & !ends_with("parent")) %>% 
    apply(1, function(x) sum(x>0)) %>% 
    as.data.frame()
  names(orfs_filtered) <- "rowsum"
  return(orfs_filtered)
}

message(paste(Sys.time(),"Preparing ORF reference ..."))
gtf_data <- fread(gtf_file, skip = 5)
colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute")
tid_gid <- gtf_data %>% 
    mutate(transcript_id = as.character(str_match(attribute, "ENST\\d+|MSTRG[.]\\d+[.]\\d+")),
           gene_id = as.character(str_match(attribute, "ENST\\d+|MSTRG[.]\\d+[.]\\d+"))) %>% 
  select(transcript_id, gene_id) %>% 
  filter(!is.na(transcript_id)) %>% distinct()

gtf_exon <- gtf_data %>% 
  mutate(transcript_id = as.character(str_match(attribute, "ENST\\d+|MSTRG[.]\\d+[.]\\d+")),
         gene_id = as.character(str_match(attribute, "ENST\\d+|MSTRG[.]\\d+[.]\\d+"))) %>% 
  filter(feature == "exon") %>% select(seqname, start, end, strand, 
                                       transcript_id, gene_id)

lnc_genes <- gtf_data %>% 
  mutate(gene_id = as.character(str_match(attribute, "ENST\\d+|MSTRG[.]\\d+[.]\\d+"))) %>% 
  filter(feature == "gene") %>% 
  mutate(biotype = as.character(
    str_match(attribute,"(?<=gene_biotype \").*(?=\")"))) %>% 
  filter(biotype == "lncRNA") %>% 
  select(gene_id, biotype)

#Load GTF rannot file from ORFquant
load(gtf_rannot_file)
gtf_ref <- GTF_annotation$cds_txs %>% as.data.frame() %>% 
    group_by(group_name) %>% 
    summarise(gtf_sta = min(start), gtf_end = max(end)) %>% 
    dplyr::rename(transcript_id = group_name)

message(paste(Sys.time(),"Preparing separate ORFs ..."))

#orfquant_orfs <- prepare_orfs(orfquant_files, orf_methods[1], 2)

#save(orfquant_orfs, 
#     file = paste0(savedir, "/orfs_separate.RData"))

message(paste(Sys.time(),"Preparing combined ORFs ..."))

load(paste0(savedir, "/orfs_separate.RData"))

orfquant_orfs_combined <- get_combined_orfs(orfquant_orfs, orfquant_files)

save(orfquant_orfs_combined, 
     file = paste0(savedir, "/orfs_combined.RData"))

message(paste(Sys.time(),"Creating overlap matrix ..."))

orfquant_overlap <- get_partial_overlap(orfquant_orfs_combined, thresholds)

save(orfquant_overlap, 
     file = paste0(savedir, "/orfs_overlap.RData"))

message(paste(Sys.time(),"ORF analysis finished"))
