##########
# Load CDS annotations CHECK
# Subset CDS annotations on in-frame P-sites
##
# Load intORF genomic loci CHECK
# Subset intORF on in-frame P-sites CHECK
##
# Check overlap between intORF P-sites & CDS P-sites
##########

library(magrittr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

basedir = "/hpc/pmc_vanheesch/projects/jvandinter/rms_analysis"

rna_dir = paste(basedir,
                "01_rnaseq",sep="/")
ribo_dir = paste(basedir,
                 "02_riboseq", sep = "/")

txdb_loc <- paste0(rna_dir,"/analysis/rnaseq_pipeline/customannotation/",
                   "RMS_container/",
                   "RMS_full_novel_filtered_corrected.gtf_TxDb")

orfquant_results_location <- paste(ribo_dir,
                                   "analysis/ORFquant/RMS_merged_psites",
                                   "RMS_merged_psites_final_ORFquant_results",sep="/")

orf_annotation <- read.delim(paste(ribo_dir,
                        "results/orf_reannotation",
                        "RMS_ORFquant_table.csv",
                        sep = "/"), sep = ",") %>%
  dplyr::filter(orf_category_new == "nested_ORF")

# P-sites
psites = read.delim(paste(ribo_dir,
                                   "analysis/p_site_quantification/bedfiles",
                                   "RMS_merged_psites_Detected_ORFs.gtf_psites_p0.sorted.bed", sep = "/"), header = F) %>%
  dplyr::filter(V4 %in% orf_annotation$ORF_id_tr & V5 == "p0")

intorf_granges <- psites %>% 
  dplyr::rename(seqnames = V1,
                start = V2,
                stop = V3,
                orf_id = V4,
                codon_pos = V5,
                strand = V6)%>% 
  GenomicRanges::makeGRangesFromDataFrame()

names(intorf_granges) <- psites$V4

intorf_grl <- split(intorf_granges, names(intorf_granges))

# Load ORF definitions
rms_orfquant_orfs <- get(load(orfquant_results_location))
# Extract ORF genomic ranges
orfs_gen <- rms_orfquant_orfs$ORFs_gen
nested_gen <- orfs_gen[names(orfs_gen) %in% orf_annotation$ORF_id_tr]
orf_ranges <- split(nested_gen, names(nested_gen))

txdb <- AnnotationDbi::loadDb(txdb_loc)

# Extract CDS regions from txdb
cds_tx <- cdsBy(txdb, by = "tx", use.names = T)  # not used currently
cds_gene <- cdsBy(txdb, "gene")
cds_gene_unlist <- unlist(cds_gene)

# Find CDS relevant for intORF

overlaps <- GenomicRanges::findOverlaps(orf_ranges, cds_gene)

overlap_width <- sum(width(GenomicRanges::intersect(orf_ranges[queryHits(overlaps)], cds_gene[subjectHits(overlaps)])))

overlap_df <- data.frame(queryIdx = queryHits(overlaps), 
                         subjectIdx = subjectHits(overlaps),
                         overlapWidth = overlap_width)

max_overlaps <- overlap_df[order(overlap_df$queryIdx, -overlap_df$overlapWidth),]
max_overlaps <- max_overlaps[!duplicated(max_overlaps$queryIdx),]

query_idx <- max_overlaps$queryIdx
subject_idx <- max_overlaps$subjectIdx

selected_overlaps <- data.frame(
  queryHits = 1:length(orf_ranges),
  subjectHits = rep(NA, length(orf_ranges))
)

selected_overlaps$subjectHits[selected_overlaps$queryHits %in% query_idx] <- subject_idx

result_list <- GRangesList(rep(list(GRanges()), length(orf_ranges)))
names(result_list) <- names(orf_ranges)

non_na_indices <- !is.na(selected_overlaps$subjectHits)
result_list[selected_overlaps$queryHits[non_na_indices]] <- cds_gene[selected_overlaps$subjectHits[non_na_indices]]
no_overlap_idx <- lengths(result_list) == 0
no_overlap_names <- names(which(no_overlap_idx))

result_list[no_overlap_idx] <- GRangesList(pbapply::pblapply(no_overlap_names, function(name) {
  orf_parent_gene <- orf_annotation$gene_id[match(name, orf_annotation$ORF_id_tr)]  # You need orf_table here, which contains mappings between ORF IDs and parent gene IDs
  cds_parent_gene <- cds_gene_unlist[which(names(cds_gene_unlist) == orf_parent_gene)] # Turns out I didn't need to find nearest CDS regions using `nearest()`, I could just use the parent gene ID -> an ORF can't be a dORF or uORF if it's in a different gene
  return(cds_parent_gene)
}))

cds_matches_grl <- result_list
