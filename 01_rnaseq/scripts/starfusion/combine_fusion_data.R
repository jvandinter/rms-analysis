library(magrittr)
library(ggplot2)

# Check Star-Fusion output
basedir="/hpc/pmc_vanheesch/projects/Jip"

setwd(paste(
  basedir,
  "rms_analysis",
  sep = "/")
)

rms_meta <-
  read.delim(paste(
    basedir,
    "custom_transcriptomes/20221020_JD_quant_tumor_cohorts/analysis/RMS/metadata/metadata_annot_RMS.txt",
    sep = "/"),
    sep = ";"
  )

sj_meta <-
  read.delim(paste(
    basedir,
    "/custom_transcriptomes/20221020_JD_quant_tumor_cohorts/analysis/metadata/stjude_metadata.txt",
    sep = "/")
  )

fusion_files <-
  list.files(paste(
    basedir,
    "rms_analysis/01_rnaseq/analysis/starfusion",
    sep = "/"),
    full.names = T,
    recursive = T,
    pattern = "star-fusion.fusion_predictions.tsv"
  )
sample_names <-
  gsub(paste(
    basedir,
    "rms_analysis/01_rnaseq/analysis/starfusion",
    sep = "/"),
    "",
    gsub("/star-fusion.fusion_predictions.tsv", "", fusion_files)
  )

# Wrangle metadata
meta <- sj_meta %>%
  dplyr::filter(sample_name %in% sample_names &
                  !(sample_name %in% rms_meta$sample_id)) %>%
  dplyr::select(sample_name,
                attr_subtype_biomarkers,
                attr_oncotree_disease_code) %>%
  dplyr::rename(sample_id = sample_name,
                type = attr_subtype_biomarkers,
                condition = attr_oncotree_disease_code) %>%
  rbind(rms_meta[, c("sample_id", "type", "condition")]) %>%
  tidyr::replace_na(list(type = "Not Available")) %>%
  dplyr::mutate(type = ifelse(
    type == "FOXO1-PAX3",
    "PAX3-FOXO1",
    ifelse(type == "other", "None", type)
  ))

# Wrangle fusion tool output
fusion_df <- list()

fusion_list <- lapply(fusion_files, function(x) {
  sample_id = gsub(paste(
    basedir,
    "rms_analysis/01_rnaseq/analysis/starfusion",
    sep = "/"),
    "",
    gsub("/star-fusion.fusion_predictions.tsv", "", x)
  )
  print(paste0("running:", sample_id))
  df = read.delim(x)
  df = df[grepl(pattern = "PAX.--FOXO1", df$X.FusionName), c(
    "X.FusionName",
    "JunctionReadCount",
    "SpanningFragCount",
    "LeftBreakpoint",
    "RightBreakpoint"
  )]
  if (nrow(df) > 0) {
    df$sample_id = sample_id
    return(df)
  } else {
    df = data.frame(
      X.FusionName = NA,
      JunctionReadCount = 0,
      SpanningFragCount = 0,
      LeftBreakpoint = NA,
      RightBreakpoint = NA,
      sample_id = sample_id
    )
    return(df)
  }
  
})

# Combine metadata and fusion calls

fusion_df <- do.call(rbind, fusion_list) %>%
  dplyr::mutate(X.FusionName = gsub("--", "-", X.FusionName)) %>%
  tidyr::replace_na(list(X.FusionName = "None")) %>%
  dplyr::left_join(meta) %>%
  dplyr::mutate(
    fusion_check = ifelse(
      X.FusionName == type,
      "same",
      ifelse(
        X.FusionName %in% c("PAX3-FOXO1", "PAX7-FOXO1", "None") &
          type %in% c("Not Available", "None"),
        "info_gain",
        ifelse(condition == "ARMS", "not FP", "same")
      )
    ),
    condition_check = ifelse(
      X.FusionName %in% c("PAX3-FOXO1", "PAX7-FOXO1") &
        condition == "ARMS",
      "harmony", ifelse(X.FusionName == "None" & condition %in% c("ERMS","SCRMS","RMS"),"harmony",
      "mismatch")
    )
  ) %>%
  dplyr::mutate(fusion_check = factor(fusion_check, levels = c("same", "info_gain", "not FP")))

# Plot differences
p1 =ggplot(fusion_df, aes(fusion_check, fill = X.FusionName)) + 
  geom_bar(position = "dodge") + 
  theme_classic()

# Plot type overlap
p2 = ggplot(fusion_df, aes(condition_check, fill = X.FusionName)) +
  geom_bar(position = "dodge") +
  theme_classic()

ggsave(file = paste(basedir,"rms_analysis/01_rnaseq/results/starfusion/fusion_check.pdf", sep ="/"),
  plot = ggpubr::ggarrange(p1,p2, common.legend = T)
)

write.table(fusion_df, file = paste(basedir,"rms_analysis/01_rnaseq/results/starfusion/fusion_annotation.txt", sep ="/"),
            quote = F, row.names = F, sep = "\t")
