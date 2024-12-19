###########################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
#
# R script to generate PCA plot from featureCounts counts
# matrix.
#
###########################################################################

# Load required packages --------------------------------------------------
message("Loading required libraries ...")
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
})

# Global arguments --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

wd = args[1]
count_file = args[2]
pdf_plot_name = args[3]

# Load the RNA-seq data
message("Loading RNA-seq counts ...")
featurecount_outfile <- read.table(file = count_file,
                                   header = T)
count_df <- featurecount_outfile[,(ncol(featurecount_outfile)-(ncol(featurecount_outfile) - 7)):ncol(featurecount_outfile)] 
rownames(count_df) <- featurecount_outfile[,1]

sample_info <- data.frame(sample_ID = colnames(count_df),
                          condition = 1:ncol(count_df))

# Create normalised DESEQ count matrix
dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = sample_info,
                              design = ~ condition)

vsd <- vst(dds, blind=T)

# Generate PCA
ntop <- 500
intgroup <- colnames(sample_info)[2]

rv <- rowVars(assay(vsd))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

intgroup.df <- as.data.frame(colData(vsd)[, intgroup, drop=FALSE])
  
# add the intgroup factors together to create a new grouping factor
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(vsd)[[intgroup]]
}

# assembly the data for the plot
d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(vsd))

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = intgroup)) + 
    geom_point(size = 4) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    theme_classic() +
    coord_fixed()

# Save the plot
ggsave(device = "pdf",
       filename = paste(pdf_plot_name,"pdf",sep = "."),
       width = 8,
       height = 8,
       path = paste(wd, "data/processed/", sep = "/"))

message("Creating PCA plot ... Done!")