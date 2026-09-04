#Analysis and Visualization of Functional Prediction (PICRUST2) Output
#############################################
#Load packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(ANCOMBC)
library(dplyr)
library(pheatmap)
#######################################
#Load PICRUSt2 output and metadata
#######################################
#Create output directory 
dir.create("3-results/Functional-Pred-Analysis/figures", recursive = TRUE)

#Load KO predicted abundance matrix
ko_counts <- read.table("3-results/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv", 
                        header=TRUE, 
                        sep="\t", 
                        row.names=1, 
                        check.names=FALSE)

#Load KO descriptions
ko_descrip <- read.table("3-results/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv", 
                         header=TRUE, 
                         sep="\t", 
                         row.names=1,
                         check.names=FALSE, 
                         quote="")
#Keep only description column
ko_descrip <- ko_descrip[, 1, drop = FALSE]

#Check the description column name
colnames(ko_descrip)

dim(ko_counts)
dim(ko_descrip)

##################################################
#Load MetaCyc Pathways predicted abundance matrix
path_counts <- read.table("3-results/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv", 
                          header=TRUE, 
                          sep="\t", 
                          row.names=1, 
                          check.names=FALSE)

#Load MetaCyc Pathways descriptions
path_descrip <- read.table("3-results/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv", 
                           header=TRUE, 
                           sep="\t", 
                           row.names=1,
                           check.names=FALSE, 
                           quote="")

#Keep only description column
path_descrip <- path_descrip[, 1, drop = FALSE]

#Check the description column name
colnames(path_descrip)

dim(path_counts)
dim(path_descrip)
#############################################
#Load subjects metadata
metadata <- read.csv("1-data/paired-metadata.csv", 
                     header=TRUE, 
                     row.names=1)

rownames(metadata) <- metadata$Run
metadata <- metadata[,-1]
##########################################################
#Construct Phyloseq objects for KOs and MetaCyc pathways
##########################################################
#Construct phyloseq objects out of PICRUSt2 output for easier handling of data
#KO phyloseq
ko_otu  <- otu_table(as.matrix(ko_counts), taxa_are_rows = TRUE)
ko_tax  <- tax_table(as.matrix(ko_descrip))
ko_meta <- sample_data(metadata)
ko_ps <- phyloseq(ko_otu, ko_tax, ko_meta)

#MetaCyc pathways phyloseq
path_otu  <- otu_table(as.matrix(path_counts), taxa_are_rows = TRUE)
path_tax  <- tax_table(as.matrix(path_descrip))
path_meta <- sample_data(metadata)
path_ps <- phyloseq(path_otu, path_tax, path_meta)

###############################################################
#Beta diversity and Statistical Analysis (on KOs and Pathways)
###############################################################
#KO/Bray-Curtis Index
ko_bc <- phyloseq::distance(ko_ps, method="bray")
as.matrix(ko_bc)[1:5, 1:5]

#KO/Plot
bray_ordination <- ordinate(ko_ps, method = "PCoA", distance = "bray")

plot_ko <- plot_ordination(ko_ps, 
                           bray_ordination, 
                           color="SampleType", 
                           title="PCoA of KOs - Bray Curtis") +
  stat_ellipse( )

plot_ko

#Save KO beta diversity plot
ggsave(
  "3-results/Functional-Pred-Analysis/figures/PCoA_KOs.png",
  plot_ko,
  width = 10,
  height = 5,
  dpi = 300
)

#KO/Statistical test: PERMANOVA
#Match the order of metadata to ko_bc
metadata_ko <- metadata[labels(ko_bc), , drop = FALSE]
all(rownames(metadata_ko) == labels(ko_bc))

set.seed(42)
ko_stat <- adonis2(
  ko_bc ~ SampleType,
  data = metadata_ko,
  permutations = 999,
  strata = metadata_ko$SubjectID
)
print(ko_stat)

#Save KO PERMANOVA results
write.csv(
  as.data.frame(ko_stat),
  "3-results/Functional-Pred-Analysis/KO_PERMANOVA_results.csv"
)
#################################
#Pathways/Bray-Curtis Index
path_bc <- phyloseq::distance(path_ps, method="bray")
as.matrix(path_bc)[1:5, 1:5]

#Pathways/Plot
bray_ordination <- ordinate(path_ps, method = "PCoA", distance = "bray")

plot_path <- plot_ordination(path_ps, 
                             bray_ordination, 
                             color="SampleType", 
                             title="PCoA of MetaCyc Pathways - Bray Curtis") +
  stat_ellipse( )
plot_path

#Save MetaCyc pathways beta diversity plot
ggsave(
  "3-results/Functional-Pred-Analysis/figures/PCoA_MetaCyc_Pathways.png",
  plot_path,
  width = 10,
  height = 5,
  dpi = 300
)

#Pathways/Statistical test: PERMANOVA
#Match the order of metadata to path_bc
metadata_path <- metadata[labels(path_bc), , drop = FALSE]
all(rownames(metadata_path) == labels(path_bc))

set.seed(42)
path_stat <- adonis2(
  path_bc ~ SampleType,
  data = metadata_path,
  permutations = 999,
  strata = metadata_path$SubjectID
)
print(path_stat)

#Save MetaCyc pathways PERMANOVA results
write.csv(
  as.data.frame(path_stat),
  "3-results/Functional-Pred-Analysis/MetaCyc_PERMANOVA_results.csv"
)
#############################################
#Differential abundance analysis(ANCOM-BC2)
#############################################
#Differentially Adbundant MetaCyc Pathways 
pathways_ancom <- ancombc2(
  data = path_ps, 
  assay_name = "counts", 
  tax_level = "description",
  fix_formula = "SampleType",
  group = "SampleType",
  p_adj_method = "BH",
  prv_cut = 0.20,
  lib_cut = 0,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05, 
)

#Extract results from ANCOM-BC2
pathways_res <- pathways_ancom$res

#Save Differentially Abundant MetaCyc Pathways results
write.csv(
  pathways_res,
  "3-results/Functional-Pred-Analysis/MetaCyc_ANCOMBC2_results.csv",
  row.names = FALSE
)

#Extract pathways with valid lfc and q-values
pathways_res <- pathways_res[
  complete.cases(
    pathways_res$`lfc_SampleTypeVaginal swab`, 
    pathways_res$`q_SampleTypeVaginal swab`),
]

dim(pathways_res)

#Subset significant pathways only (adj-p value < 0.05)
sig_pathways <- pathways_res[
  pathways_res$`q_SampleTypeVaginal swab`< 0.05,
]

#Save Differentially Abundant MetaCyc Pathways results with adj-pvalue<0.05
write.csv(
  sig_pathways,
  "3-results/Functional-Pred-Analysis/MetaCyc_ANCOMBC2_sig.csv",
  row.names = FALSE
)

#Subset and order top 10 pathways enriched in each sample type
top20_pathways <- rbind(
  slice_max(sig_pathways, `lfc_SampleTypeVaginal swab`, n = 10),
  slice_min(sig_pathways, `lfc_SampleTypeVaginal swab`, n = 10)
) %>%
  mutate(
    Direction = ifelse(`lfc_SampleTypeVaginal swab` > 0, "Vaginal swab", "Endometrial brushing"),
)

top20_pathways$pathway <- path_descrip[
  as.character(top20_pathways$taxon), "description"
]
top20_pathways <- top20_pathways %>%
  mutate(pathway = reorder(pathway, `lfc_SampleTypeVaginal swab`))
#######################################
#Barplot: Top 20 differentially abundant MetaCyc pathways
#######################################
#Pathways/Plot1
#Bar plot
pathways_plot1<-ggplot(
  top20_pathways,
  aes(
    x = `lfc_SampleTypeVaginal swab`,
    y = pathway,
    fill = Direction
  )
) +
  geom_col() +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +
  labs(
    title = "Top 20 Differentially Abundant MetaCyc Pathways by Sample Type",
    x = "Log-Fold Change (ANCOM-BC2)",
    y = "MetaCyc Pathways",
    fill = "Direction"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.y = element_text(size = 9.5),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom"
  )

pathways_plot1

#Save barplot
ggsave(
  "3-results/Functional-Pred-Analysis/figures/Top20_MetaCyc_diff_abun_barplot.png",
  pathways_plot1,
  width = 10,
  height = 5,
  dpi = 300
)

#######################################
#Heatmap: Top 20 differentially abundant MetaCyc pathways
#######################################
#Pathways/Plot2
#top20_pathways contains: Top 10 Vaginal and Top 10 Endometrial brushing pathways
#Create a vector with top pathways from top20_pathways and extract their counts
matched_rows <- match(top20_pathways$taxon, rownames(path_counts))
heat_counts <- path_counts[matched_rows, , drop = FALSE]

#Set the rownames to readable descriptions in path_descrip
rownames(heat_counts) <- as.character(top20_pathways$pathway)

#Perfrom log transformation on heat matrix counts
heat_matrix <- log10(as.matrix(heat_counts) + 1)

#Align the metadata rows to match the heat matrix column order
ordered_metadata <- metadata[colnames(heat_matrix), , drop = FALSE]

#Create an annotation column using the ordered metadata
annotation_col <- data.frame(
  SampleType = ordered_metadata$SampleType,
  row.names = colnames(heat_matrix)
)
#Create an annotation row to clarify what the order of metacyc pathways represent
annotation_row <- data.frame(
  Enriched_In = top20_pathways$Direction, 
  row.names = rownames(heat_matrix)
)

#Match the colors in annotation col and row
matched_colors <- list(
  SampleType = c("Vaginal swab" = "#D660AC", "Endometrial brushing" = "#24A293"),
  Enriched_In = c("Vaginal swab" = "#D660AC", "Endometrial brushing" = "#24A293")
)

#Plot the heatmap
pathways_plot2 <- pheatmap(
  heat_matrix,
  scale = "row",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = matched_colors,
  clustering_method = "ward.D2",
  main = "Top 20 Differentially Abundant MetaCyc Pathways",
  fontsize_row = 8,
  show_colnames = FALSE, 
  cluster_rows = FALSE
)

#Save the heatmap
pathways_plot2 <- pheatmap(
  heat_matrix,
  scale = "row",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = matched_colors,
  main = "Top 20 Differentially Abundant MetaCyc Pathways",
  fontsize_row = 8,
  show_colnames = FALSE, 
  cluster_rows = FALSE,
  filename = "3-results/Functional-Pred-Analysis/figures/Top20_MetaCyc_diff_abun_heatmap.pdf", 
  width = 12,
  height = 8,
  border_color = NA
)

graphics.off()
