#Load required packages
library(phyloseq)
library(ggplot2)
library(patchwork)
library(vegan)

#Create ouput directory
dir.create("3-results/Postprocessing-diversity/figures", recursive = TRUE)

load(file = "3-results/DADA2/phyloseq.RData")
samdf <- as.data.frame(sample_data(ps))
##############################################
#Alpha Diversity: Vaginal vs. Endometrial brushing
alpha_diversity <- estimate_richness(ps, measures = c("Observed", "Shannon"))
alpha <- merge(alpha_diversity, samdf, by = "row.names")

##Alpha diversity: Plots
plot_alpha_1 <- ggplot(alpha, aes(x = SampleType, y = Observed, fill=SampleType)) +
  theme_bw() +
  geom_boxplot() +
  labs(x="Sample Type",
       y="Observed Richness",
       title="Observed Richness")
plot_alpha_1

plot_alpha_2 <- ggplot(alpha, aes(x = SampleType, y = Shannon, fill=SampleType)) +
  theme_bw() +
  geom_boxplot() +
  labs(x="Sample Type",
       y="Shannon Index",
       title="Shannon Diversity Index")
plot_alpha_2

alpha_plots <- (plot_alpha_1+plot_alpha_2)

#save alpha diversity boxplots
ggsave(
  "3-results/Postprocessing-diversity/figures/alpha_diversity_boxplot.png",
  alpha_plots,
  width = 10,
  height = 5,
  dpi = 300
)

##Statistical test
##Test normality for the Vaginal Swab group
shapiro.test(alpha$Shannon[alpha$SampleType == "Vaginal swab"])

##Test normality for the Endometrial Brushing group
shapiro.test(alpha$Shannon[alpha$SampleType == "Endometrial brushing"])

#Statistical test: Paired Samples 
alpha_wide <- reshape(
  alpha,
  idvar = "SubjectID",
  timevar = "SampleType",
  direction = "wide"
)

colnames(alpha_wide) <- make.names(colnames(alpha_wide))

alpha_paired_test <- wilcox.test(
  alpha_wide$'Shannon.Vaginal.swab',
  alpha_wide$'Shannon.Endometrial.brushing',
  paired = TRUE,
  exact = FALSE
)

alpha_paired_test

#####################################################
#Beta diversity using Bray-Curtis metric
##Bray-Curtis 
bc <- phyloseq::distance(ps, method="bray")
as.matrix(bc)[1:5, 1:5]

##Beta diversity: Plots
bray_ordination <- ordinate(ps, method = "PCoA", distance = "bray")

beta_bc_plot <- plot_ordination(ps, 
                                bray_ordination, 
                                color="SampleType", 
                                title="PCoA by Sample Type - Bray Curtis") + stat_ellipse( )
beta_bc_plot

#save beta diversity PCOA
ggsave(
  "3-results/Postprocessing-diversity/figures/bray_curtis_pcoa.png",
  beta_bc_plot,
  width = 7,
  height = 5,
  dpi = 300
)

##Statistical test: PERMANOVA
#Paired
set.seed(42)
paired <- adonis2(
  bc ~ SampleType,
  data = samdf,
  permutations = 999,
  strata = samdf$SubjectID
)

print(paired)


