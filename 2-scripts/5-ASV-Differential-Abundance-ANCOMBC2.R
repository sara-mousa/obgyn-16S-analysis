#ASVs Differential Abundance (ANCOM-BC2)
##########################################
#Load required packages 
library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(gt)
#######################################################
#Create output directory 
dir.create("3-results/ASV-ANCOMBC2")

#Load phyloseq object
load("3-results/phyloseq.RData")

#Set Endometrial brushing as the reference group
sample_data(ps)$SampleType <- factor(
  sample_data(ps)$SampleType,
  levels = c("Endometrial brushing", "Vaginal swab")
)

#######################################################
#ASVs Differential Abundance analysis using ANCOM-BC2
out <- ancombc2(
  data = ps,
  fix_formula = "SampleType",
  group = "SampleType",
  rand_formula = "(1|SubjectID)",
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 0,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05
)

#Inspect results columns
colnames(out$res)

#Extract results 
res <- out$res

#Extract significant differentially abundant ASVs
sig_asv <- res[
  res$`q_SampleTypeVaginal swab` < 0.05,
]

#Save ANCOM-BC2 results
write.csv(
  res,
  "3-results/ASV-ANCOMBC2/ANCOMBC2_results.csv",
  row.names = FALSE
)

#Save significant ASVs
write.csv(
  sig_asv,
  "3-results/ASV-ANCOMBC2/significant_ASVs.csv",
  row.names = FALSE
)

#LFC
tab_coef <- res[, c(
  "taxon",
  "lfc_(Intercept)",
  "lfc_SampleTypeVaginal swab"
)]

tab_html <- tab_coef %>%
  gt(
    rowname_col = "taxon",
    caption = "Coefficients from the Primary Result"
  ) %>%
  as_raw_html()

write(tab_html, "3-results/ASV-ANCOMBC2/dataCoeff.html")

#SE
tab_SE <- res[, c(
  "taxon",
  "se_(Intercept)",
  "se_SampleTypeVaginal swab"
)]

tab_html <- tab_SE %>%
  gt(
    rowname_col = "taxon",
    caption = "SE from the Primary Result"
  ) %>%
  as_raw_html()

write(tab_html, "3-results/ASV-ANCOMBC2/SEs.html")


#Test Statistics 
tab_w <- res[, c(
  "taxon",
  "W_(Intercept)",
  "W_SampleTypeVaginal swab"
)]

tab_html <- tab_w %>%
  gt(
    rowname_col = "taxon",
    caption = "Test Statistics from the Primary Result"
  ) %>%
  as_raw_html()

write(tab_html, "3-results/ASV-ANCOMBC2/Test-Statistics.html")

#P-value
tab_p <- res[, c(
  "taxon",
  "p_(Intercept)",
  "p_SampleTypeVaginal swab"
)]

tab_html <- tab_p %>%
  gt(
    rowname_col = "taxon",
    caption = "P-value from the Primary Result"
  ) %>%
  as_raw_html()

write(tab_html, "3-results/ASV-ANCOMBC2/p-value.html")

#q-value
tab_q <- res[, c(
  "taxon",
  "q_(Intercept)",
  "q_SampleTypeVaginal swab"
)]

tab_html <- tab_q %>%
  gt(
    rowname_col = "taxon",
    caption = " Adjusted p-value from the Primary Result"
  ) %>%
  as_raw_html()

write(tab_html, "3-results/ASV-ANCOMBC2/adj-p-value.html")

#Differentially abundant ASVs
tab_diff <- res[, c(
  "taxon",
  "diff_(Intercept)",
  "diff_SampleTypeVaginal swab"
)]

tab_html <- tab_diff %>%
  gt(
    rowname_col = "taxon",
    caption = " Differentially Abundant ASVs from the Primary Result"
  ) %>%
  as_raw_html()

write(tab_html, "3-results/ASV-ANCOMBC2/differentially-abundant-ASVs.html")

