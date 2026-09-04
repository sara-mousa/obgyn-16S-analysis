#Convert ASV table to BIOM format
library(biomformat)
################################
#Read the ASV table
asv <- read.csv("3-results/DADA2/asv_table.csv",row.names = 1)

#Transpose the table so ASVs are rows and samples are columns
asv <- as.data.frame(t(asv))

#Save Transposed ASV table as csv
write.csv(asv, "3-results/transposed_ASVs.csv", row.names = TRUE)

#Convert to BIOM format
asv_biom <- make_biom(asv)

#Write BIOM file
write_biom(asv_biom, "3-results/otu_table.biom")

