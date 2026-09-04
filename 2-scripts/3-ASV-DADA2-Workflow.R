#DADA2 processing and phyloseq construction
#DADA2 processing follows the workflow described in the
# "DADA2 Pipeline Tutorial (1.16)" by Callahan et al.
#############################################
#Load required packages
library(readxl)
library(dplyr)
library(dada2)
library(ggplot2)
library(phyloseq)
library(Biostrings)
#############################################
#Create output directory 
dir.create("3-results/DADA2/figures", recursive = TRUE)

#Import METADATA of the selected samples
run_data <-read.csv("1-data/SraRunTable.csv")
metadata <- run_data[, c("Run","Sample.Name")]

metadata$SampleType <- ifelse(
  grepl("DNA-V", metadata$Sample.Name),
  "Vaginal swab",
  ifelse(
    grepl("DNA-U", metadata$Sample.Name),
    "Endometrial brushing",
    NA
  )
)

metadata$SubjectID <- sapply(strsplit(metadata$Sample.Name, "-"), function(x) x[3])

supp_data <- read_xlsx("1-data/Supp-table-1.xlsx", skip=3)
supp_data <- supp_data[, c(1:4)]
colnames(supp_data) <- c("SubjectID", "InfertilityDx", "Age", "BMI")


#Exclude subjects with missing paired samples or with missing demographic/clinical data
metadata <- left_join(
  metadata,
  supp_data,
  by = c("SubjectID")
)

paired_metadata <- metadata %>%
  group_by(SubjectID) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  arrange(SubjectID, SampleType)

paired_metadata <- filter(paired_metadata, !SubjectID %in% c("G083", "G113", "G132"))

write.csv(paired_metadata, file = "1-data/paired-metadata.csv")
#######################################################################
#Run the DADA2 Pipeline to produce ASV table 

path <- "1-data/fastq"
list.files(path)

#Forward and reverse fastq filenames have format
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

#Extract sample names as accession number
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

#Samples to process
keep <- sample.names %in% paired_metadata$Run
fnFs <- fnFs[keep]
fnRs <- fnRs[keep]
sample.names <- sample.names[keep]

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
tail(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#Inspect the merger data
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Save essential DADA2 outputs (FilterAndTrim + errF/errR + mergers)
save(out, errF, errR, mergers,file = "3-results/DADA2/dada2_output.RData")

#Assign taxonomy using SILVA NR99 v138.1
taxa <- assignTaxonomy(seqtab.nochim,
                       "1-data/taxa-database/silva_nr99_v138.1_train_set.fa",
                       multithread=TRUE)
taxa.print <- taxa

#Removing sequence rownames
rownames(taxa.print) <- NULL
head(taxa.print)

#Construct phyloseq
samdf <- paired_metadata
str(samdf)
samdf <- as.data.frame(samdf)
rownames(samdf) <- paired_metadata$Run
samdf<- samdf[rownames(seqtab.nochim), , drop = FALSE]
head(samdf)


ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
  sample_data(samdf), 
  tax_table(taxa)
)

print(ps)

###################################################
#Store reference sequences in the phyloseq refseq
#Extract sequences from the column names
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

#Merge the reference sequences into the object
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

##################################################
#Extract the DNAStringSet from your phyloseq object
asv_sequences <- refseq(ps)

#Write the sequences to a FASTA file
writeXStringSet(asv_sequences, filepath = "3-results/DADA2/asv_sequences.fasta", format = "fasta")
write.csv(as.data.frame(otu_table(ps)), file = "3-results/DADA2/asv_table.csv")
length(asv_sequences)

##########################################################################
#Identify the top 20 most abundant ASVs
top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]

#Transform absolute read counts to relative abundance
ps_abun <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))

#Prune the object down to contain only top 20 ASVs
ps_top20 <- prune_taxa(top20, ps_abun)

#Melt into a data frame
ps_df <- psmelt(ps_top20)

#Average the relative abundance for each SampleType
ps_summary <- ps_df %>%
  group_by(SampleType, Genus) %>%
  summarize(Mean_Abundance = mean(Abundance), .groups = "drop")

#Build the plot
top20_asv <- ggplot(ps_summary, aes(x = SampleType, y = Mean_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option = "mako") +
  labs(x = "Sample Type",
       y = "Mean Relative Abundance",
       title = "Top 20 ASVs: Genus-Level Relative Abundance")

top20_asv

#Save barplot
ggsave(
  "3-results/DADA2/figures/top20_ASVs_Genus.png",
  top20_asv,
  width = 10,
  height = 5,
  dpi = 300
)

#####################################################
#Extract the taxonomy table from the phyloseq object
taxa_matrix <- as(tax_table(ps), "matrix")
taxa_df <- as.data.frame(taxa_matrix, stringsAsFactors = FALSE)
taxa_df$ASV <- rownames(taxa_df)
write.csv(taxa_df[, c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
          file = "3-results/DADA2/taxonomy_table.csv",
          row.names = FALSE)

save(ps, file = "3-results/DADA2/phyloseq.RData")

