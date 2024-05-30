library(phyloseq)
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)

###############################################################################################################
### set workding directory
setwd("/home/folder")
options(scipen=999) #no scientific notation
merged_metagenomes <- import_biom("rumen_microbiome_dna_seq.biom")
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)

colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
###############################################################################################################
#import the metadata file
metadata <- read.csv("metadata.tsv", sep = "\t", header = T, row.names = 1)

id_addition <- function(input_matrix) {
  for (i in seq_len(nrow(input_matrix))) {
    input_matrix[i, "unique_id"] <- paste0(unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[3], unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[4])
  }
  return(input_matrix)
}

metadata <- id_addition(metadata)
sampledata <- sample_data(metadata)


merged_rumen_kit = merge_phyloseq(merged_metagenomes, sampledata)

#subset the data to dna extraction
merged_rumen_seq_kit = subset_samples(merged_rumen_kit, !(extraction_kit == "Puregene" | extraction_kit == "Dneasy"))
merged_rumen_dna_kit = subset_samples(merged_rumen_kit, location=="Australia" & !(sequencing_kit == "Rapid_kit"))

mr.rarefied_seq_kit = rarefy_even_depth(merged_rumen_seq_kit, rngseed=7127, sample.size=min(sample_sums(merged_rumen_seq_kit)), replace=F)
mr.rarefied_dna_kit = rarefy_even_depth(merged_rumen_dna_kit, rngseed=7127, sample.size=min(sample_sums(merged_rumen_dna_kit)), replace=F)
mr.rarefied_basecall = rarefy_even_depth(merged_rumen_kit, rngseed=7127, sample.size=min(sample_sums(merged_rumen_kit)), replace=F)

###############################################################################################################

set.seed(7127)

###############################################################################################################
pcoa_caulculate <- function(p_object) {
  dist_bc = phyloseq::distance(p_object, method="bray")
  
  #beta diversity using principal coordinate analysis (PCoA) and Bray-Curtis distance
  pcoa_bc = ordinate(p_object, method = "PCoA", distance = dist_bc)
  return(pcoa_bc)
}

pcoa_seq_kit <- pcoa_caulculate(mr.rarefied_seq_kit)
pcoa_dna_kit <- pcoa_caulculate(mr.rarefied_dna_kit)

plot_ordination(mr.rarefied_seq_kit, pcoa_seq_kit, color = "sequencing_kit", shape = "animal_id") + 
  geom_point(aes(size = basecall_mode)) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  scale_shape_manual("Animal ID", breaks=c("AUS01", "SPA01", "SPA02", "SPA03", "SPA04"),
                     values = c(16, 17, 15, 8, 3)) +
  scale_size_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                    values = c(2, 4, 6)) +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.tag = element_text(face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))

plot_ordination(mr.rarefied_dna_kit, pcoa_dna_kit, color = "extraction_kit", axes = 1:2) + 
  geom_point(aes(color = extraction_kit, fill=extraction_kit), alpha = 0.5) +
  geom_point(aes(size = basecall_mode)) +
  scale_size_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                    values = c(2, 4, 6)) +
  scale_fill_manual("Extraction Kit", breaks=c("Dneasy", "PowerFecal", "Puregene"),
                      values = c("#549AAB", "#a73c5a", "#ff7954")) +
  scale_colour_manual("Extraction Kit", breaks=c("Dneasy", "PowerFecal", "Puregene"),
                      values = c("#549AAB", "#a73c5a", "#ff7954")) +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.tag = element_text(face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))
