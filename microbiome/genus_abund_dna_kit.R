library(phyloseq)
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(readr)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(deeptime)
library(gridExtra)

###############################################################################################################
### set workding directory
setwd("/home/folder")
options(scipen=999) #no scientific notation
merged_metagenomes <- import_biom("rumen_microbiome_dna_seq.biom")
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)

colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#import the metadata file
metadata <- read.csv("metadata.tsv", sep = "\t", header = T, row.names = 1)
sampledata <- sample_data(metadata)
merged_rumen_kit = merge_phyloseq(merged_metagenomes, sampledata)

#subset the data to dna extraction
merged_rumen_dna_kit = subset_samples(merged_rumen_kit, 
                                      location=="Australia" & sequencing_kit=="Ligation_kit" & basecall_mode=="SUP")
mr.rarefied_dna_kit = rarefy_even_depth(merged_rumen_dna_kit, rngseed=7127, sample.size=min(sample_sums(merged_rumen_dna_kit)), replace=F)

###############################################################################################################

set.seed(7127)

###############################################################################################################
# get rid of some empty taxonomy levels
# kingdom_level is the kngdom level for input
kingdom_subset <- function(p_object, k_name) {
  kingdom_subset_mr <- subset_taxa(p_object, Genus != "") %>%
    subset_taxa(Kingdom == k_name)
  return(kingdom_subset_mr)
}

###############################################################################################################
##### otu manipulation

otu_relative <- function(kingdom_subset_mr, metadata) {
  otu <- rownames(tax_table(kingdom_subset_mr))
  taxonomy <- cbind(otu, tax_table(kingdom_subset_mr)) %>%
    as.data.frame()
  otu_counts <- t(otu_table(kingdom_subset_mr))
  sample_id <- rownames(otu_counts)
  otu_counts <- cbind(sample_id, otu_counts) %>%
    as.data.frame()
  otu_counts <- pivot_longer(otu_counts, -sample_id, names_to="otu", values_to = "count")
  otu_counts$count <- as.numeric(otu_counts$count)
  
  rownames(metadata) <- NULL
  
  otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>% #joint the data by "Replicate_ID"
    inner_join(taxonomy, by="otu") %>%
    group_by(sample_id) %>%
    mutate(rel_abund = count / sum(count)) %>%
    ungroup() %>%
    select(-count) %>%
    pivot_longer(cols =  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "otu"), 
                 names_to = "level",
                 values_to = "taxon") %>%
    mutate(Extraction_kit = factor(extraction_kit, levels=c("Dneasy", "PowerFecal", "Puregene")))
  return(otu_rel_abund)
}

###############################################################################################################
#calculate the relative abundance
###############################################################################################################

relative_calculation <- function(otu_rel_abund) {
  genus_rel_abund <-  otu_rel_abund %>%
    filter(level=='Genus') %>%
    group_by(location, extraction_kit, sample_id, taxon) %>%
    summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
    ungroup()
  return(genus_rel_abund)
}

###############################################################################################################
#filter some low abundance taxons
###############################################################################################################

taxon_pool_check <- function(genus_rel_abund) {
  genus_pool_check <- genus_rel_abund %>%
    group_by(location, extraction_kit, taxon) %>%
    summarise(mean = mean(rel_abund), .groups = "drop") %>%
    group_by(taxon)
  return(genus_pool_check)
}

taxon_pool <- function(genus_rel_abund, threhold) {
  genus_pool <- genus_rel_abund %>%
    group_by(location, extraction_kit, taxon) %>%
    summarise(mean = mean(rel_abund), .groups = "drop") %>%
    group_by(taxon) %>%
    summarise(pool = max(mean) < threhold, 
              mean = mean(mean),
              .group="drop")
  return(genus_pool)
}

###############################################################################################################
#combine the data frame
###############################################################################################################

taxon_pool_join <- function(genus_rel_abund, genus_pool) {
  genus_join_lr_replicate <- inner_join(genus_rel_abund, genus_pool, by="taxon") %>%
    mutate(taxon = if_else(pool, "Other", taxon)) %>%
    group_by(location, sample_id, extraction_kit, taxon) %>%
    summarise(rel_abund = sum(rel_abund), mean = min(mean), .groups = "drop") %>%
    mutate(taxon = factor(taxon),
           taxon = fct_reorder(taxon, mean, .desc = FALSE)) %>%
    filter(rel_abund != 0)
  return(genus_join_lr_replicate)
}

###############################################################################################################
#relative abundance csv file
###############################################################################################################

rel_stat <- function(genus_join_lr_replicate) {
  relative_lr_csv <- genus_join_lr_replicate %>%
    group_by(location, taxon, extraction_kit) %>%
    summarise(mean = mean(rel_abund), sd = sd(rel_abund), median = median(rel_abund), .groups = "drop") %>%
    arrange(taxon)
  return(relative_lr_csv)
}

###############################################################################################################
#### perform comparison test
###############################################################################################################

rel_compare <- function(metadata, genus_join_lr_replicate) {
  print(by(genus_join_lr_replicate$rel_abund, genus_join_lr_replicate$extraction_kit, shapiro.test))
  
  genus_join_lr_replicate_comparison <- genus_join_lr_replicate %>%
    group_by(location, taxon) %>%
    t_test(rel_abund ~ extraction_kit, paired = F) %>%
    adjust_pvalue() %>%
    add_significance
  
  genus_join_lr_replicate_comparison <- genus_join_lr_replicate_comparison %>%
    add_xy_position(x = "taxon", dodge = 0.8)
  return(genus_join_lr_replicate_comparison)
}

###############################################################################################################
#plot the relative abundance
###############################################################################################################

taxon_plot <- function(genus_join_lr_replicate, genus_join_lr_replicate_comparison) {
  taxon_graph <- genus_join_lr_replicate %>%
    ggplot(aes(x=fct_reorder(taxon, mean, .desc = F), y=rel_abund)) +
    geom_boxplot(aes(color=extraction_kit, x=taxon), position = "dodge", lwd=1, width=0.4) +
    stat_pvalue_manual(genus_join_lr_replicate_comparison, hide.ns = TRUE, coord.flip = TRUE, 
                       label = "{p.adj.signif}", size=9, bracket.size=1) +
    scale_y_continuous(limits = c(NA, NA),
                       breaks = c(0.1, 1, 10, 100),
                       labels = c(0.1, 1, 10, 100)) +
    coord_trans_flip(y = "log10") +
    scale_color_manual("Extraction Kit", breaks=c("Dneasy", "PowerFecal", "Puregene"),
                       values = c("#549AAB", "#a73c5a", "#ff7954")) +
    scale_shape_manual("Location",
                       breaks=c("Australia", "Spain"),
                       values = c(0, 2)) +
    labs(x = "Genus", y="Relative abundance (%)") +
    theme_bw() +
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          strip.text.y = element_text(angle = 0),
          strip.background = element_blank()) +
    theme(text = element_text(size = 28))
  return(taxon_graph)
}

###############################################################################################################
### using the above functions to do analysis
###############################################################################################################
### bacteria analysis in dna kits
k_name <- "Bacteria"
bacteria_mr <- kingdom_subset(mr.rarefied_dna_kit, k_name)
bacteria_otu_relative <- otu_relative(bacteria_mr, metadata)
bacteria_relative_calculation <- relative_calculation(bacteria_otu_relative)

bacteria_taxon_pool_check <- taxon_pool_check(bacteria_relative_calculation)
bacteria_taxon_pool <- taxon_pool(bacteria_relative_calculation, 1.2)
bacteria_taxon_pool_join <- taxon_pool_join(bacteria_taxon_pool, bacteria_relative_calculation)

bacteria_rel_stat <- rel_stat(bacteria_taxon_pool_join)
bacteria_rel_compare <- rel_compare(metadata, bacteria_taxon_pool_join)

bacteria_taxon_dna_plot <- taxon_plot(bacteria_taxon_pool_join, bacteria_rel_compare)
bacteria_taxon_dna_plot

### archaea analysis in dna kits
k_name <- "Archaea"
archaea_mr <- kingdom_subset(mr.rarefied_dna_kit, k_name)
archaea_otu_relative <- otu_relative(archaea_mr, metadata)
archaea_relative_calculation <- relative_calculation(archaea_otu_relative)

archaea_taxon_pool_check <- taxon_pool_check(archaea_relative_calculation)
archaea_taxon_pool <- taxon_pool(archaea_relative_calculation, 2.1)
archaea_taxon_pool_join <- taxon_pool_join(archaea_taxon_pool, archaea_relative_calculation)

archaea_rel_stat <- rel_stat(archaea_taxon_pool_join)
archaea_rel_compare <- rel_compare(metadata, archaea_taxon_pool_join)

archaea_taxon_dna_plot <- taxon_plot(archaea_taxon_pool_join, archaea_rel_compare)
archaea_taxon_dna_plot

### eukaryota analysis in dna kits
k_name <- "Eukaryota"
eukaryota_mr <- kingdom_subset(mr.rarefied_dna_kit, k_name)
eukaryota_otu_relative <- otu_relative(eukaryota_mr, metadata)
eukaryota_relative_calculation <- relative_calculation(eukaryota_otu_relative)

eukaryota_taxon_pool_check <- taxon_pool_check(eukaryota_relative_calculation)
eukaryota_taxon_pool <- taxon_pool(eukaryota_relative_calculation, 2.27)
eukaryota_taxon_pool_join <- taxon_pool_join(eukaryota_taxon_pool, eukaryota_relative_calculation)

eukaryota_rel_stat <- rel_stat(eukaryota_taxon_pool_join)
eukaryota_rel_compare <- rel_compare(metadata, eukaryota_taxon_pool_join)

eukaryota_taxon_dna_plot <- taxon_plot(eukaryota_taxon_pool_join, eukaryota_rel_compare)
eukaryota_taxon_dna_plot