library(phyloseq)
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(deeptime)

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
merged_rumen_seq_kit = merge_phyloseq(merged_metagenomes, sampledata)

#subset the data to dna extraction
merged_rumen_dna_kit = subset_samples(merged_rumen_seq_kit, location=="Australia" & sequencing_kit=="Ligation_kit" & basecall_mode=="SUP")

mr.rarefied_dna_kit = rarefy_even_depth(merged_rumen_dna_kit, rngseed=7127, sample.size=min(sample_sums(merged_rumen_dna_kit)), replace=F)

###############################################################################################################

set.seed(7127)

###############################################################################################################
#get rid of some empty taxonomy levels

Kingdom_level <- mr.rarefied_dna_kit

###############################################################################################################

otu <- rownames(tax_table(Kingdom_level))
taxonomy <- cbind(otu, tax_table(Kingdom_level)) %>%
  as.data.frame()

otu_counts <- t(otu_table(Kingdom_level))
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

###############################################################################################################
#calculate the relative abundance
###############################################################################################################

genus_rel_abund_lr <-  otu_rel_abund %>%
  filter(level=='Kingdom') %>%
  group_by(location, extraction_kit, sample_id, taxon) %>%
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  ungroup()

###############################################################################################################
#filter some low abundance taxons
###############################################################################################################

genus_pool_lr <- genus_rel_abund_lr %>%
  group_by(location, extraction_kit, taxon) %>%
  summarise(mean = mean(rel_abund), .groups = "drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) < 0, 
            mean = mean(mean),
            .group="drop")

###############################################################################################################
#combine the data frame
###############################################################################################################

genus_join_lr_replicate <- inner_join(genus_rel_abund_lr, genus_pool_lr, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(location, sample_id, extraction_kit, taxon) %>%
  summarise(rel_abund = sum(rel_abund), mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = FALSE)) %>%
  filter(rel_abund != 0)

###############################################################################################################
#relative abundance csv file
###############################################################################################################
relative_lr_csv <- genus_join_lr_replicate %>%
  group_by(sample_id, taxon) %>%
  summarise(mean = mean(rel_abund), sd = sd(rel_abund), median = median(rel_abund), .groups = "drop") %>%
  arrange(taxon)

###############################################################################################################
#perform comparison test
###############################################################################################################
by(genus_join_lr_replicate$rel_abund, genus_join_lr_replicate$extraction_kit, shapiro.test)
genus_join_lr_replicate_comparison <- genus_join_lr_replicate %>%
  group_by(location, taxon) %>%
  t_test(rel_abund ~ extraction_kit, paired = F) %>%
  adjust_pvalue() %>%
  add_significance

genus_join_lr_replicate_comparison <- genus_join_lr_replicate_comparison %>%
  add_xy_position(x = "taxon", dodge = 0.8)

###############################################################################################################
#plot the relative abundance
###############################################################################################################

kingdom_dna <- genus_join_lr_replicate %>%
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
  labs(x = "Kingdom", y="Relative abundance (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 28))
