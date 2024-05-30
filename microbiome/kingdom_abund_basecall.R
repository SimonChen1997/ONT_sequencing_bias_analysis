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
library(grid)
library(cowplot)

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
merged_rumen_basecall = subset_samples(merged_rumen_seq_kit, 
                                       extraction_kit=="PowerFecal" | extraction_kit=="PowerSoil")

mr.rarefied_basecall = rarefy_even_depth(merged_rumen_basecall, rngseed=7127, sample.size=min(sample_sums(merged_rumen_seq_kit)), replace=F)

###############################################################################################################

set.seed(7127)

###############################################################################################################
#get rid of some empty taxonomy levels

Kingdom_level <- mr.rarefied_basecall

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
  mutate(Basecall_Mode = factor(basecall_mode, levels=c("FAST", "HAC", "SUP")))

###############################################################################################################
#calculate the relative abundance
###############################################################################################################

genus_rel_abund_lr <-  otu_rel_abund %>%
  filter(level=='Kingdom') %>%
  group_by(location, sequencing_kit, basecall_mode, sample_id, taxon) %>%
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  ungroup()

###############################################################################################################
#filter some low abundance taxons
###############################################################################################################

genus_pool_lr <- genus_rel_abund_lr %>%
  group_by(location, sequencing_kit, basecall_mode, taxon) %>%
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
  group_by(location, sequencing_kit, sample_id, basecall_mode, taxon) %>%
  summarise(rel_abund = sum(rel_abund), mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = FALSE)) %>%
  filter(rel_abund != 0)

###############################################################################################################
#relative abundance csv file
###############################################################################################################
relative_lr_csv <- genus_join_lr_replicate %>%
  group_by(location, sequencing_kit, basecall_mode, taxon) %>%
  summarise(mean = mean(rel_abund), sd = sd(rel_abund), median = median(rel_abund), .groups = "drop") %>%
  arrange(taxon)

write.csv(relative_lr_csv, row.names = F, 
          "/Users/JimmieBoice/Desktop/2023/Research/Sequencing_kit_bias/rumen_dna_seq_location/comparison_csv/relative_kingdom_basecall.csv")

###############################################################################################################
#perform comparison test
###############################################################################################################

genus_join_lr_replicate_comparison_aus <- genus_join_lr_replicate %>%
  filter(location=="Australia") %>%
  group_by(location, sequencing_kit, taxon) %>%
  t_test(rel_abund ~ basecall_mode, paired = F) %>%
  add_significance %>%
  mutate(location = "Australia") %>%
  select(location, everything())

genus_join_lr_replicate_comparison_aus <- genus_join_lr_replicate_comparison_aus %>%
  add_xy_position(x = "taxon", dodge = 0.8)

genus_join_lr_replicate_comparison_es <- genus_join_lr_replicate %>%
  filter(location=="Spain") %>%
  group_by(location, sequencing_kit, taxon) %>%
  t_test(rel_abund ~ basecall_mode, paired = T) %>%
  add_significance %>%
  mutate(location = "Spain") %>%
  select(location, everything())

genus_join_lr_replicate_comparison_es <- genus_join_lr_replicate_comparison_es %>%
  add_xy_position(x = "taxon", dodge = 0.8)

###############################################################################################################
#plot the relative abundance
###############################################################################################################

kingdom_basecall_aus <- genus_join_lr_replicate[genus_join_lr_replicate$location=="Australia",] %>%
  ggplot(aes(x=fct_reorder(taxon, mean, .desc = F), y=rel_abund)) +
  geom_boxplot(aes(color=basecall_mode, x=taxon), position = "dodge", lwd=1, width=0.4) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  stat_pvalue_manual(genus_join_lr_replicate_comparison_aus, hide.ns = TRUE, coord.flip = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(0.1, 1, 10, 100),
                     labels = c(0.1, 1, 10, 100)) +
  coord_trans_flip(y = "log10") +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  scale_shape_manual("Location",
                     breaks=c("Australia", "Spain"),
                     values = c(0, 2)) +
  labs(x = "Kingdom", y="Relative abundance (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))

kingdom_basecall_es <-genus_join_lr_replicate[genus_join_lr_replicate$location=="Spain",] %>%
  ggplot(aes(x=fct_reorder(taxon, mean, .desc = F), y=rel_abund)) +
  geom_boxplot(aes(color=basecall_mode, x=taxon), position = "dodge", lwd=1, width=0.4) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  stat_pvalue_manual(genus_join_lr_replicate_comparison_es, hide.ns = TRUE, coord.flip = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(0.1, 1, 10, 100),
                     labels = c(0.1, 1, 10, 100)) +
  coord_trans_flip(y = "log10") +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  scale_shape_manual("Location",
                     breaks=c("Australia", "Spain"),
                     values = c(0, 2)) +
  labs(x = "Kingdom", y="Relative abundance (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position="bottom",
        legend.title = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(text = element_text(size = 24))
kingdom_basecall_es
y.grob <- textGrob("Kingdom", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)

legend_basecall <- get_legend(kingdom_basecall_es + theme(legend.position="bottom",
                                                          legend.title = element_text(size = 24)))
kingdom_basecall <- grid.arrange(arrangeGrob(plot_grid(kingdom_basecall_aus, kingdom_basecall_es, ncol=1, align = 'v', rel_heights = c(1, 1.2)), left = y.grob))

