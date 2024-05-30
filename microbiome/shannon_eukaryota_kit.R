library(phyloseq)
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(readr)
library(rstatix)
library(tidyverse)
library(ggpubr)
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
id_addition <- function(input_matrix) {
  for (i in seq_len(nrow(input_matrix))) {
    input_matrix[i, "animal_id"] <- paste0(unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[3])
    input_matrix[i, "unique_id"] <- paste0(unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[3], unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[4])
  }
  return(input_matrix)
}

metadata <- id_addition(metadata)
sampledata <- sample_data(metadata)

merged_rumen_kit = merge_phyloseq(merged_metagenomes, sampledata)

#subset the data
merged_rumen_kit <- subset_taxa(merged_rumen_kit, Kingdom == "Eukaryota")
merged_rumen_seq_kit = subset_samples(merged_rumen_kit, 
                                      (extraction_kit=="PowerFecal" | extraction_kit=="PowerSoil") & basecall_mode=="SUP")
merged_rumen_dna_kit = subset_samples(merged_rumen_kit, 
                                      location=="Australia" & sequencing_kit=="Ligation_kit" & basecall_mode=="SUP")
merged_rumen_basecall = subset_samples(merged_rumen_kit, 
                                       extraction_kit=="PowerFecal" | extraction_kit=="PowerSoil")

mr.rarefied_seq_kit = rarefy_even_depth(merged_rumen_seq_kit, rngseed=7127, sample.size=min(sample_sums(merged_rumen_seq_kit)), replace=F)
mr.rarefied_dna_kit = rarefy_even_depth(merged_rumen_dna_kit, rngseed=7127, sample.size=min(sample_sums(merged_rumen_dna_kit)), replace=F)
mr.rarefied_basecall = rarefy_even_depth(merged_rumen_basecall, rngseed=7127, sample.size=min(sample_sums(merged_rumen_kit)), replace=F)

###############################################################################################################

set.seed(7127)

###############################################################################################################

#calculation of index related to alpha diversity.
rich_overall_seq_kit = estimate_richness(mr.rarefied_seq_kit)
rich_overall_dna_kit = estimate_richness(mr.rarefied_dna_kit)
rich_overall_basecall = estimate_richness(mr.rarefied_basecall)

#combine the calculation results with rarefied meta data.
alpha_meta_seq_kit <- cbind(sample_data(mr.rarefied_seq_kit), rich_overall_seq_kit)
alpha_meta_dna_kit <- cbind(sample_data(mr.rarefied_dna_kit), rich_overall_dna_kit)
alpha_meta_basecall <- cbind(sample_data(mr.rarefied_basecall), rich_overall_basecall)

###############################################################################################################
### statistic analysis
#seq kit analysis
shannon_seq <- alpha_meta_seq_kit %>%
  select(sample_id, sequencing_kit, extraction_kit, Shannon, location) %>%
  group_by(location, extraction_kit, sequencing_kit) %>%
  summarise(mean = mean(Shannon), sd = sd(Shannon), median = median(Shannon), .groups = "drop") %>%
  arrange(sequencing_kit)

shannon_compare_seq_kit_aus <- alpha_meta_seq_kit %>%
  filter(location=="Australia") %>%
  group_by(location) %>%
  t_test(data =., Shannon ~ sequencing_kit, paired = F) %>%
  add_significance %>%
  mutate(location = "Australia") %>%
  select(location, everything())

shannon_compare_seq_kit_aus <- shannon_compare_seq_kit_aus %>%
  add_xy_position(dodge = 0.8)

shannon_compare_seq_kit_es <- alpha_meta_seq_kit %>%
  filter(location=="Spain") %>%
  group_by(location) %>%
  t_test(data =., Shannon ~ sequencing_kit, paired = T) %>%
  add_significance %>%
  mutate(location = "Spain") %>%
  select(location, everything())

shannon_compare_seq_kit_es <- shannon_compare_seq_kit_es %>%
  add_xy_position(dodge = 0.8)

#dna kit analysis
shannon_dna <- alpha_meta_dna_kit %>%
  select(sample_id, extraction_kit, sequencing_kit, Shannon, location) %>%
  group_by(location, sequencing_kit, extraction_kit) %>%
  summarise(mean = mean(Shannon), sd = sd(Shannon), median = median(Shannon), .groups = "drop") %>%
  arrange(extraction_kit)

shannon_compare_dna_kit <- alpha_meta_dna_kit %>%
  group_by(location) %>%
  t_test(data =., Shannon ~ extraction_kit, paired = F) %>%
  adjust_pvalue() %>%
  add_significance

shannon_compare_dna_kit <- shannon_compare_dna_kit %>%
  add_xy_position(dodge = 0.8)

#basecall analysis
shannon_basecall <- alpha_meta_basecall %>%
  select(sample_id, basecall_mode, sequencing_kit, Shannon, location) %>%
  group_by(location, sequencing_kit, basecall_mode) %>%
  summarise(mean = mean(Shannon), sd = sd(Shannon), median = median(Shannon), .groups = "drop") %>%
  arrange(basecall_mode)

shannon_compare_basecall_mode_aus <- alpha_meta_basecall %>%
  filter(location=="Australia") %>%
  group_by(location, sequencing_kit) %>%
  t_test(data =., Shannon ~ basecall_mode, paired = F) %>%
  add_significance %>%
  mutate(location = "Australia") %>%
  select(location, everything())

shannon_compare_basecall_mode_aus <- shannon_compare_basecall_mode_aus %>%
  add_xy_position(dodge = 0.8)

shannon_compare_basecall_mode_es <- alpha_meta_basecall %>%
  filter(location=="Spain") %>%
  group_by(location, sequencing_kit) %>%
  t_test(data =., Shannon ~ basecall_mode, paired = T) %>%
  add_significance %>%
  mutate(location = "Spain") %>%
  select(location, everything())

shannon_compare_basecall_mode_es <- shannon_compare_basecall_mode_es %>%
  add_xy_position(dodge = 0.8)

###############################################################################################################
#plot the Shannon graph
###############################################################################################################
##seq kit
shannon_seq_aus <- ggplot(alpha_meta_seq_kit[alpha_meta_seq_kit$location=="Australia",], aes(x=sequencing_kit, y=Shannon, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location)) +
  labs(y="Shannon index", x=NULL) +
  scale_shape_manual("Location",
                     breaks=c("Australia", "Spain"),
                     values = c(16, 17)) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  stat_pvalue_manual(shannon_compare_seq_kit_aus, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))

shannon_seq_es <- ggplot(alpha_meta_seq_kit[alpha_meta_seq_kit$location=="Spain",], aes(x=sequencing_kit, y=Shannon, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location)) +
  labs(y="Shannon index", x=NULL) +
  scale_shape_manual("Location",
                     breaks=c("Australia", "Spain"),
                     values = c(16, 17)) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  stat_pvalue_manual(shannon_compare_seq_kit_es, hide.ns = TRUE,
                     label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))

y.grob <- textGrob("Shannon index", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
shannon_seq <- grid.arrange(arrangeGrob(plot_grid(shannon_seq_aus, shannon_seq_es, ncol=1, align = 'v', 
                                                  labels="A", label_size=28, hjust = 1, rel_heights = c(0.9, 1)), left = y.grob))

## dna kit
shannon_dna <- ggplot(alpha_meta_dna_kit, aes(x=extraction_kit, y=Shannon, color=extraction_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  labs(y="Shannon index", x=NULL) +
  scale_colour_manual("Extraction Kit", breaks=c("Dneasy", "PowerFecal", "Puregene", "Unknown"),
                      values = c("#549AAB", "#a73c5a", "#ff7954")) +
  stat_pvalue_manual(shannon_compare_dna_kit, hide.ns = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))
shannon_dna

shannon_dna <- shannon_dna +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

grid.arrange(shannon_seq, shannon_dna, ncol=2)

## basecall mode
shannon_basecall_aus <- ggplot(alpha_meta_basecall[alpha_meta_basecall$location=="Australia",], aes(x=basecall_mode, y=Shannon, color=basecall_mode)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  geom_line(aes(group=animal_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  labs(y="Shannon index", x=NULL) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  stat_pvalue_manual(shannon_compare_basecall_mode_aus, hide.ns = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))

shannon_basecall_es <- ggplot(alpha_meta_basecall[alpha_meta_basecall$location=="Spain",], aes(x=basecall_mode, y=Shannon, color=basecall_mode)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  geom_line(aes(group=animal_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  labs(y="Shannon index", x=NULL) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  stat_pvalue_manual(shannon_compare_basecall_mode_es, hide.ns = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(text = element_text(size = 24))

y.grob <- textGrob("Shannon index", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
grid.arrange(arrangeGrob(plot_grid(shannon_basecall_aus, shannon_basecall_es, ncol=1, align = 'v'), left = y.grob))
