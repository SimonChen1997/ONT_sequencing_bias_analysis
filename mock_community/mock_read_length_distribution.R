library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(grid)
library(data.table)

###############################################################################################################
### set workding directory
setwd("/home/folder")

metadata_length <- fread("mock_seq_bias_quality_length_distribution.tsv", sep = "\t", header = T)
percentage_read_map <- fread("mock_seq_bias_stat_read_primary_bam.tsv", sep = "\t", header = T)
percentage_base_map <- fread("mock_seq_bias_stat_base_primary_bam.tsv", sep = "\t", header = T)
###############################################################################################################

ggplot(metadata_length, aes(x=length, color=sequencing_kit, fill=sequencing_kit)) +
  facet_grid(vars(genus), vars(mock_type)) +
  geom_density(alpha = .5) +
  scale_color_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  scale_fill_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                    values = c("#EFC000FF", "#0073C2FF")) +
  theme_bw() +
  labs(y="Density", x="Read length") +
  scale_x_log10() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.tag = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        panel.spacing = unit(2, "line")) +
  theme(text = element_text(size = 24))

ggsave("mock_seq_dna_length_distribusion.tiff", width = 12, height =8, dpi = 300)