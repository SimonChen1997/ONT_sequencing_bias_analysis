library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(patchwork)
library(gridExtra)
library(cowplot)
library(grid)

#####################################################################################################
### set workding directory
setwd("/home/folder")
options(scipen=999) #no scientific notation

percentage_base_map_normalize_all <-  read.csv("percentage_base_map_normalize_all.tsv", sep = "\t", header=T)
mock_seq_n50_stat <-  read.csv("mock_seq_n50_stat.tsv", sep = "\t", header=T)

#####################################################################################################
######
merged_normalized_percentage_n50 <- percentage_base_map_normalize_all[,c(1:4,9:10)] %>%
  inner_join(mock_seq_n50_stat[,c(1,6)], by="sample_id")

#####################################################################################################
######
ggplot(merged_normalized_percentage_n50, aes(x=n50,
                                        differential_normalzed_primary_base_map_percentage, color=sequencing_kit)) +
  # facet_grid(vars(genus)) +
  geom_point(size=3, alpha=0.7, aes(shape = genus)) +
  stat_cor(method = "pearson", p.accuracy=0.001, size=7) +
  geom_smooth(method = lm) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(x="Read length N50", y="Differential base percentage (%)", shape="Genus") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(colour = "black"),
        strip.background = element_blank(),
        legend.position = "right") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 28)) +
  guides(colour = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))

ggsave("mock_seq_n50_differential_relative.tiff", width = 12, height =8, dpi = 100)