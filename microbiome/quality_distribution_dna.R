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

metadata_quality_dna <- fread("rumen_fastq_dna_quality_distribution.tsv", sep = "\t", header = T)

ggplot(metadata_quality_dna, aes(x=quality, color=basecall_mode, fill=basecall_mode)) +
  geom_density(alpha = .2) +
  facet_grid(.~extraction_kit) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  scale_fill_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                    values = c("#444574", "#7da66a", "#328cb0")) +
  theme_bw() +
  labs(y="Density", x="Read quality") +
  scale_x_continuous(limits = c(5, 55),
                     breaks = c(5, 15, 25, 35, 45, 55),
                     labels = c(5, 15, 25, 35, 45, 55)) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "bottom",
        plot.tag = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))
