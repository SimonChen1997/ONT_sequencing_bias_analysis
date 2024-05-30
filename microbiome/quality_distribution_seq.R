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

metadata_quality <- fread("rumen_fastq_seq_quality_distribution.tsv", sep = "\t", header = T)

#qualtiy aus
qualtiy_aus <- ggplot(metadata_quality[metadata_quality$location=="Australia",], aes(x=quality, color=basecall_mode, fill=basecall_mode)) +
  geom_density(alpha = .2) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  scale_fill_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  theme_bw() +
  labs(y="Density", x="Read quality") +
  scale_x_continuous(limits = c(5, 40),
                     breaks = c(5, 10, 15, 20, 25, 30, 35, 40),
                     labels = c(5, 10, 15, 20, 25, 30, 35, 40)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))


#qualtiy es
qualtiy_es <- ggplot(metadata_quality[metadata_quality$location=="Spain",], aes(x=quality, color=basecall_mode, fill=basecall_mode)) +
  geom_density(alpha = .2) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  scale_fill_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                    values = c("#444574", "#7da66a", "#328cb0")) +
  theme_bw() +
  labs(y="Density", x="Read quality") +
  scale_x_continuous(limits = c(5, 40),
                     breaks = c(5, 10, 15, 20, 25, 30, 35, 40),
                     labels = c(5, 10, 15, 20, 25, 30, 35, 40)) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        strip.text.x.top = element_blank(),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))


y.grob <- textGrob("Density", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)

quality_legend <- get_legend(qualtiy_es + theme(legend.position="bottom",
                                                legend.title = element_text(size = 24, face = "bold")))

quality_density_plot <- grid.arrange(arrangeGrob(plot_grid(qualtiy_aus, 
                                                              qualtiy_es, ncol=1, align = 'v'), left=y.grob))
quality_density_plot_legend <- plot_grid(quality_density_plot, quality_legend, ncol = 1, nrow = 2, rel_heights = c(1, 0.1))
quality_density_plot_legend