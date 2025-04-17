library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(patchwork)
library(gridExtra)
library(cowplot)
library(grid)
library(rstatix)

#####################################################################################################
### set workding directory
setwd("/home/folder")
options(scipen=999) #no scientific notation

mock_seq_n50_stat <-  read.csv("mock_seq_n50_stat.tsv", sep = "\t", header=T)

#####################################################################################################
mock_seq_n50_stat_summary <- mock_seq_n50_stat %>%
  group_by(mock_type, sequencing_kit, genus) %>%
  summarise(mean_n50=mean(n50),
            sd_n50=sd(n50), .groups = "drop")

mock_seq_n50_stat$genus <- factor(mock_seq_n50_stat$genus, 
                                  levels = c("Lactobacillus", "Escherichia", "Bifidobacterium"))

mock_seq_n50_compare_seq_kit <- mock_seq_n50_stat %>%
  group_by(mock_type, genus) %>%
  t_test(data = ., n50 ~ sequencing_kit, paired = F) %>%
  add_significance()

mock_seq_n50_compare_seq_kit <- mock_seq_n50_compare_seq_kit %>%
  add_xy_position(x=factor("genus"), dodge = 0.8)

mock_seq_n50_compare_genus <- mock_seq_n50_stat %>%
  group_by(mock_type) %>%
  t_test(data = ., n50 ~ genus, paired = F) %>%
  adjust_pvalue() %>%
  add_significance()

mock_seq_n50_compare_genus <- mock_seq_n50_compare_genus %>%
  add_xy_position(x=factor("genus"), step.increase=0.2)

#####################################################################################################
ggplot(mock_seq_n50_stat_summary, aes(x=factor(genus),
                                      y=mean_n50)) +
  facet_wrap(vars(mock_type)) +
  # geom_boxplot() +
  geom_col(aes(fill=sequencing_kit, color=sequencing_kit), position = position_dodge(width = 0.6), width = 0.6, alpha=0.6) +
  geom_errorbar(aes(ymin= mean_n50 - sd_n50,
                    ymax = mean_n50 + sd_n50,
                    fill=sequencing_kit, color=sequencing_kit),
                position = position_dodge(0.6), width = 0.3, colour = "black", alpha=0.6) +
  stat_pvalue_manual(mock_seq_n50_compare_seq_kit, hide.ns = TRUE,
                     label = "{p.signif}", size=9, bracket.size=1) +
  stat_pvalue_manual(mock_seq_n50_compare_genus, hide.ns = TRUE,
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  scale_color_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  scale_fill_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                    values = c("#EFC000FF", "#0073C2FF")) +
  # coord_trans(y="log1p") +
  # scale_y_continuous(limits = c(NA, NA),
  #                    breaks = c(0, 10, 100, 1000, 10000),
  #                    labels = c(0, 10, 100, 1000, 10000)) +
  labs(y="Read length N50") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(colour = "black"),
        strip.background = element_blank(),
        legend.position = "bottom") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 28))

ggsave("mock_seq_n50.tiff", width = 15, height =8, dpi = 300)