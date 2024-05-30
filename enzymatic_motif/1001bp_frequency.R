library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(patchwork)
library(gridExtra)

setwd("/Users/JimmieBoice/Desktop/2023/Research/Sequencing_kit_bias/results/TSU_1001bp_frequency")

#####################################
# Background Frequency
#A	0.290006
#C	0.209623
#G	0.209716
#T	0.290644
#
#####################################

##################################################################################
### read files
lsk_rbk_1001bp <- fread("lsk_rbk_tsu_1001bp_base_frequency.csv", sep = ",", header = T)

Background_A <- data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-500:500)) {Background_A <- rbind(Background_A, BA=c("Background", "A", i, 29.0006, 0))}
colnames(Background_A) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

Background_C = data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-500:500)) {Background_C <- rbind(Background_C, BA=c("Background", "C", i, 20.9623, 0))}
colnames(Background_C) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

Background_G = data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-500:500)) {Background_G <- rbind(Background_G, BA=c("Background", "G", i, 20.9716, 0))}
colnames(Background_G) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

Background_T = data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-500:500)) {Background_T <- rbind(Background_T, BA=c("Background", "T", i, 29.0644, 0))}
colnames(Background_T) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

##################################################################################
### summary file

lsk_rbk_1001bp_summary <- lsk_rbk_1001bp %>%
  filter(lsk_rbk_1001bp$bases != "N") %>%
  group_by(sequencing_kit, bases, positions) %>%
  summarise(mean_percentages=mean(percentages),
            sd_percentages=sd(percentages), .groups = "drop")

lsk_rbk_1001bp_summary$positions <- lsk_rbk_1001bp_summary$positions - 501

lsk_rbk_1001bp_summary <- lsk_rbk_1001bp_summary %>%
  rbind(Background_A, Background_C, Background_G, Background_T)

lsk_rbk_1001bp_summary$positions <- as.integer(lsk_rbk_1001bp_summary$positions)
lsk_rbk_1001bp_summary$mean_percentages <- as.numeric(lsk_rbk_1001bp_summary$mean_percentages)

lsk_rbk_1001bp_overall_summary <- lsk_rbk_1001bp %>%
  filter(lsk_rbk_1001bp$bases != "N") %>%
  group_by(sequencing_kit, bases) %>%
  summarise(mean_percentages=mean(percentages),
            sd_percentages=sd(percentages), .groups = "drop")

write_csv(lsk_rbk_1001bp_overall_summary,
          "/Users/JimmieBoice/Desktop/2023/Research/Sequencing_kit_bias/results/TSU_1001bp_frequency/1001bp_nucleotide_overal_comparison.csv")

##################################################################################
### comparison file
lsk_rbk_1001bp_a <- lsk_rbk_1001bp %>%
  filter(bases=="A")

lsk_rbk_1001bp_t <- lsk_rbk_1001bp %>%
  filter(bases=="T")

lsk_rbk_1001bp_a_comparison <- compare_means(percentages ~ sequencing_kit, 
                                           data=lsk_rbk_1001bp_a, method = "wilcox.test", paired = FALSE)

lsk_rbk_1001bp_t_comparison <- compare_means(percentages ~ sequencing_kit, 
                                           data=lsk_rbk_1001bp_t, method = "wilcox.test", paired = FALSE)

##################################################################################
### plot data
ggplot(lsk_rbk_1001bp_summary, aes(x=positions , y=mean_percentages, col=sequencing_kit, group=sequencing_kit)) + 
  geom_line(linewidth=0.8) +
  geom_point(size=0.8) + 
  facet_grid(bases~. ) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted")) +
  theme_bw() +
  scale_color_manual(breaks=c("Ligation_kit", "Rapid_kit", "Background"),
                     values = c("#EFC000FF", "#0073C2FF", "#868686FF")) +
  labs(x="Positions", y="Percentages") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(plot.tag = element_text(face = "bold")) +
  # scale_x_continuous(limits = c(-15, 15),
  #                    breaks = c(-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
  #                    labels = c(-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) +
  theme(text = element_text(size = 28))
ggsave("lsk_rbk_1001bp_frequency.tiff", width = 18, height = 9, dpi = 300)
