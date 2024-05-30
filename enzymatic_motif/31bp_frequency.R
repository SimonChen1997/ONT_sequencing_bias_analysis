library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(patchwork)
library(gridExtra)

##################################################################################
# Background Frequency
#A	0.290006
#C	0.209623
#G	0.209716
#T	0.290644
#
##################################################################################
### set workding directory
setwd("/home/folder")

##################################################################################
### read files
lsk_rbk_31bp <- fread("lsk_rbk_tsu_31bp_base_frequency.csv", sep = ",", header = T)

Background_A <- data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-15:15)) {Background_A <- rbind(Background_A, BA=c("Background", "A", i, 29.0006, 0))}
colnames(Background_A) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

Background_C = data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-15:15)) {Background_C <- rbind(Background_C, BA=c("Background", "C", i, 20.9623, 0))}
colnames(Background_C) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

Background_G = data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-15:15)) {Background_G <- rbind(Background_G, BA=c("Background", "G", i, 20.9716, 0))}
colnames(Background_G) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

Background_T = data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(-15:15)) {Background_T <- rbind(Background_T, BA=c("Background", "T", i, 29.0644, 0))}
colnames(Background_T) <- c("sequencing_kit", "bases", "positions", "mean_percentages", "sd_percentages")

##################################################################################
### summary file

lsk_rbk_31bp_summary <- lsk_rbk_31bp %>%
  filter(lsk_rbk_31bp$bases != "N") %>%
  group_by(sequencing_kit, bases, positions) %>%
  summarise(mean_percentages=mean(percentages),
            sd_percentages=sd(percentages), .groups = "drop")

### minus 16 to allow the interaction site show "0" in the data
lsk_rbk_31bp_summary$positions <- lsk_rbk_31bp_summary$positions - 16

lsk_rbk_31bp_summary <- lsk_rbk_31bp_summary %>%
  rbind(Background_A, Background_C, Background_G, Background_T)

lsk_rbk_31bp_summary$positions <- as.integer(lsk_rbk_31bp_summary$positions)
lsk_rbk_31bp_summary$mean_percentages <- as.numeric(lsk_rbk_31bp_summary$mean_percentages)

lsk_rbk_31bp_overall_summary <- lsk_rbk_31bp %>%
  filter(lsk_rbk_31bp$bases != "N") %>%
  group_by(sequencing_kit, bases) %>%
  summarise(mean_percentages=mean(percentages),
            sd_percentages=sd(percentages), .groups = "drop")

##################################################################################
### plot data
ggplot(lsk_rbk_31bp_summary, aes(x=positions , y=mean_percentages, col=sequencing_kit, group=sequencing_kit)) + 
  geom_line(linewidth=2) +
  geom_point(size=3) + 
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
  scale_x_continuous(limits = c(-15, 15),
                     breaks = c(-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                     labels = c(-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) +
  theme(text = element_text(size = 28))
