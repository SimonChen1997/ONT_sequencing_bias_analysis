library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(patchwork)
library(gridExtra)
library(ggfortify)

########################################################################################################################################
### set workding directory
setwd("/home/folder")
options(scipen=999) #no scientific notation

########################################################################################################################################
## read files
insertion_region_gc <- fread("region_chrom_insertion_gc.tsv", sep = "\t", header=T)
depth_region_gc <- fread("region_chrom_depth_gc.tsv", sep = "\t", header=T)
sample_length <- fread("sample_length_subsampled.tsv", sep = "\t", header=T)
primary_mapped <- fread("primary_mapped.tsv", sep = "\t", header=T)

########################################################################################################################################
## mergee file

merged_info <- insertion_region_gc %>%
  inner_join(depth_region_gc[,c(9,10)], by="id", keep = F)

merged_info <- merged_info %>%
  inner_join(sample_length[,c(2,4)], by="sample_id", keep=F) %>%
  inner_join( primary_mapped[,c(2,4)], by="sample_id", keep=F)

########################################################################################################################################
## normalized correlation
corr_depth_insertion <- ggplot(merged_info, aes(x=normalized_probability.x, 
                               normalized_probability.y, color=sequencing_kit)) +
  geom_point(size=3, alpha=0.7) +
  stat_cor(method = "pearson", p.accuracy=0.001, size=7) +
  geom_smooth(method = lm) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(x="Normalized interaction", y="Normalized coverage") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 28)) +
  guides(colour = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))