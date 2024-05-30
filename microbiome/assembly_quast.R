library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(rstatix)
library(gridExtra)
library(grid)
library(cowplot)

##################################################################################
### set workding directory
setwd("/home/folder")

##################################################################################
### read files
sample_quast <- fread("lsk_rbk_quast_info.tsv", sep = "\t", header = F, col.names = c("sample_id", "sequencing_kit", "contig_number", "N50", "location"))
read_n50 <- fread("/home/folder/rumen_dna_seq_location/fastq_stats/rumen_sub_read_stats.tsv", sep = "\t", header = T)

##################################################################################
### add the source of samples

merged_stat <- sample_quast %>%
  inner_join(read_n50, by="sample_id", keep = F)


id_addition <- function(input_matrix) {
  for (i in seq_len(nrow(input_matrix))) {
    input_matrix[i, "unique_id"] <- unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[3]
  }
  return(input_matrix)
}

sample_quast <- id_addition(sample_quast)

##################################################################################
#### contig stat and comparison
contig_number_summary <- sample_quast %>%
  group_by(location, sequencing_kit) %>%
  summarise(mean_contig_number = mean(contig_number),
            sd_contig_number = sd(contig_number))

contig_number_comparison_aus <- sample_quast %>%
  filter(location=="Australia") %>%
  group_by(location) %>%
  t_test(contig_number ~ sequencing_kit, paired = F) %>%
  add_significance %>%
  mutate(location = "Australia") %>%
  select(location, everything())

contig_number_comparison_aus <- contig_number_comparison_aus %>%
  add_xy_position(dodge = 0.8)

contig_number_comparison_es <- sample_quast %>%
  filter(location=="Spain") %>%
  group_by(location) %>%
  t_test(contig_number ~ sequencing_kit, paired = T) %>%
  add_significance %>%
  mutate(location = "Spain") %>%
  select(location, everything())

contig_number_comparison_es <- contig_number_comparison_es %>%
  add_xy_position(dodge = 0.8)

#########################
#### N50 stat and comparison
N50_summary <- sample_quast %>%
  group_by(location, sequencing_kit) %>%
  summarise(mean_N50 = mean(N50),
            sd_N50 = sd(N50), .groups = "drop")

N50_comparison_aus <- sample_quast %>%
  filter(location=="Australia") %>%
  group_by(location) %>%
  t_test(N50 ~ sequencing_kit, paired = F) %>%
  add_significance %>%
  mutate(location = "Australia") %>%
  select(location, everything())

N50_comparison_aus <- N50_comparison_aus %>%
  add_xy_position(dodge = 0.8)

N50_comparison_es <- sample_quast %>%
  filter(location=="Spain") %>%
  group_by(location) %>%
  t_test(N50 ~ sequencing_kit, paired = T) %>%
  add_significance %>%
  mutate(location = "Spain") %>%
  select(location, everything())

N50_comparison_es <- N50_comparison_es %>%
  add_xy_position(dodge = 0.8)

##################################################################################
### plot of contig data

contig_plot_aus <- ggplot(sample_quast[sample_quast$location=="Australia",], aes(x=sequencing_kit, y=contig_number, color=sequencing_kit)) +
  facet_grid(vars(location)) +
  geom_boxplot(lwd=1, width=0.5) +
  geom_point(shape=16, size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(y="Contig numbers") +
  stat_pvalue_manual(contig_number_comparison_aus, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

contig_plot_es <- ggplot(sample_quast[sample_quast$location=="Spain"], aes(x=sequencing_kit, y=contig_number, color=sequencing_kit)) +
  facet_grid(vars(location)) +
  geom_boxplot(lwd=1, width=0.5) +
  geom_point(shape=16, size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(y="Contig numbers") +
  stat_pvalue_manual(contig_number_comparison_es, hide.ns = TRUE,
                     label = "{p.signif}", size=9, bracket.size=1) +
  # stat_compare_means(
  #                    label.y = 5500, size=6, aes(label = paste0("p=", after_stat(p.format))),
  #                    comparisons = sample_length_comparisons,
  #                    inherit.aes= F) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

y.grob <- textGrob("Contig numbers", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)

contig_assembly_seq <- grid.arrange(arrangeGrob(plot_grid(contig_plot_aus, contig_plot_es, 
                                                          ncol=1, align = 'v', labels="A", label_size=24, hjust = 1, rel_heights = c(0.8,1)), left = y.grob))

##################################################################################
### plot of N50 data
N50_plot_aus <- ggplot(sample_quast[sample_quast$location=="Australia",], aes(x=sequencing_kit, y=N50, color=sequencing_kit)) +
  facet_grid(vars(location)) +
  geom_boxplot(lwd=1, width=0.5) +
  geom_point(shape=16, size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(y="Contig N50") +
  stat_pvalue_manual(N50_comparison_aus, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

N50_plot_es <- ggplot(sample_quast[sample_quast$location=="Spain"], aes(x=sequencing_kit, y=N50, color=sequencing_kit)) +
  facet_grid(vars(location)) +
  geom_boxplot(lwd=1, width=0.5) +
  geom_point(shape=16, size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(y="Contig N50") +
  stat_pvalue_manual(N50_comparison_es, hide.ns = TRUE,
                     label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

y.grob <- textGrob("Contig N50", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
contig_N50_seq <- grid.arrange(arrangeGrob(plot_grid(N50_plot_aus, N50_plot_es, 
                                                          ncol=1, align = 'v',  labels="B", label_size=24, hjust = 1, rel_heights = c(0.9,1)), left = y.grob))

grid.arrange(contig_assembly_seq, contig_N50_seq, ncol=2)

##################################################################################
### plot of N50 and contig data

contig_N50_plot <- ggplot(data = merged_stat, aes(x=N50_sub, y=contig_number, color=sequencing_kit)) +
  labs(x = "Read length N50", y = "Contig number") +
  geom_point(aes(color = sequencing_kit, shape=location), size=5) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  scale_shape_manual("Location", breaks=c("Australia", "Spain"),
                     values = c(16, 17)) +
  stat_cor(method = "pearson", p.accuracy=0.001, size=9) +
  geom_smooth(method=lm, se=FALSE, size =1) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 28)) +
  guides(colour = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))
contig_N50_plot

N50_N50_plot <- ggplot(data = merged_stat, aes(x=N50_sub, y=N50, color=sequencing_kit)) +
  labs(x = "Read length N50", y = "Contig N50") +
  geom_point(aes(color = sequencing_kit, shape=location), size=5) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  scale_shape_manual("Location", breaks=c("Australia", "Spain"),
                     values = c(16, 17)) +
  stat_cor(method = "pearson", p.accuracy=0.001, size=9) +
  geom_smooth(method=lm, se=FALSE, size =1) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 28)) +
  guides(colour = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))
N50_N50_plot

contig_N50_plot <- contig_N50_plot +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 28))

N50_N50_plot <- N50_N50_plot +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 28))

plot_grid(contig_N50_plot, N50_N50_plot, ncol=2, rel_widths = c(0.8,1))