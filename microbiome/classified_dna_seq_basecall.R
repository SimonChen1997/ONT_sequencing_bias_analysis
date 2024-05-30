library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(grid)

###############################################################################################################
### set workding directory
setwd("/home/folder")

metadata <- read.csv("metadata.tsv", sep = "\t", header = T, row.names = 1)
read_n50 <- read.csv("rumen_read_stats.tsv", sep="\t", header = T)

###############################################################################################################

meta_n50 <- metadata %>%
  inner_join(read_n50, by = "sample_id", keep = F)

id_addition <- function(input_matrix) {
  for (i in seq_len(nrow(input_matrix))) {
    input_matrix[i, "animal_id"] <- paste0(unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[3])
    input_matrix[i, "unique_id"] <- paste0(unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[3], unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[4])
  }
  return(input_matrix)
}

meta_n50 <- id_addition(meta_n50)

meta_n50_seq_kit <- meta_n50[!(meta_n50$extraction_kit %in% "Puregene" | meta_n50$extraction_kit %in% "Dneasy") & meta_n50$basecall_mode %in% "SUP",]
meta_n50_dna_kit <- meta_n50[(meta_n50$location %in% "Australia" & !(meta_n50$sequencing %in% "Rapid_kit")) & meta_n50$basecall_mode %in% "SUP",]
meta_n50_basecall <- meta_n50[!(meta_n50$extraction_kit %in% "Puregene" | meta_n50$extraction_kit %in% "Dneasy"),]

###############################################################################################################

n50_classified_stat_seq_kit <- meta_n50_seq_kit %>%
  group_by(sequencing_kit, location) %>%
  summarise(mean_n50 = mean(N50_trim),
            mean_classified = mean(classified_read),
            sd_n50 = sd(N50_trim),
            sd_classified = sd(classified_read),
            median_n50=median(N50_trim),
            median_classified=median(classified_read))

n50_classified_stat_dna_kit <- meta_n50_dna_kit %>%
  group_by(extraction_kit, location) %>%
  summarise(mean_n50 = mean(N50_trim),
            mean_classified = mean(classified_read),
            sd_n50 = sd(N50_trim),
            sd_classified = sd(classified_read),
            median_n50=median(N50_trim),
            median_classified=median(classified_read))

n50_classified_stat_basecall_mode <- meta_n50_basecall %>%
  group_by(location, sequencing_kit, basecall_mode) %>%
  summarise(mean_n50 = mean(N50_trim),
            mean_classified = mean(classified_read),
            sd_n50 = sd(N50_trim),
            sd_classified = sd(classified_read),
            median_n50=median(N50_trim),
            median_classified=median(classified_read))

####################################################################################################
#### classified comparison seq

classified_compare_seq_kit_aus <- meta_n50_seq_kit %>%
  filter(location=="Australia") %>%
  t_test(data = ., classified_read ~ sequencing_kit, paired = F) %>%
  add_significance()%>%
  mutate(location = "Australia") %>%
  select(location, everything())

classified_compare_seq_kit_aus <- classified_compare_seq_kit_aus %>%
  add_xy_position(dodge = 0.8)


classified_compare_seq_kit_es <- meta_n50_seq_kit %>%
  filter(location=="Spain") %>%
  t_test(data = ., classified_read ~ sequencing_kit, paired = T) %>%
  add_significance() %>%
  mutate(location = "Spain") %>%
  select(location, everything())

classified_compare_seq_kit_es <- classified_compare_seq_kit_es %>%
  add_xy_position(dodge = 0.8)

#### classified comparison dna

by(meta_n50_dna_kit$classified_read, meta_n50_dna_kit$extraction_kit, shapiro.test)
classified_compare_dna_kit <- meta_n50_dna_kit %>%
  pairwise_t_test(data =., classified_read ~ extraction_kit, paired = F) %>%
  adjust_pvalue() %>%
  add_significance

classified_compare_dna_kit <- classified_compare_dna_kit %>%
  add_xy_position(dodge = 0.8)

classified_compare_basecall_mode <- meta_n50_basecall %>%
  group_by(location, sequencing_kit) %>%
  pairwise_t_test(data =., classified_read ~ basecall_mode, paired = T) %>%
  adjust_pvalue() %>%
  add_significance

#### classified comparison basecall

classified_compare_basecall_mode_aus <- meta_n50_basecall %>%
  filter(location=="Australia") %>%
  group_by(sequencing_kit) %>%
  t_test(data = ., classified_read ~ basecall_mode, paired = F) %>%
  add_significance()%>%
  mutate(location = "Australia") %>%
  select(location, everything())

classified_compare_basecall_mode_aus <- classified_compare_basecall_mode_aus %>%
  add_xy_position(dodge = 0.8)

classified_compare_basecall_mode_es <- meta_n50_basecall %>%
  filter(location=="Spain") %>%
  group_by(sequencing_kit) %>%
  t_test(data = ., classified_read ~ basecall_mode, paired = T) %>%
  add_significance() %>%
  mutate(location = "Spain") %>%
  select(location, everything())

classified_compare_basecall_mode_es <- classified_compare_basecall_mode_es %>%
  add_xy_position(dodge = 0.8)

####################################################################################################
#### N50 comparison seq

n50_compare_seq_kit_aus <- meta_n50_seq_kit %>%
  filter(location=="Australia") %>%
  t_test(data = ., N50_trim ~ sequencing_kit, paired = F) %>%
  add_significance()%>%
  mutate(location = "Australia") %>%
  select(location, everything())

n50_compare_seq_kit_aus <- n50_compare_seq_kit_aus %>%
  add_xy_position(dodge = 0.8)


n50_compare_seq_kit_es <- meta_n50_seq_kit %>%
  filter(location=="Spain") %>%
  t_test(data = ., N50_trim ~ sequencing_kit, paired = T) %>%
  add_significance() %>%
  mutate(location = "Spian") %>%
  select(location, everything())

n50_compare_seq_kit_es <- n50_compare_seq_kit_es %>%
  add_xy_position(dodge = 0.8)

by(meta_n50_dna_kit$N50_trim, meta_n50_dna_kit$extraction_kit, shapiro.test)
n50_compare_dna_kit <- meta_n50_dna_kit %>%
  group_by(location) %>%
  t_test(data =., N50_trim ~ extraction_kit, paired = F) %>%
  adjust_pvalue() %>%
  add_significance

n50_compare_dna_kit <- n50_compare_dna_kit %>%
  add_xy_position(dodge = 0.8)

#### N50 comparison basecall

n50_compare_basecall_mode_aus <- meta_n50_basecall %>%
  filter(location=="Australia") %>%
  group_by(sequencing_kit) %>%
  t_test(data = ., N50_trim ~ basecall_mode, paired = F) %>%
  add_significance()%>%
  mutate(location = "Australia") %>%
  select(location, everything())

n50_compare_basecall_mode_aus <- n50_compare_basecall_mode_aus %>%
  add_xy_position(dodge = 0.8)

n50_compare_basecall_mode_es <- meta_n50_basecall %>%
  filter(location=="Spain") %>%
  group_by(sequencing_kit) %>%
  t_test(data = ., N50_trim ~ basecall_mode, paired = T) %>%
  add_significance() %>%
  mutate(location = "Spain") %>%
  select(location, everything())

n50_compare_basecall_mode_es <- n50_compare_basecall_mode_es %>%
  add_xy_position(dodge = 0.8)


###############################################################################################################
#### N50 classified
lm_n50_classified <- ggplot(data = meta_n50[(meta_n50$location=="Australia"), ], aes(x=N50_trim, y=classified_read)) +
  labs(x = "Read length N50", y = "Percentage of classified reads") +
  geom_point(aes(color = extraction_kit, shape = sequencing_kit, size = basecall_mode)) +
  scale_colour_manual("Extraction Kit", breaks=c("Dneasy", "PowerFecal", "Puregene"),
                      values = c("#549AAB", "#a73c5a", "#ff7954")) +
  scale_size_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                    values = c(2, 4, 6)) +
  scale_shape_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c(16, 17)) +
  stat_cor(method = "pearson", p.accuracy=0.001, size=7) +
  geom_smooth(method=lm, se=FALSE, color="black", size =1) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(text = element_text(size = 24)) +
  theme(plot.tag = element_text(face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)),
         shape = guide_legend(override.aes = list(size=4)))
lm_n50_classified

###############################################################################################################
###### n50 seq aus
n50_seq_aus <- ggplot(meta_n50_seq_kit[meta_n50_seq_kit$location=="Australia",], aes(x=sequencing_kit, y=N50_trim, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  labs(y="Read length N50", x=NULL) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  stat_pvalue_manual(n50_compare_seq_kit_aus, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  labs(tag = "Australia") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        plot.tag.position = "right",
        plot.tag = element_text(size=22)) +
  theme(text = element_text(size = 24))
n50_seq_aus

###### n50 seq es
n50_seq_es <- ggplot(meta_n50_seq_kit[meta_n50_seq_kit$location=="Spain",], aes(x=sequencing_kit, y=N50_trim, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  labs(y="Read length N50", x=NULL) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  stat_pvalue_manual(n50_compare_seq_kit_es, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  labs(tag = "Spain") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        plot.tag.position = "right",
        plot.tag = element_text(size=22)) +
  theme(text = element_text(size = 24))
n50_seq_es

y.grob <- textGrob("Read length N50", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
n50_seq <- grid.arrange(arrangeGrob(plot_grid(n50_seq_aus, n50_seq_es, ncol=1, align = 'v', labels="A", label_size=28, hjust = 1), left = y.grob))

###############################################################################################################
###### n50 dna aus
n50_dna <- ggplot(meta_n50_dna_kit, aes(x=extraction_kit, y=N50_trim, color=extraction_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  labs(y="Read length N50", x=NULL) +
  scale_colour_manual("Extraction Kit", breaks=c("Dneasy", "PowerFecal", "Puregene"),
                      values = c("#549AAB", "#a73c5a", "#ff7954")) +
  stat_pvalue_manual(n50_compare_dna_kit, hide.ns = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))

n50_dna <- n50_dna +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

grid.arrange(n50_seq, n50_dna, ncol=2)

lm_n50_classified <- lm_n50_classified +
  labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

grid.arrange(grid.arrange(n50_seq, n50_dna, ncol=2), lm_n50_classified, nrow=2)

###############################################################################################################
###### n50 basecall aus
n50_basecall_aus <- ggplot(meta_n50_basecall[meta_n50_basecall$location=="Australia",], aes(x=basecall_mode, y=N50_trim, color=basecall_mode)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=animal_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  labs(y="Read length N50", x=NULL) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  stat_pvalue_manual(n50_compare_basecall_mode_aus, hide.ns = TRUE, 
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
n50_basecall_aus

###### n50 basecall aus
n50_basecall_es <- ggplot(meta_n50_basecall[meta_n50_basecall$location=="Spain",], aes(x=basecall_mode, y=N50_trim, color=basecall_mode)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=animal_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  labs(y="Read length N50", x=NULL) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  stat_pvalue_manual(n50_compare_basecall_mode_es, hide.ns = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(text = element_text(size = 24))
n50_basecall_es

y.grob <- textGrob("Read length N50", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
grid.arrange(arrangeGrob(plot_grid(n50_basecall_aus, n50_basecall_es, ncol=1, align = 'v'), left = y.grob))

###############################################################################################################
###### classified seq aus
classified_seq_aus <- ggplot(meta_n50_seq_kit[meta_n50_seq_kit$location=="Australia",], aes(x=sequencing_kit, y=classified_read, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  labs(y="Percentage of classified reads", x=NULL) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  stat_pvalue_manual(classified_compare_seq_kit_aus, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  labs(tag = "Australia") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        plot.tag.position = "right",
        plot.tag = element_text(size=22)) +
  theme(text = element_text(size = 22))
classified_seq_aus

###### classified seq es
classified_seq_es <- ggplot(meta_n50_seq_kit[meta_n50_seq_kit$location=="Spain",], aes(x=sequencing_kit, y=classified_read, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=unique_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  labs(y="Percentage of classified reads", x=NULL) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  stat_pvalue_manual(classified_compare_seq_kit_es, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  labs(tag = "Spain") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        plot.tag.position = "right",
        plot.tag = element_text(size=22)) +
  theme(text = element_text(size = 22))
classified_seq_es
y.grob <- textGrob("Percentage of classified reads", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)

grid.arrange(arrangeGrob(plot_grid(classified_seq_aus, classified_seq_es, ncol=1, align = 'v',  labels="A", label_size=24, hjust = 1), left = y.grob))
#################################################################################################
###### classified dna
classified_dna <- ggplot(meta_n50_dna_kit, aes(x=extraction_kit, y=classified_read, color=extraction_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  labs(y="Percentage of classified reads", x=NULL) +
  scale_colour_manual("Extraction Kit", breaks=c("Dneasy", "PowerFecal", "Puregene"),
                      values = c("#549AAB", "#a73c5a", "#ff7954")) +
  stat_pvalue_manual(classified_compare_dna_kit, hide.ns = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 24))

grid.arrange(arrangeGrob(plot_grid(classified_seq_aus, classified_seq_es, ncol=1, align = 'v',  labels="A", label_size=28), left = y.grob))
classified_seq <- grid.arrange(arrangeGrob(plot_grid(classified_seq_aus, classified_seq_es, ncol=1, align = 'v',  labels="A", label_size=28, hjust = 1), left = y.grob))
classified_dna <- classified_dna +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 24))

###### classified basecall aus
classified_basecall_aus <- ggplot(meta_n50_basecall[meta_n50_basecall$location=="Australia",], aes(x=basecall_mode, y=classified_read, color=basecall_mode)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=animal_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  labs(y="Percentage of classified reads", x=NULL) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  stat_pvalue_manual(classified_compare_basecall_mode_aus, hide.ns = TRUE, 
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
classified_basecall_aus

###### classified basecall es
classified_basecall_es <- ggplot(meta_n50_basecall[meta_n50_basecall$location=="Spain",], aes(x=basecall_mode, y=classified_read, color=basecall_mode)) +
  geom_boxplot(lwd=1, width=0.4) +
  geom_point(size=3) +
  geom_line(aes(group=animal_id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location), vars(sequencing_kit)) +
  labs(y="Percentage of classified reads", x=NULL) +
  scale_color_manual("Basecall Mode", breaks=c("FAST", "HAC", "SUP"),
                     values = c("#444574", "#7da66a", "#328cb0")) +
  stat_pvalue_manual(classified_compare_basecall_mode_es, hide.ns = TRUE, 
                     label = "{p.adj.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.x.top = element_blank()) +
  theme(text = element_text(size = 24))
classified_basecall_es

y.grob <- textGrob("Percentage of classified reads", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)

grid.arrange(arrangeGrob(plot_grid(classified_basecall_aus, classified_basecall_es, ncol=1, align = 'v'), left = y.grob))