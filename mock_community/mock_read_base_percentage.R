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

percentage_read_map <- fread("mock_seq_bias_stat_read_primary_bam.tsv", sep = "\t", header = T)
percentage_base_map <- fread("mock_seq_bias_stat_base_primary_bam.tsv", sep = "\t", header = T)
mock_seq_n50_stat <-  read.csv("mock_seq_n50_stat.tsv", sep = "\t", header=T)

###############################################################################################################
#### normalization based on the mock sample DNA ratios
###############################################################################################################

### in mock 1, the ratio is Lactobacillus:Escherichia:Bifidobacterium = 1:1:1
### in mock 2, the ratio is Lactobacillus:Escherichia:Bifidobacterium = 1:2:3

mock_1_probability <- data.frame(genus=c("Lactobacillus", "Escherichia", "Bifidobacterium"),
                                 probability=c(100*1/3, 100*1/3, 100*1/3),
                                 mock_type=c("mock_1", "mock_1", "mock_1"))

mock_2_probability <- data.frame(genus=c("Lactobacillus", "Escherichia", "Bifidobacterium"),
                                 probability=c(100*1/6, 100*2/6, 100*3/6),
                                 mock_type=c("mock_2", "mock_2", "mock_2"))

mock_ratio_fold_change <- function(relative_abundance_form, probability_file, normalize_col_name){
  normalized_probability_form <- as.data.frame(matrix(nrow = 0, ncol=6))
  colnames(normalized_probability_form) <- c("sample_id", "technical_replicate", "mock_type", "sequencing_kit", "taxon", normalize_col_name)
  
  mutate_col_name=paste0("differential_", normalize_col_name)
  for (i in seq_len(nrow(probability_file))) {
    
    genus_type <- probability_file[i,"genus"]
    probability <- probability_file[i,"probability"]
    mock_catogery <- probability_file[i,"mock_type"]
    
    inter_form <- subset(relative_abundance_form, relative_abundance_form$genus == genus_type & relative_abundance_form$mock_type == mock_catogery) %>%
      mutate(!!mutate_col_name := (!!sym(normalize_col_name) - probability))
    
    
    normalized_probability_form <- rbind(normalized_probability_form, inter_form)
  }
  return(normalized_probability_form)
}

###############################################################################################################
#### normalization based on the mock sample DNA ratios
percentage_read_map_normalize <- percentage_read_map %>%
  group_by(sequencing_kit, mock_type, technical_replicate) %>%
  mutate(normalzed_map_percentage=100*map_percentage/sum(map_percentage),
         normalzed_primary_map_percentage=100*primary_map_percentage/sum(primary_map_percentage))

percentage_base_map_normalize <- percentage_base_map %>%
  group_by(sequencing_kit, mock_type, technical_replicate) %>%
  mutate(normalzed_primary_base_map_percentage=100*primary_base_map_percentage/sum(primary_base_map_percentage))

###############################################################################################################
#### normalization based on the mock sample DNA ratios
read_mock_1_test <- mock_ratio_fold_change(percentage_read_map_normalize, mock_1_probability, "normalzed_map_percentage")
read_mock_1_test <- mock_ratio_fold_change(read_mock_1_test, mock_1_probability, "normalzed_primary_map_percentage")

read_mock_2_test <- mock_ratio_fold_change(percentage_read_map_normalize, mock_2_probability, "normalzed_map_percentage")
read_mock_2_test <- mock_ratio_fold_change(read_mock_2_test, mock_2_probability, "normalzed_primary_map_percentage")

read_mock_seq_bias_normalized <- rbind(read_mock_1_test, read_mock_2_test)

base_mock_1_test <- mock_ratio_fold_change(percentage_base_map_normalize, mock_1_probability, "normalzed_primary_base_map_percentage")
base_mock_2_test <- mock_ratio_fold_change(percentage_base_map_normalize, mock_2_probability, "normalzed_primary_base_map_percentage")

base_mock_seq_bias_normalized <- rbind(base_mock_1_test, base_mock_2_test)

###############################################################################################################
#### normalization based on the mock sample DNA ratios
read_mock_seq_bias_normalized_summary <- read_mock_seq_bias_normalized %>%
  group_by(sequencing_kit, mock_type, genus) %>%
  summarise(mean_differential_normalzed_primary_map_percentage=mean(differential_normalzed_primary_map_percentage),
            sd_differential_normalzed_primary_map_percentage=sd(differential_normalzed_primary_map_percentage), .groups = "drop")

base_mock_seq_bias_normalized_summary <- base_mock_seq_bias_normalized %>%
  group_by(sequencing_kit, mock_type, genus) %>%
  summarise(mean_differential_normalzed_primary_base_map_percentage=mean(differential_normalzed_primary_base_map_percentage),
            mean_normalzed_primary_base_map_percentage=mean(normalzed_primary_base_map_percentage),
            sd_differential_normalzed_primary_base_map_percentage=sd(differential_normalzed_primary_base_map_percentage), 
            sd_normalzed_primary_base_map_percentage=sd(normalzed_primary_base_map_percentage), .groups = "drop")

base_mock_seq_bias_normalized_no_mock_summary <- base_mock_seq_bias_normalized %>%
  group_by(sequencing_kit, genus) %>%
  summarise(mean_differential_normalzed_primary_base_map_percentage=mean(differential_normalzed_primary_base_map_percentage),
            sd_differential_normalzed_primary_base_map_percentage=sd(differential_normalzed_primary_base_map_percentage), .groups = "drop")

base_mock_seq_bias_normalized_summary$genus <- factor(base_mock_seq_bias_normalized_summary$genus, 
                                                      levels = c("Lactobacillus", "Escherichia", "Bifidobacterium"))

mock_seq_n50_stat_summary <- mock_seq_n50_stat %>%
  group_by(mock_type, sequencing_kit, genus) %>%
  summarise(mean_n50=mean(n50),
            sd_n50=sd(n50), .groups = "drop")

base_mock_seq_bias_normalized_summary <- base_mock_seq_bias_normalized_summary %>%
  mutate(id=paste(sequencing_kit,mock_type,genus, sep = "_"))

mock_seq_n50_stat_summary <- mock_seq_n50_stat_summary %>%
  mutate(id=paste(sequencing_kit,mock_type,genus, sep = "_"))

base_mock_seq_bias_normalized_summary <- base_mock_seq_bias_normalized_summary %>%
  inner_join(mock_seq_n50_stat_summary[,c(4:6)], by="id", keep=F)

###############################################################################################################
base_mock_seq_bias_normalized_compare <- base_mock_seq_bias_normalized %>%
  group_by(mock_type, genus) %>%
  t_test(data = ., differential_normalzed_primary_base_map_percentage ~ sequencing_kit, paired = F) %>%
  add_significance()

base_mock_seq_bias_normalized_compare <- base_mock_seq_bias_normalized_compare %>%
  add_xy_position(dodge = 0.6)

base_mock_seq_bias_normalized_compare$y.position <- base_mock_seq_bias_normalized_compare$y.position + 9

###############################################################################################################
ggplot(base_mock_seq_bias_normalized_summary, aes(x=sequencing_kit, y=mean_differential_normalzed_primary_base_map_percentage)) +
  facet_grid(vars(mock_type), vars(as.factor(genus)), scales = "free_y") +
  geom_point(position = position_dodge(width = 0.4), size = 2, aes(colour = sequencing_kit, fill = sequencing_kit)) +
  geom_col(position = position_dodge(width = 0.4), width = 0.5, alpha=0.7, aes(colour = sequencing_kit, fill = sequencing_kit)) +
  geom_errorbar(aes(ymin= mean_differential_normalzed_primary_base_map_percentage - sd_differential_normalzed_primary_base_map_percentage, 
                    ymax = mean_differential_normalzed_primary_base_map_percentage + sd_differential_normalzed_primary_base_map_percentage), position = position_dodge(0.4), width = 0.3, colour = "black", alpha=0.8) +
  scale_color_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  scale_fill_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                    values = c("#EFC000FF", "#0073C2FF")) +
  stat_pvalue_manual(base_mock_seq_bias_normalized_compare, hide.ns = TRUE,
                     label = "{p.signif}", size=9, bracket.size=1) +
  labs(y="Differential base percentage (%)",
       tag = "mock1: Lactobacillus (34.5% GC) : Escherichia (50% GC) : Bifidobacterium (64% GC) = 1 : 1 : 1\nmock2: Lactobacillus (34.5% GC) : Escherichia (50% GC) : Bifidobacterium (64% GC) = 1 : 2 : 3") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 0.9),
        legend.title = element_text(face = "bold"),
        legend.position = "none",
        plot.tag.position = "top",
        plot.tag = element_text(face = "bold", size = 18),
        strip.text = element_text(colour = "black"),
        strip.background = element_blank()) +
  theme(text = element_text(size = 22))

ggsave("mock_seq_dna_primary_mapped_base_proportion.tiff", width = 13, height =9, dpi = 300)

ggplot(read_mock_seq_bias_normalized_summary, aes(x=factor(genus, levels = c("Lactobacillus", "Escherichia", "Bifidobacterium")), y=mean_differential_normalzed_primary_map_percentage, 
                                                  colour = sequencing_kit, fill = sequencing_kit)) +
  facet_wrap(vars(mock_type), scales = "free_x") +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6, alpha=0.6) +
  geom_errorbar(aes(ymin= mean_differential_normalzed_primary_map_percentage - sd_differential_normalzed_primary_map_percentage, 
                    ymax = mean_differential_normalzed_primary_map_percentage + sd_differential_normalzed_primary_map_percentage), 
                position = position_dodge(0.6), width = 0.3, colour = "black", alpha=0.8) +
  scale_color_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  scale_fill_manual("Sequencing kit", breaks=c("Ligation_kit", "Rapid_kit"),
                    values = c("#EFC000FF", "#0073C2FF")) +
  labs(y="Differential read percentage (%)",
       tag = "mock1: Lactobacillus : Escherichia : Bifidobacterium = 1 : 1 : 1\nmock2: Lactobacillus : Escherichia : Bifidobacterium = 1 : 2 : 3") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 0.9),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        plot.tag.position = "top",
        plot.tag = element_text(face = "bold", size = 16),
        strip.text = element_text(colour = "black"),
        strip.background = element_blank()) +
  theme(text = element_text(size = 22))

ggsave("mock_seq_dna_primary_mapped_read_proportion.tiff", width = 9, height =8.5, dpi = 300)