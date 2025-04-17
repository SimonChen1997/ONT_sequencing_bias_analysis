library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(patchwork)
library(ggfortify)
library(reshape2)
library(grid)
library(ggvenn)
library(gridExtra)
library(rstatix)

##########################################################################################################
### set workding directory
setwd("/home/folder")

##########################################################################################################
## function to clean and manipulate the data
read_mcaller_file <- function(path, lsk_file, rbk_file) {
  setwd(path)
  Spain_replicate_id <- c("lsk_code01", "lsk_code02", "lsk_code03", "lsk_code04", 
                          "rbk_code01", "rbk_code02", "rbk_code03", "rbk_code04")
  Australia_replicate_id <- c("lsk_16702", "lsk_16703", "lsk_16717", 
                              "rbk_16702", "rbk_16703", "rbk_16717")
  
  lsk_m6a <- read.csv(lsk_file, sep='\t', header=F)
  colnames(lsk_m6a) <- c("sample_id", "sequencing_kit", "contig", "seq", "position", "k-mer_context", "features", "dir", "label", "probability")
  lsk_m6a <- lsk_m6a %>%
    filter(label=="m6A") %>%
    mutate(identifier=paste0(contig,"_",position, "_", dir))
  
  lsk_m6a_distinct <- lsk_m6a %>%
    filter(label=="m6A") %>%
    mutate(identifier=paste0(contig,"_",position, "_", dir)) %>%
    distinct(identifier, .keep_all = T)
  
  rbk_m6a <- read.csv(rbk_file, sep='\t', header=F)
  colnames(rbk_m6a) <- c("sample_id", "sequencing_kit", "contig", "seq", "position", "k-mer_context", "features", "dir", "label", "probability")
  rbk_m6a <- rbk_m6a %>%
    filter(label=="m6A") %>%
    mutate(identifier=paste0(contig,"_",position, "_", dir))
  rbk_m6a_distinct <- rbk_m6a %>%
    filter(label=="m6A") %>%
    mutate(identifier=paste0(contig,"_",position, "_", dir)) %>%
    distinct(identifier, .keep_all = T)
  lsk_rbk_m6a <- rbind(lsk_m6a, rbk_m6a)
  lsk_rbk_m6a_distinct <- rbind(lsk_m6a_distinct, rbk_m6a_distinct)
  
  for (i in seq_len(nrow(lsk_rbk_m6a))) {
    if (lsk_rbk_m6a[i, "sample_id"] %in% Spain_replicate_id) {
      lsk_rbk_m6a[i,"location"] <- "Spain"
    }
    else {
      lsk_rbk_m6a[i,"location"] <- "Australia"
    }
  }
  
  for (i in seq_len(nrow(lsk_rbk_m6a_distinct))) {
    if (lsk_rbk_m6a_distinct[i, "sample_id"] %in% Spain_replicate_id) {
      lsk_rbk_m6a_distinct[i,"location"] <- "Spain"
    }
    else {
      lsk_rbk_m6a_distinct[i,"location"] <- "Australia"
    }
  }
  
  return(list(lsk = lsk_m6a, lsk_distinct=lsk_m6a_distinct, rbk = rbk_m6a, rbk_distinct=rbk_m6a_distinct, merged = lsk_rbk_m6a, merged_distinct = lsk_rbk_m6a_distinct))
}
##########################################################################################################
###### prevotella_ruminicola
prevotella_ruminicola_mcaller <- read_mcaller_file("/home/folder", 
                                                   "lsk_rumen_m6a.tsv", "rbk_rumen_m6a.tsv")
p.ruminicola_lsk <- prevotella_ruminicola_mcaller$lsk
p.ruminicola_rbk <- prevotella_ruminicola_mcaller$rbk
p.ruminicola_merged <- prevotella_ruminicola_mcaller$merged
p.ruminicola_distinct_merged <- prevotella_ruminicola_mcaller$merged_distinct
p.ruminicola_lsk_distinct <- prevotella_ruminicola_mcaller$lsk_distinct
p.ruminicola_rbk_distinct <- prevotella_ruminicola_mcaller$rbk_distinct

###### prevotella_bryantii
prevotella_bryantii_mcaller <- read_mcaller_file("/home/folder", 
                                                 "lsk_prevotella_bryantii.tsv", "rbk_prevotella_bryantii.tsv")
p.bryantii_lsk <- prevotella_bryantii_mcaller$lsk
p.bryantii_rbk <- prevotella_bryantii_mcaller$rbk
p.bryantii_merged <- prevotella_bryantii_mcaller$merged
p.bryantii_distinct_merged <- prevotella_bryantii_mcaller$merged_distinct
p.bryantii_lsk_distinct <- prevotella_bryantii_mcaller$lsk_distinct
p.bryantii_rbk_distinct <- prevotella_bryantii_mcaller$rbk_distinct

###### xanthomonas_hortorum
xanthomonas_hortorum_mcaller <- read_mcaller_file("/home/folder", 
                                                  "lsk_xanthomonas_hortorum.tsv", "rbk_xanthomonas_hortorum.tsv")
x.hortorum_lsk <- xanthomonas_hortorum_mcaller$lsk
x.hortorum_rbk <- xanthomonas_hortorum_mcaller$rbk
x.hortorum_merged <- xanthomonas_hortorum_mcaller$merged
x.hortorum_distinct_merged <- xanthomonas_hortorum_mcaller$merged_distinct
x.hortorum_lsk_distinct <- xanthomonas_hortorum_mcaller$lsk_distinct
x.hortorum_rbk_distinct <- xanthomonas_hortorum_mcaller$rbk_distinct

###### anabaena_sp_YBS01
anabaena_sp_YBS01_mcaller <- read_mcaller_file("/home/folder", 
                                               "lsk_anabaena_sp_YBS01.tsv", "rbk_anabaena_sp_YBS01.tsv")
a.sp_YBS01_lsk <- anabaena_sp_YBS01_mcaller$lsk
a.sp_YBS01_rbk <- anabaena_sp_YBS01_mcaller$rbk
a.sp_YBS01_merged <- anabaena_sp_YBS01_mcaller$merged
a.sp_YBS01_distinct_merged <- anabaena_sp_YBS01_mcaller$merged_distinct
a.sp_YBS01_lsk_distinct <- anabaena_sp_YBS01_mcaller$lsk_distinct
a.sp_YBS01_rbk_distinct <- anabaena_sp_YBS01_mcaller$rbk_distinct

###### fibrobacter_succinogenes
fibrobacter_succinogenes_mcaller <- read_mcaller_file("/home/folder", 
                                                      "lsk_fibrobacter_succinogenes.tsv", "rbk_fibrobacter_succinogenes.tsv")
f.succinogenes_lsk <- fibrobacter_succinogenes_mcaller$lsk
f.succinogenes_rbk <- fibrobacter_succinogenes_mcaller$rbk
f.succinogenes_merged <- fibrobacter_succinogenes_mcaller$merged
f.succinogenes_distinct_merged <- fibrobacter_succinogenes_mcaller$merged_distinct
f.succinogenes_lsk_distinct <- fibrobacter_succinogenes_mcaller$lsk_distinct
f.succinogenes_rbk_distinct <- fibrobacter_succinogenes_mcaller$rbk_distinct

lsk_rbk_m6a <- rbind(p.ruminicola_merged, p.bryantii_merged, x.hortorum_merged, a.sp_YBS01_merged, f.succinogenes_merged)
lsk_rbk_m6a_distinct <- rbind(p.ruminicola_distinct_merged, p.bryantii_distinct_merged, x.hortorum_distinct_merged, a.sp_YBS01_distinct_merged, f.succinogenes_distinct_merged)

##########################################################################################################
### a function for generateing venn graph on defined elements
setwd("/home/folder")
ggvenn_graphing <- function(lsk_object, rbk_object, species_name) {
  element_methylation_kits <- list(
    Ligation_kit = lsk_object$identifier,
    Rapid_kit = rbk_object$identifier
  )
  
  ggvenn(element_methylation_kits, fill_color = c("#EFC000FF", "#0073C2FF"),
         stroke_size = 0.5, set_name_size = 8, text_size = 8) +
    theme_void() +
    theme(text = element_text(face = "bold")) +
    theme(text = element_text(size = 16))
}

distinct_ggvenn <- ggvenn_graphing(lsk_rbk_m6a_distinct[lsk_rbk_m6a_distinct$sequencing_kit=="Ligation_kit",], 
                                   lsk_rbk_m6a_distinct[lsk_rbk_m6a_distinct$sequencing_kit== "Rapid_kit",], 
                                   "p.ruminicola - p.bryantii - x.hortorum - a.sp_YBS01 - f.succinogenes")
distinct_ggvenn

##########################################################################################################
### cahnge the direction to number
lsk_rbk_m6a$dir[lsk_rbk_m6a$dir == "+"] <- 0
lsk_rbk_m6a$dir[lsk_rbk_m6a$dir == "-"] <- 1
lsk_rbk_m6a$dir <- as.numeric(lsk_rbk_m6a$dir)
lsk_rbk_m6a$position <- as.numeric(lsk_rbk_m6a$position)

#### a function for generateing venn graph on defined elements
lsk_rbk_m6a_statistic <- function(merged_file){
  merged_statistic <- merged_file %>%
    group_by(location, sequencing_kit, sample_id) %>%
    summarise(methylated_num = n(), .groups = "drop")
  return(merged_statistic)
}

merged_distinct_statistic <- lsk_rbk_m6a_statistic(lsk_rbk_m6a_distinct)

merged_distinct_statistic <- merged_distinct_statistic %>%
  mutate(id=substring(sample_id, 5))
##########################################################################################################
### perform statistic analysis
by(merged_distinct_statistic$methylated_num, merged_distinct_statistic$sequencing_kit, shapiro.test)

merged_distinct_statistic_comparison_aus <- merged_distinct_statistic %>%
  filter(location=="Australia") %>%
  group_by(location) %>%
  t_test(data =., methylated_num ~ sequencing_kit, paired = F) %>%
  add_significance

merged_distinct_statistic_comparison_aus <- merged_distinct_statistic_comparison_aus %>%
  add_xy_position(dodge = 0.8)

merged_distinct_statistic_comparison_es <- merged_distinct_statistic %>%
  filter(location=="Spain") %>%
  group_by(location) %>%
  t_test(data =., methylated_num ~ sequencing_kit, paired = T) %>%
  add_significance

merged_distinct_statistic_comparison_es <- merged_distinct_statistic_comparison_es %>%
  add_xy_position(dodge = 0.8)

merged_distinct_statistic_summary <- merged_distinct_statistic %>%
  group_by(location, sequencing_kit) %>%
  summarise(mean_methylated_sites=mean(methylated_num), 
            sd_methylated_sites=sd(methylated_num), 
            median_methylated_sites=median(methylated_num), .groups = "drop")

##########################################################################################################
### plot data through oxplot
##########################################################################################################
### aus data

disctinct_boxplot_aus <- ggplot(merged_distinct_statistic[merged_distinct_statistic$location=="Australia", ], aes(x=sequencing_kit, y=methylated_num, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.5) +
  geom_point(shape=16, size=3) +
  geom_line(aes(group=id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location)) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(y="Numbers of Methylated Sites") +
  stat_pvalue_manual(merged_distinct_statistic_comparison_aus, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 26))

### spain data
disctinct_boxplot_es <- ggplot(merged_distinct_statistic[merged_distinct_statistic$location=="Spain", ], aes(x=sequencing_kit, y=methylated_num, color=sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.5) +
  geom_point(shape=16, size=3) +
  geom_line(aes(group=id), color="#3B3B3BFF", linetype="dotted", linewidth=1) +
  facet_grid(vars(location)) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  labs(y="Numbers of Methylated Sites") +
  stat_pvalue_manual(merged_distinct_statistic_comparison_es, hide.ns = TRUE, 
                     label = "{p.signif}", size=9, bracket.size=1) +
  stat_pvalue_manual(merged_distinct_statistic_comparison_es, hide.ns = TRUE,
                     y.position = 0.6, label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank()) +
  theme(text = element_text(size = 26))

##########################################################################################################
### fisher's exact test
##########################################################################################################
contingency_table_1 <- matrix(c(314,1005,697, 0), nrow = 2)
fisher.test(contingency_table_1, alternative = "less")


contingency_table_2 <- matrix(c(1848,3135,1912, 0), nrow = 2)
fisher.test(contingency_table_2, alternative = "less")

contingency_table_3 <- matrix(c(6424,2229,1881, 0), nrow = 2)
fisher.test(contingency_table_3, alternative = "less")

contingency_table_1 <- matrix(c(314,1005,697, 1996346), nrow = 2)
fisher.test(contingency_table_1, alternative = "greater")


contingency_table_2 <- matrix(c(1848,3135,1912, 2092643), nrow = 2)
fisher.test(contingency_table_2, alternative = "greater")

contingency_table_3 <- matrix(c(6424,2229,1881, 1777923), nrow = 2)
fisher.test(contingency_table_3, alternative = "greater")

