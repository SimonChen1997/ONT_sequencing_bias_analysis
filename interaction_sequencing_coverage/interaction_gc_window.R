library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggridges)
library(patchwork)
library(gridExtra)

########################################################################################################################################
### set workding directory
setwd("/home/folder")
options(scipen=999) #no scientific notation

insertion_gc <- fread("lsk_rbk_insertion_gc_window.tsv", sep="\t", header = T)
insertion_gc <- insertion_gc %>%
  filter(insertion_number<100)

chromosome_gc <- fread("reference_10kb_GC.tsv", sep="\t", header = T)

########################################################################################################################################
## chromosome gc density
chromosome_gc$chromosome <- factor(chromosome_gc$chromosome, 
                                   levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                                              17,18,19,20,21,22,23,24,25,26,27,28,29,"X"))

ggplot(chromosome_gc, aes(x = gc, y=chromosome)) + 
  geom_density_ridges(fill="#868686FF", color="#4F4E52", alpha = .3, rel_min_height = 0.01) +
  labs(y="Chromosome", x="GC%") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 24))

ggplot(chromosome_gc, aes(x = gc)) + 
  geom_density(fill="#868686FF", color="#4F4E52", alpha = .3) +
  labs(y="Overall density", x="GC%") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 28))

########################################################################################################################################
### calculate the frequency of different gc contents
region_chrom_gc <- function(start_frame, end_frame, chromosome_gc_file) {
  region_frame=data.frame(start=start_frame,
                          end=end_frame)
  
  chrom_gc_column_name <- c("start", "end", "occurence")
  region_chrom_gc = data.frame(matrix(nrow=0, ncol=length(chrom_gc_column_name)))
  colnames(region_chrom_gc) = chrom_gc_column_name
  
  for (i in seq_len(nrow(region_frame))){
    frame_start <- region_frame[i,"start"]
    frame_end <- region_frame[i,"end"]
    region_data_chrom <- subset(chromosome_gc_file, 
                                chromosome_gc_file$gc >= frame_start &
                                  chromosome_gc_file$gc < frame_end)
    region_summary_chrom <- region_data_chrom %>%
      summarise(occurence=nrow(region_data_chrom), .groups = "drop")
    region_summary_chrom["start"] <- frame_start
    region_summary_chrom["end"] <- frame_end
    region_chrom_gc <- rbind(region_chrom_gc, region_summary_chrom)
  }
  
  region_chrom_gc <- region_chrom_gc %>%
    mutate(probability=100*(occurence/sum(occurence))+1) %>%
    ungroup() %>%
    mutate(region=paste0(start, "~", end, "%"))
  
  return(region_chrom_gc)
}

########################################################################################################################################
### calculate the frequency of insertion in different gc contents
start_frame=seq(0, 90, by=10)
end_frame=seq(10,100, by=10)

region_insertion_gc <- function(start_frame, end_frame, insertion_gc_file, chromosome_gc_file) {
  region_frame=data.frame(start=start_frame,
                          end=end_frame)
  
  column_name <- c("sequencing_kit", "sample_id", "total_insertion_number", "start", "end")
  region_insertion_gc = data.frame(matrix(nrow=0, ncol=length(column_name)))
  colnames(region_insertion_gc) = column_name
  
  for (i in seq_len(nrow(region_frame))){
    frame_start <- region_frame[i,"start"]
    frame_end <- region_frame[i,"end"]
    region_data <- subset(insertion_gc_file, 
                          insertion_gc_file$gc >= frame_start &
                            insertion_gc_file$gc < frame_end)
    region_summary <- region_data %>%
      group_by(sequencing_kit,sample_id) %>%
      summarise(total_insertion_number=sum(insertion_number), .groups = "drop")
    region_summary["start"] <- frame_start
    region_summary["end"] <- frame_end
    region_insertion_gc <- rbind(region_insertion_gc, region_summary)
  }
  
  region_insertion_gc <- region_insertion_gc %>%
    group_by(sequencing_kit, sample_id) %>%
    mutate(insertion_proportion=100*(total_insertion_number/sum(total_insertion_number))+1) %>%
    ungroup() %>%
    mutate(region=paste0(start, "~", end, "%"))
  
  ## calculate the normalized insertion frequency
  region_chrom_insertion_gc <- inner_join(region_insertion_gc, chromosome_gc_file[,c(4,5)], by="region", keep = F)
  
  region_chrom_insertion_gc <- region_chrom_insertion_gc %>%
    mutate(normalized_probability=log2(insertion_proportion/probability))
  return(region_chrom_insertion_gc)
}

########################################################################################################################################
## calculate the background probability and relative insertion number of in 0.5% interval gc windows
start_frame=seq(0, 99.5, by=0.5)
end_frame=seq(0.5,100, by=0.5)

region_chrom_gc_2 <- region_chrom_gc(start_frame, end_frame, chromosome_gc)
region_chrom_insertion_gc_2 <- region_insertion_gc(start_frame, end_frame, insertion_gc, region_chrom_gc_2)

region_chrom_insertion_gc_2 <- region_chrom_insertion_gc_2 %>%
  mutate(id = paste0(sample_id, "_",region))

region_chrom_insertion_gc_2_summary <- region_chrom_insertion_gc_2 %>%
  group_by(sequencing_kit, end, region) %>%
  summarise(mean_normalized_probability=mean(normalized_probability), 
            sd_normalized_probability=sd(normalized_probability), .groups = "drop")

########################################################################################################################################
## plot data for 0.5% interval
insertion_line_0.5 <- ggplot(region_chrom_insertion_gc_2_summary, aes(x=end, y=mean_normalized_probability, color=sequencing_kit)) +
  geom_line(linewidth=1, alpha=0.7) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  theme_bw() +
  scale_x_continuous(limits = c(25, 75),
                     breaks = seq(25, 75, by=5),
                     labels = seq(25, 75, by=5)) +
  labs(x="GC%", y="Normalized interaction") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 30))
insertion_line_0.5

insertion_point <- ggplot(region_chrom_insertion_gc_2, aes(x=end, y=normalized_probability, color=sequencing_kit)) +
  geom_point(size=3, alpha=0.6) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  theme_bw() +
  labs(x="GC%", y="Normalized interaction") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 30))

########################################################################################################################################
### calculate the insertion frequency of each gc windows (low and high gc regions)

start_frame=seq(0, 80, by=20)
end_frame=seq(20, 100, by=20)
region_chrom_gc_40 <- region_chrom_gc(start_frame, end_frame, chromosome_gc)
region_chrom_insertion_gc_40 <- region_insertion_gc(start_frame, end_frame, insertion_gc, region_chrom_gc_40)

region_chrom_insertion_gc_40 <- region_chrom_insertion_gc_40 %>%
  filter(region != "40~60%")

region_chrom_insertion_gc_40_summary <- region_chrom_insertion_gc_40 %>%
  group_by(sequencing_kit, region) %>%
  summarise(mean_depth=mean(normalized_probability), .groups = "drop")

########################################################################################################################################
## comparison files for low and high GC analysis

by(region_chrom_insertion_gc_40$normalized_probability, region_chrom_insertion_gc_40$sequencing_kit, shapiro.test)

region_chrom_insertion_gc_40_comparisons <- compare_means(normalized_probability ~ region, group.by = "sequencing_kit",
                                                      data=region_chrom_insertion_gc_40, method = "t.test", paired = FALSE)

region_chrom_insertion_gc_40_summary <- region_chrom_insertion_gc_40 %>%
  group_by(sequencing_kit, region) %>%
  summarise(mean_normalized_probability=mean(normalized_probability), 
            cd_normalized_probability=sd(normalized_probability), .groups = "drop")

########################################################################################################################################
## boxplot for low and high GC analysis
box_insertion_40 <- ggplot(region_chrom_insertion_gc_40, aes(x=region, y=normalized_probability, color=region)) +
  facet_grid(cols=vars(sequencing_kit)) +
  geom_boxplot(lwd=1, width=0.4) +
  scale_color_manual(breaks=c("20~40%", "60~80%"),
                     values = c("#CD534CFF", "#7AA6DCFF")) +
  labs(y="Normalized interaction", x="GC region") +
  stat_pvalue_manual(region_chrom_insertion_gc_40_comparisons, hide.ns = TRUE,
                     y.position = 0.6, label = "{p.signif}", size=9, bracket.size=1) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank()) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 30))

########################################################################################################################################
## boxplot for gc contents of the genome
overal_density <- ggplot(chromosome_gc, aes(x = gc)) + 
  geom_density(fill="#868686FF", color="#4F4E52", alpha = .3) +
  labs(y="Overall density", x="GC%") +
  theme_void() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 28))

insertion_line_0.5 <- ggplot(region_chrom_insertion_gc_2_summary, aes(x=end, y=mean_normalized_probability, color=sequencing_kit)) +
  geom_line(linewidth=1, alpha=0.7) +
  scale_color_manual("Sequencing Kit", breaks=c("Ligation_kit", "Rapid_kit"),
                     values = c("#EFC000FF", "#0073C2FF")) +
  theme_bw() +
  scale_x_continuous(limits = c(25, 75),
                     breaks = seq(25, 75, by=5),
                     labels = seq(25, 75, by=5)) +
  labs(x="GC%", y="Normalized interaction") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  theme(plot.tag = element_text(face = "bold")) +
  theme(text = element_text(size = 30))
