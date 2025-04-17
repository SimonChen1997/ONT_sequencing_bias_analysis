library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(gridExtra)
library(lme4)

###############################################################################################################
setwd("/home/folder")

metadata <- read.csv("metadata.tsv", sep = "\t", header = T, row.names = 1)
read_n50 <- read.csv("rumen_read_stats.tsv", sep="\t", header = T)

###############################################################################################################
### function to handle dataset
id_addition <- function(input_matrix) {
  for (i in seq_len(nrow(input_matrix))) {
    input_matrix[i, "technical_replicate"] <- paste0(unlist(strsplit(as.character(input_matrix[i, "sample_id"]), "_"))[3])
  }
  return(input_matrix)
}

basecall_chr_to_num <- function(metafile){
  for (i in seq_len(nrow(metafile))) {
    if (metafile[i, "basecall_mode"] == "FAST") {
      metafile[i, "basecall_mode"] <- 1
    }
    else if (metafile[i, "basecall_mode"] == "HAC") {
      metafile[i, "basecall_mode"] <- 2
    }
    else {
      metafile[i, "basecall_mode"] <- 3
    }
  }
  for (i in seq_len(nrow(metafile))) {
    if (metafile[i, "animal_id"] == "SPA01") {
      metafile[i, "animal_id"] <- 1
    }
    else if (metafile[i, "animal_id"] == "SPA02") {
      metafile[i, "animal_id"] <- 2
    }
    else if (metafile[i, "animal_id"] == "SPA03") {
      metafile[i, "animal_id"] <- 3
    }
    else if (metafile[i, "animal_id"] == "SPA04") {
      metafile[i, "animal_id"] <- 4
    }
    else {
      metafile[i, "animal_id"] <- 5
    }
  }
  for (i in seq_len(nrow(metafile))) {
    if (metafile[i, "technical_replicate"] %in% c("16715", "16702", "16705")) {
      metafile[i, "technical_replicate"] <- "replicate_1"
    }
    else if (metafile[i, "technical_replicate"] %in% c("16718", "16703", "16723")) {
      metafile[i, "technical_replicate"] <- "replicate_2"
    }
    else if (metafile[i, "technical_replicate"] %in% c("16733", "16717", "16724")) {
      metafile[i, "technical_replicate"] <- "replicate_3"
    }
  }
  return(metafile)
}

###############################################################################################################
### subset the dataset
meta_n50 <- metadata %>%
  inner_join(read_n50, by = "sample_id", keep = F)

meta_n50 <- id_addition(meta_n50)
meta_n50 <- basecall_chr_to_num(meta_n50)

meta_n50_aus <- meta_n50[meta_n50$location=="Australia",]
meta_n50_es <- meta_n50[meta_n50$location=="Spain",]

meta_n50_aus[,c(2:3,7,10)] <- mutate_if(meta_n50_aus[, c(2:3,7,10)], is.character, factor)
meta_n50_es[,c(2:3,7,10)] <- mutate_if(meta_n50_es[, c(2:3,7,10)], is.character, factor)

###############################################################################################################
##### linear mix model on Australia dataset
### aus all variables
lm_classified_basecall_sequencing_extraction_technical_australia <- lmer(classified_read~basecall_mode+sequencing_kit+extraction_kit+(1|technical_replicate),
                                                                         data=meta_n50_aus)

Anova(lm_classified_basecall_sequencing_extraction_technical_australia)

### aus no basecall
lm_classified_sequencing_extraction_technical_australia <- lmer(classified_read~sequencing_kit+extraction_kit+(1|technical_replicate), 
                                                                data=meta_n50_aus)

anova(lm_classified_basecall_sequencing_extraction_technical_australia, 
      lm_classified_sequencing_extraction_technical_australia)

### aus sequencing
lm_classified_basecall_extraction_technical_australia <- lmer(classified_read~basecall_mode+extraction_kit+(1|technical_replicate), 
                                                              data=meta_n50_aus)

anova(lm_classified_basecall_sequencing_extraction_technical_australia, 
      lm_classified_basecall_extraction_technical_australia)

### aus extraction
lm_classified_basecall_sequencing_technical_australia <- lmer(classified_read~basecall_mode+sequencing_kit+(1|technical_replicate), 
                                                              data=meta_n50_aus)

anova(lm_classified_basecall_sequencing_extraction_technical_australia, 
      lm_classified_basecall_sequencing_technical_australia)

### aus no technical
lm_classified_basecall_mode_sequencing_extraction_australia <- lmer(classified_read~basecall_mode+sequencing_kit+extraction_kit, 
                                                                    data=meta_n50_aus)
anova(lm_classified_basecall_sequencing_extraction_technical_australia, 
      lm_classified_basecall_mode_sequencing_extraction_australia)

##### linear mix model on Spain dataset
### es all variables
lm_classified_basecall_sequencing_extraction_technical_spain <- lmer(classified_read~basecall_mode+sequencing_kit+(1|animal_id),
                                                                     data=meta_n50_es)

Anova(lm_classified_basecall_sequencing_extraction_technical_spain)

### es no basecall
lm_classified_sequencing_extraction_technical_spain <- lmer(classified_read~sequencing_kit+(1|animal_id), 
                                                            data=meta_n50_es)

anova(lm_classified_basecall_sequencing_extraction_technical_spain, 
      lm_classified_sequencing_extraction_technical_spain)

### es sequencing
lm_classified_basecall_extraction_technical_spain <- lmer(classified_read~basecall_mode+(1|animal_id), 
                                                          data=meta_n50_es)

anova(lm_classified_basecall_sequencing_extraction_technical_spain, 
      lm_classified_basecall_extraction_technical_spain)

