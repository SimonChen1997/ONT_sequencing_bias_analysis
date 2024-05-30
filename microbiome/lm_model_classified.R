library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(gridExtra)

###############################################################################################################
### set workding directory
setwd("/home/folder")
metadata <- read.csv("metadata.tsv", sep = "\t", header = T, row.names = 1)
read_n50 <- read.csv("rumen_read_stats.tsv", sep="\t", header = T)

###############################################################################################################

meta_n50 <- metadata %>%
  inner_join(read_n50, by = "sample_id", keep = F)
meta_n50_seq_kit <- meta_n50[!(meta_n50$extraction_kit %in% "Puregene" | meta_n50$extraction_kit %in% "Dneasy") & meta_n50$basecall_mode %in% "SUP",]
meta_n50_dna_kit <- meta_n50[(meta_n50$location %in% "Australia" & !(meta_n50$sequencing %in% "Rapid_kit")) & meta_n50$basecall_mode %in% "SUP",]
meta_n50_basecall <- meta_n50[!(meta_n50$extraction_kit %in% "Puregene" | meta_n50$extraction_kit %in% "Dneasy"),]

###############################################################################################################

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
    if (metafile[i, "sequencing_kit"] == "Ligation_kit") {
      metafile[i, "sequencing_kit"] <- 1
    }
    else {
      metafile[i, "sequencing_kit"] <- 2
    }
  }
  for (i in seq_len(nrow(metafile))) {
    if (metafile[i, "extraction_kit"] == "Dneasy") {
      metafile[i, "extraction_kit"] <- 1
    }
    else if (metafile[i, "extraction_kit"] == "PowerFecal") {
      metafile[i, "extraction_kit"] <- 2
    }
    else if (metafile[i, "extraction_kit"] == "Puregene") {
      metafile[i, "extraction_kit"] <- 3
    }
    else {
      metafile[i, "extraction_kit"] <- 4
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
  return(metafile)
}

meta_n50_basecall <- basecall_chr_to_num(meta_n50_basecall)
meta_n50_dna_kit <- basecall_chr_to_num(meta_n50_dna_kit)
meta_n50_seq_kit <- basecall_chr_to_num(meta_n50_seq_kit)
meta_n50 <- basecall_chr_to_num(meta_n50)

meta_n50_basecall$basecall_mode <- as.numeric(meta_n50_basecall$basecall_mode)
meta_n50_dna_kit$basecall_mode <- as.numeric(meta_n50_dna_kit$basecall_mode)
meta_n50_seq_kit$basecall_mode <- as.numeric(meta_n50_seq_kit$basecall_mode)
meta_n50$basecall_mode <- as.numeric(meta_n50$basecall_mode)

#### change the other category variales to factor except basecalling
meta_n50[,c(2:3,9)] <- mutate_if(meta_n50[, c(2:3,9)], is.character, factor)
meta_n50_basecall[,c(2:3,9)] <- mutate_if(meta_n50_basecall[, c(2:3,9)], is.character, factor)
meta_n50_seq_kit[,c(2:3,9)] <- mutate_if(meta_n50_seq_kit[, c(2:3,9)], is.character, factor)
meta_n50_dna_kit[,c(2:3,9)] <- mutate_if(meta_n50_dna_kit[, c(2:3,9)], is.character, factor)

###########################################
### lm models
###########################################
### aus all variables
lm_classified_basecalling_sequencing_extraction_australia = lm(classified_read~basecall_mode+sequencing_kit+extraction_kit, data=meta_n50[meta_n50$location %in% "Australia",])
### aus no basecalling
lm_classified_sequencing_extraction_australia = lm(classified_read~sequencing_kit+extraction_kit, data=meta_n50[meta_n50$location %in% "Australia",])
### aus no sequencing
lm_classified_basecalling_extraction_australia = lm(classified_read~basecall_mode+extraction_kit, data=meta_n50[meta_n50$location %in% "Australia",])
### aus no extraction
lm_classified_basecalling_sequencing_australia = lm(classified_read~basecall_mode+sequencing_kit, data=meta_n50[meta_n50$location %in% "Australia",])

###########################################
### es all variables
lm_classified_basecalling_sequencing_animal_spain = lm(classified_read~basecall_mode+sequencing_kit+animal_id, data=meta_n50[meta_n50$location %in% "Spain",])
### es no basecalling
lm_classified_sequencing_animal_spain = lm(classified_read~sequencing_kit+animal_id, data=meta_n50[meta_n50$location %in% "Spain",])
### es no sequencing
lm_classified_basecalling_animal_spain = lm(classified_read~basecall_mode+animal_id, data=meta_n50[meta_n50$location %in% "Spain",])
### es no animal
lm_classified_basecalling_sequencing_spain = lm(classified_read~basecall_mode+sequencing_kit, data=meta_n50[meta_n50$location %in% "Spain",])

###########################################
### anova analysis
###########################################
### aus basecalling
anova(lm_classified_basecalling_sequencing_extraction_australia, lm_classified_sequencing_extraction_australia)

### aus sequencing
anova(lm_classified_basecalling_sequencing_extraction_australia, lm_classified_basecalling_extraction_australia)

### aus extraction
anova(lm_classified_basecalling_sequencing_extraction_australia, lm_classified_basecalling_sequencing_australia)

##################
### es basecalling
anova(lm_classified_basecalling_sequencing_animal_spain, lm_classified_sequencing_animal_spain)

### es sequencing
anova(lm_classified_basecalling_sequencing_animal_spain, lm_classified_basecalling_animal_spain)

### es animal
anova(lm_classified_basecalling_sequencing_animal_spain, lm_classified_basecalling_sequencing_spain)
