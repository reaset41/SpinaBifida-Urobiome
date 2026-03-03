
##SB sequencing

library(tidyverse)
library(dada2)
library(here)
library(phyloseq)
library(decontam)
library(vegan)
library(ggpubr)

phyloseq_unfiltered

full_tax <- tax_table(phyloseq_unfiltered)
full_ASV <- otu_table(phyloseq_unfiltered)
TableAllASVs <- cbind(full_tax,t(full_ASV))

write.csv(TableAllASVs,file="TableFullASVs.csv")

##create subset of data excluding mock community samples
Patient_Blanks_Env <- subset_samples(phyloseq_unfiltered, !SampleType =="Mock")

##apply Decontam
contam_Patient <- isNotContaminant(Patient_Blanks_Env, neg = "is.neg", method ="prevalence",threshold = 0.5, normalize = TRUE, detailed = TRUE)

##how many ASVs are not contaminants 
sum(contam_Patient$not.contaminant)

##create physloseq object without contaminants
patient_true_blanks <- prune_taxa(contam_Patient$not.contaminant, Patient_Blanks_Env)

##remove chloroplast, mitochondrial, and unspecified taxa
patient_true_final <-subset_taxa(patient_true_blanks, Kingdom == "Bacteria"& Family != "Mitochondria" & Genus != "Mitochondria" & Order!="Chloroplast")

patient_true_final <- subset_samples(patient_true_final, !is.neg =="TRUE")

patient_rel <- 

is.na(otu_table(patient_true_final))

final_tax <- tax_table(patient_true_final)
final_ASV <- otu_table(patient_true_final)
FinalASVs <- cbind(final_tax,t(final_ASV))

write.csv(FinalASVs,file="FinalASVs.csv")