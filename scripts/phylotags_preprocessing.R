library(ggplot2)
library(dplyr)
library(taxize)
library(stringr)
library(taxize)
library(tidyr)
library(RColorBrewer)
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(microbiome)
library(knitr)
library(ape)
library(ggpubr)
library(cowplot)
library(data.table)
library(broom)
library(purrr)
library(gridExtra)
#library(ggtree)
library(DESeq2)

# replace "marker" for the actual marker: "ssu", "dnak","gyrB"

## METADATA 

#Metadata files with  info for each library
marker_meta = read.table("metadata_libs/marker.files", header = F)

#edit and remove library names
marker_meta$lib <- sapply(strsplit(as.character(marker_meta[,2]), "\\_"), `[`, 1)

#metadata bmi
sample_info <- read.table("bsp_classification.txt",header=T)
rownames(sample_info) <- sapply(strsplit(as.character(sample_info$ID),"_H"), `[`, 1)

#metadata for design variables
microbio_selected <- read.table("microbio_selected.meta.txt",header=T)

#remove the seconde records (records with "H2") -> el id "MI_093_H12" está mal anotado fue manualmente corregido
microbio_selected <- microbio_selected[which(sapply(strsplit(as.character(microbio_selected$ID),"_"), `[`, 3)!="H2"),]
rownames(microbio_selected) <- sapply(strsplit(as.character(microbio_selected$ID),"_H"), `[`, 1)
microbio_selected$bmi_class <- factor(sapply(strsplit(as.character(microbio_selected$bmi_class),"-"), `[`, 1),levels = c("Lean","Overweight","Obese"))

#Split body-size-phenotypes to get the BMI index
sample_info$bmi <- factor(sapply(strsplit(as.character(sample_info$bsp_class),"-"), `[`, 1),levels = c("Normoweight","Overweight","Obese"))

microbio_selected$chs_class <- sample_info[rownames(sample_info) %in% rownames(microbio_selected),"chs_class"]
#filter to keep only the samples used in phyloTags

## FEATURES-COUNTS 

#Feature files: counts of library for each sample
#these files start with # wich R does not read because interprets it as comment
#Read the header extract it and then name the columns from the file
marker_feature_dada2 <- read.table("qiime_dada_deblur_complete/marker_Tabla_Dada2.txt", header = F)
marker_con <- file("qiime_dada_deblur_complete/marker_Tabla_Dada2.txt", "r")
marker_header <- readLines(marker_con, n=2)
colnames(marker_feature_dada2) <- strsplit(marker_header, "\\s+")[[2]][-1] 
close(marker_con)

#asv table of dada2
marker_mat_dada2 = as.matrix(marker_feature_dada2[,-1])
rownames(marker_mat_dada2) <- marker_feature_dada2[,1]

#OTU tables dada2
marker_OTU_dada2 = otu_table(marker_mat_dada2, taxa_are_rows = TRUE)
colnames(marker_OTU_dada2) <- marker_meta[match(colnames(marker_OTU_dada2),marker_meta$lib),]$V1
colnames(marker_OTU_dada2) <- gsub("-", "_", colnames(marker_OTU_dada2)) #change lib names for sample names (low stripes)
colnames(marker_OTU_dada2) <- sapply(strsplit(as.character(colnames(marker_OTU_dada2)),"_H"), `[`, 1)

## TAXONOMY 
# These files were previously created using the script phyloTAGs_preprocessing.R and deduplicated
#ssu_dada2_taxcomplete <- read.csv("/lsalazar/proyectos/phylotags/Taxonomy/ssu_dada2_taxcomplete.txt",  header=F)
#colnames(ssu_dada2_taxcomplete)<- c("seq_id","query","gtdb","domain","phylum","class","order","family","genus","species")
#write.csv(file= "/lsalazar/proyectos/phylotags/Taxonomy/ssu_dada2_taxcomplete_ed.txt",x=ssu_dada2_taxcomplete[which(!duplicated(ssu_dada2_taxcomplete)),],row.names = F)

#marker_dada2_taxcomplete<- read.delim("marker_dada2_taxcomplete_ed.txt", sep=",", header=T)
#marker_dada2_gtdb_updated <- read.delim("marker_dada2_gtdb_taxonomy.tsv", sep="\t", header=T)
#colnames(marker_dada2_gtdb_updated) <- c("seq_id","query","domain","phylum","class","order","family","genus","species")
#marker_dada2_taxcomplete <- merge(marker_dada2_taxcomplete_nas,marker_dada2_gtdb_updated,by="seq_id")

marker_dada2_taxcomplete <- read.delim("marker_dada2_gtdb_taxonomy.tsv", sep="\t", header=T)

#convert to matrix
marker_taxmat_dada2 = as.matrix(marker_dada2_taxcomplete)
rownames(marker_taxmat_dada2) <- marker_dada2_taxcomplete[,1]
marker_TAX_dada2 = tax_table(marker_taxmat_dada2) #[,-ncol(dnak_taxmat_dada2)]
#colnames(dnak_TAX_dada2) <- c("lib","query","kingdom","phylum","class","order","family","genus","species")

## TREES
marker_qza_dada2 <- read_qza("marker_rooted_tree_Dada2.qza")


## SEQs
library(Biostrings)

marker_seq <- readDNAStringSet("marker_RepSeq_Dada2.fasta") 

## PHYLOSEQ
marker_physeq_dada2 = phyloseq(marker_OTU_dada2, marker_TAX_dada2, sample_data(microbio_selected), marker_qza_dada2$data, marker_seq)

#filter out asv that are only found once 
marker_physeq <- phyloseq::subset_samples(marker_physeq_dada2, phyloseq::sample_sums(marker_physeq_dada2) > 100)
marker_physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(marker_physeq) > 2, marker_physeq)

# Rarefaction
marker_rare <- phyloseq::rarefy_even_depth(marker_physeq, rngseed = 123, replace = FALSE)

## DESeq 
marker_physeq <- prune_taxa(taxa_sums(marker_physeq_dada2) > 1, marker_physeq_dada2)
marker_physeq <- prune_samples(sample_sums(marker_physeq_dada2)>0, marker_physeq_dada2)

marker_dds = phyloseq_to_deseq2(marker_physeq, ~ bmi_class)
marker_dds.bmi <- estimateSizeFactors(marker_dds,type="poscounts") # solves zero-inflated counts
marker_dds.wald= DESeq(marker_dds.bmi, fitType="local")

#NW vs Ob
alpha = 0.05

marker_res.nw_ob <- results(marker_dds.wald,contrast=c("bmi_class","Obese","Lean"))
marker_res.nw_ob = marker_res.nw_ob[order(marker_res.nw_ob$padj, na.last=NA), ]
marker_sigtab.nw_ob = marker_res.nw_ob[(marker_res.nw_ob$padj < alpha), ]
marker_sigtab.nw_ob = cbind(as(marker_sigtab.nw_ob, "data.frame"), as(tax_table(marker_physeq)[rownames(marker_sigtab.nw_ob), ], "matrix"))
