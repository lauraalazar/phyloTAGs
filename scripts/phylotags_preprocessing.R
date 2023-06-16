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

############
# METADATA #
############

#Metadata files with  info for each library
ssu_meta = read.table("metadata_libs/16S.files", header = F)
dnak_meta = read.table("metadata_libs/dnak1.files", header = F)
gyrB_meta = read.table("metadata_libs/gyrB.files", header = F)

#edit and remove library names
ssu_meta$lib <- sapply(strsplit(as.character(ssu_meta[,2]), "\\_"), `[`, 1)
dnak_meta$lib <- sapply(strsplit(as.character(dnak_meta[,2]), "\\_"), `[`, 1)
gyrB_meta$lib <- sapply(strsplit(as.character(gyrB_meta[,2]), "\\_"), `[`, 1)

#metadata bmi
sample_info <- read.table("bsp_classification.txt",header=T)
rownames(sample_info) <- sapply(strsplit(as.character(sample_info$ID),"_H"), `[`, 1)

#metadata for design variables
microbio_selected <- read.table("metadata_libs/microbio_selected.meta.txt",header=T)
#remove the seconde records (records with "H2") -> el id "MI_093_H12" est? mal anotado fue manualmente corregido
microbio_selected <- microbio_selected[which(sapply(strsplit(as.character(microbio_selected$ID),"_"), `[`, 3)!="H2"),]
rownames(microbio_selected) <- sapply(strsplit(as.character(microbio_selected$ID),"_H"), `[`, 1)
microbio_selected$bmi_class <- factor(sapply(strsplit(as.character(microbio_selected$bmi_class),"-"), `[`, 1),levels = c("Lean","Overwe$

#filter to keep only the samples used in phyloTags

###################
# FEATURES-COUNTS #
###################

#Feature files: counts of library for each sample
#these files start with "#" wich R does not read because interprets it as comment
#Read the header extract it and then name the columns from the file

ssu_feature_dada2 <- read.table("qiime_dada_deblur_complete/16S_Tabla_Dada2.txt", header = F)
ssu_con <- file("qiime_dada_deblur_complete/16S_Tabla_Dada2.txt", "r")
ssu_header <- readLines(ssu_con, n=2)
colnames(ssu_feature_dada2) <- strsplit(ssu_header, "\\s+")[[2]][-1]
close(ssu_con)

dnak_feature_dada2 <- read.table("qiime_dada_deblur_complete/dnaK_Dada2_Table.txt", header = F)
dnak_con <- file("qiime_dada_deblur_complete/dnaK_Dada2_Table.txt", "r")
dnak_header_dada2 <- readLines(dnak_con, n=2)
colnames(dnak_feature_dada2) <- strsplit(dnak_header_dada2, "\\s+")[[2]][-1]
close(dnak_con)

gyrB_feature_dada2 <- read.table("qiime_dada_deblur_complete/gyrB_Dada2_Table.txt", header = F)
gyrB_con <- file("qiime_dada_deblur_complete/gyrB_Dada2_Table.txt", "r")
gyrB_header_dada2 <- readLines(gyrB_con, n=2)
colnames(gyrB_feature_dada2) <- strsplit(gyrB_header_dada2, "\\s+")[[2]][-1]   
close(gyrB_con)

#asv table of dada2
ssu_mat_dada2 = as.matrix(ssu_feature_dada2[,-1])
rownames(ssu_mat_dada2) <- ssu_feature_dada2[,1]

dnak_mat_dada2 = as.matrix(dnak_feature_dada2[,-1])
rownames(dnak_mat_dada2) <- dnak_feature_dada2[,1]

gyrB_mat_dada2 = as.matrix(gyrB_feature_dada2[,-1])
rownames(gyrB_mat_dada2) <- gyrB_feature_dada2[,1]
                                                                                                                         
#OTU tables dada2
ssu_OTU_dada2 = otu_table(ssu_mat_dada2, taxa_are_rows = TRUE)
colnames(ssu_OTU_dada2) <- ssu_meta[match(colnames(ssu_OTU_dada2),ssu_meta$lib),]$V1
colnames(ssu_OTU_dada2) <- gsub("-", "_", colnames(ssu_OTU_dada2)) #change lib names for sample names (low stripes)
colnames(ssu_OTU_dada2) <- sapply(strsplit(as.character(colnames(ssu_OTU_dada2)),"_H"), `[`, 1)

dnak_OTU_dada2 = otu_table(dnak_mat_dada2, taxa_are_rows = TRUE)
colnames(dnak_OTU_dada2) <- dnak_meta[match(colnames(dnak_OTU_dada2),dnak_meta$lib),1] #change lib names for sample names
dnak_OTU_dada2<-dnak_OTU_dada2[,-which(colnames(dnak_OTU_dada2)=="Bacterial_mix")] #remove bacterial mix

gyrB_OTU_dada2 = otu_table(gyrB_mat_dada2, taxa_are_rows = TRUE)
colnames(gyrB_OTU_dada2) <- gyrB_meta[match(colnames(gyrB_OTU_dada2),gyrB_meta$lib),1] #change lib names for sample names
gyrB_OTU_dada2<-gyrB_OTU_dada2[,-which(colnames(gyrB_OTU_dada2)=="Mix_bacteriano")] #remove bacterial mix
colnames(gyrB_OTU_dada2) <- sapply(strsplit(as.character(colnames(gyrB_OTU_dada2)),"_H"), `[`, 1)

############
# TAXONOMY #
############
# These files were previously created using the script phyloTAGs_preprocessing.R and deduplicated
#ssu_dada2_taxcomplete <- read.csv("/lsalazar/proyectos/phylotags/Taxonomy/ssu_dada2_taxcomplete.txt",  header=F)  
#colnames(ssu_dada2_taxcomplete)<- c("seq_id","query","gtdb","domain","phylum","class","order","family","genus","species")
#write.csv(file= "/lsalazar/proyectos/phylotags/Taxonomy/ssu_dada2_taxcomplete_ed.txt",x=ssu_dada2_taxcomplete[which(!duplicated(ssu_da$

#ssu_dada2_taxcomplete<- read.delim("/lsalazar/proyectos/phylotags/Taxonomy/ssu_dada2_taxcomplete_ed.txt", sep=",", header=T)
#ssu_dada2_gtdb_updated <- read.delim("/lsalazar/proyectos/phylotags/Taxonomy/16S_dada2_gtdb_taxonomy.tsv", sep="\t", header=T)
#colnames(ssu_dada2_gtdb_updated) <- c("seq_id","query","domain","phylum","class","order","family","genus","species")
#ssu_dada2_taxcomplete <- merge(ssu_dada2_taxcomplete_nas,ssu_dada2_gtdb_updated,by="seq_id")

ssu_dada2_taxcomplete <- read.delim("Taxonomy/16S_dada2_gtdb_taxonomy.tsv", sep="\t", header=T)
dnak_dada2_taxcomplete <- read.csv("Taxonomy/dnak_dada2_taxcomplete_ed.txt", header=T)
gyrB_dada2_taxcomplete <- read.csv("Taxonomy/gyrB_dada2_taxcomplete_ed.txt", header=T)

#convert to matrix
ssu_taxmat_dada2 = as.matrix(ssu_dada2_taxcomplete)
rownames(ssu_taxmat_dada2) <- ssu_dada2_taxcomplete[,1]
ssu_TAX_dada2 = tax_table(ssu_taxmat_dada2) #[,-ncol(ssu_taxmat_dada2)]
colnames(ssu_TAX_dada2) <- c("seq_id","query","domain","phylum","class","order","family","genus","species")

dnak_taxmat_dada2 = as.matrix(dnak_dada2_taxcomplete)
rownames(dnak_taxmat_dada2) <- dnak_dada2_taxcomplete[,1]
dnak_TAX_dada2 = tax_table(dnak_taxmat_dada2) #[,-ncol(dnak_taxmat_dada2)]
                                                                                                                         
gyrB_taxmat_dada2 = as.matrix(gyrB_dada2_taxcomplete)
rownames(gyrB_taxmat_dada2) <- gyrB_dada2_taxcomplete[,1]
gyrB_TAX_dada2 = tax_table(gyrB_taxmat_dada2)  #remove the last column of db -> NO necessary anymore!

#########
# TREES #
#########
ssu_nwk_dada2 <- read.tree("Trees/16S_rooted_tree_Dada2.nwk")   #16S_RepSeq_Dada2.qza")
dnak_qza_dada2 <- read_qza("Trees/dnaK_rooted_tree_Dada2.qza")
gyrB_qza_dada2 <- read_qza("Trees/gyrB_rooted_tree_Dada2.qza")

#Samples metadaaft
#sample_info <- read.table("bsp_classification.txt",header=T)
#rownames(sample_info) <- sapply(strsplit(as.character(sample_info$ID),"_H"), `[`, 1)
# take the subsamples that were wequenced
#subsample <- sample_data(sample_info %>%
#                           filter(rownames(sample_info) %in% colnames(dnak_OTU)))

########
# SEQs #
########
library(Biostrings)

#ssu_seq <- readDNAStringSet("/lsalazar/proyectos/phylotags/6_CD-HIT_16S/6_CD-HIT/16S_dna-sequences_100.fasta")
ssu_seq <- readDNAStringSet("struo/input_dada2/16S_RepSeq_Dada2.fasta") #chequear cuales tiene Luis
dnak_seq <- readDNAStringSet("qiime_dada_deblur_complete/Secuencias/Secuencias/dnaK1_Dada2_RepSeq.fasta")
gyrB_seq <- readDNAStringSet("qiime_dada_deblur_complete/Secuencias/Secuencias/gyrB_Dada2_RepSeq.fasta")

############
# PHYLOSEQ #
############
#metadata bmi
#ssu_physeq_dada2 = phyloseq(ssu_OTU_dada2, ssu_TAX_dada2, sample_data(sample_info), ssu_qza_dada2$data)
#dnak_physeq_dada2 = phyloseq(dnak_OTU_dada2, dnak_TAX_dada2, sample_data(sample_info), dnak_qza_dada2$data)
#gyrB_physeq_dada2 = phyloseq(gyrB_OTU_dada2, gyrB_TAX_dada2, sample_data(sample_info), gyrB_qza_dada2$data)

ssu_physeq_dada2 = phyloseq(ssu_OTU_dada2, ssu_TAX_dada2, sample_data(microbio_selected), ssu_nwk_dada2, ssu_seq)
dnak_physeq_dada2 = phyloseq(dnak_OTU_dada2, dnak_TAX_dada2, sample_data(microbio_selected), dnak_qza_dada2$data, dnak_seq)
gyrB_physeq_dada2 = phyloseq(gyrB_OTU_dada2, gyrB_TAX_dada2, sample_data(microbio_selected), gyrB_qza_dada2$data, gyrB_seq)

#filter out asv that are only found once
ssu_physeq <- phyloseq::subset_samples(ssu_physeq_dada2, phyloseq::sample_sums(ssu_physeq_dada2) > 100)  
ssu_physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(ssu_physeq) > 2, ssu_physeq)

dnak_physeq <- phyloseq::subset_samples(dnak_physeq_dada2, phyloseq::sample_sums(dnak_physeq_dada2) > 100)
dnak_physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(dnak_physeq) > 2, dnak_physeq)

gyrB_physeq <- phyloseq::subset_samples(gyrB_physeq_dada2, phyloseq::sample_sums(gyrB_physeq_dada2) > 100)
gyrB_physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(gyrB_physeq) > 2, gyrB_physeq)

# Rarefaction
ssu_rare <- phyloseq::rarefy_even_depth(ssu_physeq, rngseed = 123, replace = FALSE)
dnak_rare <- phyloseq::rarefy_even_depth(dnak_physeq, rngseed = 123, replace = FALSE)
gyrB_rare <- phyloseq::rarefy_even_depth(gyrB_physeq, rngseed = 123, replace = FALSE)
                                                                                                                         
#########
# DESEQ #
#########

ssu_physeq <- prune_taxa(taxa_sums(ssu_physeq_dada2) > 1, ssu_physeq_dada2)
ssu_physeq <- prune_samples(sample_sums(ssu_physeq_dada2)>0, ssu_physeq_dada2)

dnak_physeq <- prune_taxa(taxa_sums(dnak_physeq_dada2) > 1, dnak_physeq_dada2)
dnak_physeq <- prune_samples(sample_sums(dnak_physeq_dada2)>0,dnak_physeq_dada2)

gyrB_physeq <- prune_taxa(taxa_sums(gyrB_physeq_dada2) > 1, gyrB_physeq_dada2)
gyrB_physeq <- prune_samples(sample_sums(gyrB_physeq_dada2)>0,gyrB_physeq_dada2)
                                                                                                                         
ssu_dds = phyloseq_to_deseq2(ssu_physeq, ~ bmi_class)
ssu_dds.bmi <- estimateSizeFactors(ssu_dds,type="poscounts") # solves zero-inflated counts
ssu_dds.wald= DESeq(ssu_dds.bmi, fitType="local")

dnak_dds = phyloseq_to_deseq2(dnak_physeq, ~ bmi_class)
dnak_dds.bmi <- estimateSizeFactors(dnak_dds,type="poscounts") # solves zero-inflated counts
dnak_dds.wald= DESeq(dnak_dds.bmi, fitType="local")

gyrB_dds = phyloseq_to_deseq2(gyrB_physeq, ~ bmi_class)
gyrB_dds.bmi <- estimateSizeFactors(gyrB_dds,type="poscounts") # solves zero-inflated counts
gyrB_dds.wald <- DESeq(gyrB_dds.bmi, fitType="local")

                                                                                                                         
#NW vs Ob
alpha = 0.05

ssu_res.nw_ob <- results(ssu_dds.wald,contrast=c("bmi_class","Obese","Lean"))
ssu_res.nw_ob = ssu_res.nw_ob[order(ssu_res.nw_ob$padj, na.last=NA), ]
ssu_sigtab.nw_ob = ssu_res.nw_ob[(ssu_res.nw_ob$padj < alpha), ]
ssu_sigtab.nw_ob = cbind(as(ssu_sigtab.nw_ob, "data.frame"), as(tax_table(ssu_physeq)[rownames(ssu_sigtab.nw_ob), ], "matrix"))

dnak_res.nw_ob <- results(dnak_dds.wald,contrast=c("bmi_class","Obese","Lean"))
dnak_res.nw_ob = dnak_res.nw_ob[order(dnak_res.nw_ob$padj, na.last=NA), ]
dnak_sigtab.nw_ob = dnak_res.nw_ob[(dnak_res.nw_ob$padj < alpha), ]
dnak_sigtab.nw_ob = cbind(as(dnak_sigtab.nw_ob, "data.frame"), as(tax_table(dnak_physeq)[rownames(dnak_sigtab.nw_ob), ], "matrix"))

gyrB_res.nw_ob <- results(gyrB_dds.wald,contrast=c("bmi_class","Obese","Lean"))
gyrB_res.nw_ob = gyrB_res.nw_ob[order(gyrB_res.nw_ob$padj, na.last=NA), ]
gyrB_sigtab.nw_ob = gyrB_res.nw_ob[(gyrB_res.nw_ob$padj < alpha), ]
gyrB_sigtab.nw_ob = cbind(as(gyrB_sigtab.nw_ob, "data.frame"), as(tax_table(gyrB_physeq)[rownames(gyrB_sigtab.nw_ob), ], "matrix"))


                                                                                                                         
