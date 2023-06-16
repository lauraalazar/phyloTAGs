#source("/lsalazar/proyectos/phylotags/articulo/scripts/phylotags_preprocessing.R")

#Supplementary Figure1
fig5C <- readRDS("/lsalazar/estudiantes/Maria/reporte_proyecto_marcadores/data_scripts/figure5C_16s_dnak1_fam.rds")
fig6C <- readRDS("/lsalazar/estudiantes/Maria/reporte_proyecto_marcadores/data_scripts/figure6C_16s_gyrb_fam.rds")
fig7C <- readRDS("/lsalazar/estudiantes/Maria/reporte_proyecto_marcadores/data_scripts/figure7C_dnak1_gyrb_fam.rds")

suppcomparison_title <- ggdraw() + draw_label("Pairwise comparison among Clostridia genome assemblies within families", fontface='bold')
suppmarkers.plot <- plot_grid(fig5C, fig6C, fig7C,   ncol=1, nrow=3)
suppcomparison.plot <- plot_grid(suppcomparison_title,suppmarkers.plot, ncol=1,rel_heights=c(0.1, 1))

svg("/lsalazar/proyectos/phylotags/articulo/msystems/figures/figure_S1.svg")
suppcomparison.plot
dev.off()

# Test results when clustering at 97%
#ssu 
ssu_refseq <- refseq(ssu_physeq_dada2) 
nproc <- 1 # Increase to use multiple processors
ssu_aln <- DECIPHER::AlignSeqs(ssu_refseq, processors = nproc)
ssu_d <- DECIPHER::DistanceMatrix(ssu_aln, processors = nproc)


#dna 
dnak_refseq <- refseq(dnak_physeq_dada2)
nproc <- 1 # Increase to use multiple processors
dnak_aln <- DECIPHER::AlignSeqs(dnak_refseq, processors = nproc)
dnak_d <- DECIPHER::DistanceMatrix(dnak_aln, processors = nproc)

#gyrB 
gyrB_refseq <- refseq(gyrB_physeq_dada2)
nproc <- 1 # Increase to use multiple processors
gyrB_aln <- DECIPHER::AlignSeqs(gyrB_refseq, processors = nproc)
gyrB_d <- DECIPHER::DistanceMatrix(gyrB_aln, processors = nproc)

#Clusters
ssu_clusters_97 <- DECIPHER::IdClusters(
  ssu_d, 
  method = "complete",
  cutoff = 0.03, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)

ssu_clusters_90 <- DECIPHER::IdClusters(
  ssu_d, 
  method = "complete",
  cutoff = 0.1, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)


ssu_clusters_85 <- DECIPHER::IdClusters(
  ssu_d, 
  method = "complete",
  cutoff = 0.15, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)


dnak_clusters_97 <- DECIPHER::IdClusters(
  dnak_d, 
  method = "complete",
  cutoff = 0.03, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)

dnak_clusters_90 <- DECIPHER::IdClusters(
  dnak_d, 
  method = "complete",
  cutoff = 0.1, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)


dnak_clusters_85 <- DECIPHER::IdClusters(
  dnak_d, 
  method = "complete",
  cutoff = 0.15, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)

gyrB_clusters_97 <- DECIPHER::IdClusters(
  gyrB_d, 
  method = "complete",
  cutoff = 0.03, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)

gyrB_clusters_90 <- DECIPHER::IdClusters(
  gyrB_d, 
  method = "complete",
  cutoff = 0.1, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)


gyrB_clusters_85 <- DECIPHER::IdClusters(
  gyrB_d, 
  method = "complete",
  cutoff = 0.15, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Maria)
  processors = nproc
)

# Use speedyseq to merge taxa in phyloseq object, as explained in:
# https://github.com/mikemc/speedyseq/blob/main/NEWS.md#new-general-purpose-vectorized-merging-function 
library(speedyseq)

#ssu
ssu_physeq_clust85 <- merge_taxa_vec(ssu_physeq_dada2, group=ssu_clusters_85$cluster, tax_adjust = 0)
ssu_clusters_fulltax_clust85 <- merge(ssu_clusters_85,as.data.frame(ssu_TAX_dada2), by="row.names")
ssu_rare_85 <- phyloseq::rarefy_even_depth(ssu_physeq_clust85, rngseed = 123, replace = FALSE)

ssu_physeq_clust90 <- merge_taxa_vec(ssu_physeq_dada2, group=ssu_clusters_90$cluster, tax_adjust = 0)
ssu_clusters_fulltax_clust90 <- merge(ssu_clusters_90,as.data.frame(ssu_TAX_dada2), by="row.names")
ssu_rare_90 <- phyloseq::rarefy_even_depth(ssu_physeq_clust90, rngseed = 123, replace = FALSE)

ssu_physeq_clust97 <- merge_taxa_vec(ssu_physeq_dada2, group=ssu_clusters_97$cluster, tax_adjust = 0)
ssu_clusters_fulltax_clust97 <- merge(ssu_clusters_97,as.data.frame(ssu_TAX_dada2), by="row.names")
ssu_rare_97 <- phyloseq::rarefy_even_depth(ssu_physeq_clust97, rngseed = 123, replace = FALSE)

ssu_shannon.fam <- marker_tax_stats("ssu",85,"family")$shannon
ssu_shannon.gen <- marker_tax_stats("ssu",85,"genus")$shannon


#dnaK
dnak_physeq_clust85 <- merge_taxa_vec(dnak_physeq_dada2, group=dnak_clusters_85$cluster, tax_adjust = 0)
dnak_clusters_fulltax_clust85 <- merge(dnak_clusters_85,as.data.frame(dnak_TAX_dada2), by="row.names")
dnak_rare_85 <- phyloseq::rarefy_even_depth(dnak_physeq_clust85, rngseed = 123, replace = FALSE)

dnak_physeq_clust90 <- merge_taxa_vec(dnak_physeq_dada2, group=dnak_clusters_90$cluster, tax_adjust = 0)
dnak_clusters_fulltax_clust90 <- merge(dnak_clusters_90,as.data.frame(dnak_TAX_dada2), by="row.names")
dnak_rare_90 <- phyloseq::rarefy_even_depth(dnak_physeq_clust90, rngseed = 123, replace = FALSE)

dnak_physeq_clust97 <- merge_taxa_vec(dnak_physeq_dada2, group=dnak_clusters_97$cluster, tax_adjust = 0)
dnak_clusters_fulltax_clust97 <- merge(dnak_clusters_97,as.data.frame(dnak_TAX_dada2), by="row.names")
dnak_rare_97 <- phyloseq::rarefy_even_depth(dnak_physeq_clust97, rngseed = 123, replace = FALSE)

#gyrB
gyrB_physeq_clust85 <- merge_taxa_vec(gyrB_physeq_dada2, group=gyrB_clusters_85$cluster, tax_adjust = 0)
gyrB_clusters_fulltax_clust85 <- merge(gyrB_clusters_85,as.data.frame(gyrB_TAX_dada2), by="row.names")
gyrB_rare_85 <- phyloseq::rarefy_even_depth(gyrB_physeq_clust85, rngseed = 123, replace = FALSE, sample.size = 82)

gyrB_physeq_clust90 <- merge_taxa_vec(gyrB_physeq_dada2, group=gyrB_clusters_90$cluster, tax_adjust = 0)
gyrB_clusters_fulltax_clust90 <- merge(gyrB_clusters_90,as.data.frame(gyrB_TAX_dada2), by="row.names")
gyrB_rare_90 <- phyloseq::rarefy_even_depth(gyrB_physeq_clust90, rngseed = 123, replace = FALSE, sample.size = 82)

gyrB_physeq_clust97 <- merge_taxa_vec(gyrB_physeq_dada2, group=gyrB_clusters_97$cluster, tax_adjust = 0)
gyrB_clusters_fulltax_clust97 <- merge(gyrB_clusters_97,as.data.frame(gyrB_TAX_dada2), by="row.names")
pruned <- prune_samples(sample_sums(gyrB_clusters_fulltax_clust97) > 0, gyrB_clusters_fulltax_clust97,taxa_are_rows=TRUE)
gyrB_rare_97 <- phyloseq::rarefy_even_depth(gyrB_physeq_clust97, rngseed = 123, replace = FALSE, sample.size = 82)

gyrB_shannon.fam <- marker_tax_stats("gyrB",97,"family")$shannon
gyrB_shannon.gen <- marker_tax_stats("gyrB",97,"genus")$shannon


gyrB_shannon.fam <- marker_tax_stats("gyrB",85,"family")$shannon
gyrB_shannon.gen <- marker_tax_stats("gyrB",85,"genus")$shannon


# From here on continue with the code of phylotags_diversity.R section #ANOVA forward
# Id 97%
idlevel=97
ssu_rare_97 <- ssu_rare
dnak_rare_97 <- dnak_rare
gyrB_rare_97 <- gyrB_rare
#Sig families
ssu_glmcont_city_fam.sigtab
gyrB_glmcont_city_fam.sigtab #lachno
dnak_glmcont_city_fam.sigtab #oscillo



#phylotags_plots figure2
family_title <- ggdraw() + draw_label("Siginificant families at 97% clustering", fontface='bold')
sigfam.plot <- plot_grid(gyrB_fam.plot, oscillo_fam.plot, ncol=2, nrow=1)
family.plot <- plot_grid(family_title,sigfam.plot, ncol=1,rel_heights=c(0.1, 1))

pdf("/lsalazar/proyectos/phylotags/articulo/msystems/figures/sig_fam_id97.pdf",height=5,width=15)
family.plot
dev.off()


#Sig genera cont
gyrB_glmcont_city_gen.sigtab #g__Ruminococcus_A
dnak_glmcont_city_gen.sigtab #g__CAG_177

#Sig genera bin
dnak_glmbin_gen.sigtab #CAG_170,CAG_488,CAG_83,ER4
gyrB_glmbin_gen.sigtab #UBA11524

rumino <- gyrB_shannon.gen %>% 
  select(bmi,shannon_cont = "g__Ruminococcus_A") %>% 
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

rumino_cont.plot <-ggplot(rumino, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,1.5)) + 
  xlim(c(18,40)) + 
  geom_line(data=subset(rumino,shannon_cont>0), 
            aes(y = fitted(gyrB_glmcont_gen.df["g__Ruminococcus_A"][[1]])), size=.75) + 
  annotate("text", x = 27, y = 1.5, label = "cont.pval=0.0188") +
  ggtitle("gyrB: Ruminococcus_A (Lachnospiraceae)") +
  labs(x = bquote("BMI ("~m/s^2 ~")"), y = "Shannon\n")

cag170 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_170") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag170_bin.plot <-ggplot(cag170)+
  geom_point(size=2.5,aes(bmi, shannon_cont , color= labelling)) + #ifelse(0, '0', '>0')
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,2)) + 
  xlim(c(18,40)) + 
  geom_smooth(data = cag170, aes(x = bmi, y = shannon_bin),
              method = "glm", method.args = list(family = "binomial"), 
              se = FALSE, fullrange=TRUE, color = okabe[1]) + 
  annotate("text", x = 30, y = 2, label = "bin.pval=0.000853") +
  ggtitle("dnaK: CAG-170 (Oscillospiraceae)") +
  labs(x = bquote("BMI ("~m/s^2 ~")"), y = "Shannon\n")

cag83 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_83") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag83_bin.plot <-ggplot(cag83)+
#  geom_point(size=2.5,aes(bmi, shannon_cont , color= labelling)) +
  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,1)) + 
  xlim(c(18,40)) + 
  geom_smooth(data = cag83, aes(x = bmi, y = shannon_bin),
              method = "glm", method.args = list(family = "binomial"), 
              se = FALSE, fullrange=TRUE, color = okabe[1]) + 
  annotate("text", x = 30, y = 0.75, label = "bin.pval=0.0158") +
  ggtitle("dnaK: CAG-83 (Oscillospiraceae)") +
  labs(x = bquote("BMI ("~Kg/s^2 ~")"), y = "Probability of Shannon > 0 \n")


er4 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__ER4") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

er4_bin.plot <-ggplot(er4)+
  geom_point(size=2.5,aes(bmi, shannon_cont , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) + 
  xlim(c(18,40)) + 
  geom_smooth(data = er4, aes(x = bmi, y = shannon_bin),
              method = "glm", method.args = list(family = "binomial"), 
              se = FALSE, fullrange=TRUE, color = okabe[1]) + 
  annotate("text", x = 30, y = 4, label = "bin.pval=0.0308") +
  ggtitle("dnaK: ER4 (Oscillospiraceae)") +
  labs(x = bquote("BMI ("~m/s^2 ~")"), y = "Shannon\n")


cag488 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_488") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag488_bin.plot <-ggplot(cag488)+
  geom_point(size=2.5,aes(bmi, shannon_cont , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) + 
  xlim(c(18,40)) + 
  geom_smooth(data = cag488, aes(x = bmi, y = shannon_bin),
              method = "glm", method.args = list(family = "binomial"), 
              se = FALSE, fullrange=TRUE, color = okabe[1]) + 
  annotate("text", x = 30, y = 4, label = "bin.pval=0.0127") +
  ggtitle("dnaK: CAG-488 (Acutalibacteraceae)") +
  labs(x = bquote("BMI ("~m/s^2 ~")"), y = "Shannon\n")

sigen.plot <- plot_grid(rumino_cont.plot,cag488_bin.plot,cag170_bin.plot, cag83_bin.plot,er4_bin.plot,ncol=3)

genera_title <- ggdraw() + draw_label("Significant genera at 97% clustering", fontface='bold')

genera.plot <- plot_grid(genera_title,sigen.plot, ncol=1,rel_heights=c(0.1, 1))

pdf("/lsalazar/proyectos/phylotags/articulo/msystems/figures/sig_gen_id97.pdf",height=5,width=15)
genera.plot
dev.off()

# Id 97%
idlevel=85 # from Mar?a's data

#Supplementary Figure Differential Abundance 


#NW vs Ob
alpha = 0.05

ssu_res.nw_ow <- results(ssu_dds.wald,contrast=c("bmi_class","Overweight","Lean"))
ssu_res.nw_ow = ssu_res.nw_ow[order(ssu_res.nw_ow$padj, na.last=NA), ]
ssu_sigtab.nw_ow = ssu_res.nw_ow[(ssu_res.nw_ow$padj < alpha), ]
ssu_sigtab.nw_ow = cbind(as(ssu_sigtab.nw_ow, "data.frame"), as(tax_table(ssu_physeq)[rownames(ssu_sigtab.nw_ow), ], "matrix"))

dnak_res.nw_ow <- results(dnak_dds.wald,contrast=c("bmi_class","Overweight","Lean"))
dnak_res.nw_ow = dnak_res.nw_ow[order(dnak_res.nw_ow$padj, na.last=NA), ]
dnak_sigtab.nw_ow = dnak_res.nw_ow[(dnak_res.nw_ow$padj < alpha), ]
dnak_sigtab.nw_ow = cbind(as(dnak_sigtab.nw_ow, "data.frame"), as(tax_table(dnak_physeq)[rownames(dnak_sigtab.nw_ow), ], "matrix"))

gyrB_res.nw_ow <- results(gyrB_dds.wald,contrast=c("bmi_class","Overweight","Lean"))
gyrB_res.nw_ow = gyrB_res.nw_ow[order(gyrB_res.nw_ow$padj, na.last=NA), ]
gyrB_sigtab.nw_ow = gyrB_res.nw_ow[(gyrB_res.nw_ow$padj < alpha), ]
gyrB_sigtab.nw_ow = cbind(as(gyrB_sigtab.nw_ow, "data.frame"), as(tax_table(gyrB_physeq)[rownames(gyrB_sigtab.nw_ow), ], "matrix"))

write.table(dnak_sigtab.nw_ob,"/lsalazar/proyectos/phylotags/articulo/msystems/tables/dnak_sigtab.nw_ob.txt")
write.table(gyrB_sigtab.nw_ob,"/lsalazar/proyectos/phylotags/articulo/msystems/tables/gyrB_sigtab.nw_ob.txt")


# Fig Prevalence vs Counts
#Prevalence: 
# the number of times an ASV is observed at least once. That is, it is the number of samples in which each OTU was non-zero)
# http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#taxa-prevalence-histogram-and-fast_melt
source("taxa_summary.R", local = TRUE)

ssu_mdt <- fast_melt(ssu_physeq_dada2)
dnak_mdt <- fast_melt(dnak_physeq_dada2)
gyrB_mdt <- fast_melt(gyrB_physeq_dada2)


#Table of prevalence and counts
ssu_prevdt = ssu_mdt[, list(Prevalence = sum(count > 0), 
                            TotalCounts = sum(count)),
                     by = TaxaID]

dnak_prevdt = dnak_mdt[, list(Prevalence = sum(count > 0), 
                              TotalCounts = sum(count)),
                       by = TaxaID]

gyrB_prevdt = gyrB_mdt[, list(Prevalence = sum(count > 0), 
                              TotalCounts = sum(count)),
                       by = TaxaID]


ssu_addGenus = unique(copy(ssu_mdt[, list(TaxaID, genus)]))
dnak_addGenus = unique(copy(dnak_mdt[, list(TaxaID, genus)]))
gyrB_addGenus = unique(copy(gyrB_mdt[, list(TaxaID, genus)]))

# Join by TaxaID
setkey(ssu_addGenus, TaxaID)
ssu_prevdt <- ssu_addGenus[ssu_prevdt]
ssu_showGenus = ssu_prevdt[, sum(TotalCounts), by = genus][order(-V1)][1:10]$genus
setkey(ssu_prevdt, genus)

setkey(dnak_addGenus, TaxaID)
dnak_prevdt <- dnak_addGenus[dnak_prevdt]
dnak_prevdt$genus <-gsub("g__","",dnak_prevdt$genus)
dnak_showGenus = dnak_prevdt[, sum(TotalCounts), by = genus][order(-V1)][1:10]$genus
setkey(dnak_prevdt, genus)

setkey(gyrB_addGenus, TaxaID)
gyrB_prevdt <- gyrB_addGenus[gyrB_prevdt]
gyrB_prevdt$genus <-gsub("g__","",gyrB_prevdt$genus)
gyrB_showGenus = gyrB_prevdt[, sum(TotalCounts), by = genus][order(-V1)][1:10]$genus
setkey(gyrB_prevdt, genus)

# put together all genera to make the color map
allgenera <- c(ssu_showGenus,dnak_showGenus,gyrB_showGenus)

myCol = c("pink1", "violet", "mediumpurple1", "slateblue1", "purple","yellow" ,
          "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",
          "orange", "tomato", "coral2", "palevioletred", "violetred", "red2",
          "springgreen2", "yellowgreen", "palegreen4",
          "wheat2", "tan", "tan2", "tan3", "brown",
          "grey70", "purple3", "black")

names(myCol) <- c(levels(as.factor(ssu_showGenus)),
                  levels(as.factor(dnak_showGenus)),
                  levels(as.factor(gyrB_showGenus)))

colScale <- scale_fill_manual(name = "genus",values = myCol)


plot_grid(ggplot(ssu_prevdt[ssu_showGenus], 
                 mapping = aes(Prevalence, TotalCounts, fill = genus, color=genus)) + geom_point(size = 4, alpha = 0.75) + scale_color_manual(name = "genus",values = myCol) +  theme(legend.title=element_blank()) + scale_y_log10() + ggtitle("16S rRNA"),
          ggplot(dnak_prevdt[dnak_showGenus], 
                 mapping = aes(Prevalence, TotalCounts, fill = genus, color=genus)) + geom_point(size = 4, alpha = 0.75) + scale_color_manual(name = "genus",values = myCol) +   scale_y_log10() + ggtitle("dnaK"),
          ggplot(gyrB_prevdt[gyrB_showGenus], 
                 mapping = aes(Prevalence, TotalCounts, fill = genus, color=genus)) + geom_point(size = 4, alpha = 0.75) + scale_color_manual(name = "genus",values = myCol) +  theme(legend.title=element_blank()) + scale_y_log10() + ggtitle("gyrB"),
          nrow = 3, ncol = 1, label_size = 10)

#Lean vs OW (Supplementary)


# SUPPLEMETARY MATERIAL

### Prevalences vs TotalCounts
ssu_prevdt #redo

### Beta diversity
#Plot for dnaK
pdf("/lsalazar/proyectos/phylotags/reports/beta_clr_allmarkers.pdf")
title <- ggdraw() + draw_label("PC of CLR transformation with euclidean distance", fontface='bold')
bdiv_clr_allmarkers <-plot_grid(
  ssu_bdiv_city,
  ssu_bdiv_bmi,  
  dnak_bdiv_city,
  dnak_bdiv_bmi,
  gyrB_bdiv_city,
  gyrB_bdiv_bmi,
  nrow = 3, ncol = 2, label_size = 10)
plot_grid(title,bdiv_clr_allmarkers, ncol=1, rel_heights=c(0.1, 1))

#Plot unifrac (unweighted and weighted) for city and bmi 
dnak_unif_bmi <- plot_ordination(dnak_rare, dnak_ord_unifrac_un, color = "city") + geom_point(size = 2)
dnak_wunif_bmi <- plot_ordination(dnak_rare, dnak_ord_wunifrac_un, color = "city") + geom_point(size = 2)

#Plots
pdf("/lsalazar/proyectos/phylotags/reports/beta_unifrac_allmarkers.pdf")
title <- ggdraw() + draw_label("Unifrac Distance", fontface='bold')
bdiv_unif_allmarkers <-plot_grid(
  ssu_unif_city,
  ssu_wunif_city,
  ssu_unif_bmi,
  ssu_wunif_bmi,
  dnak_unif_city,
  dnak_wunif_city,
  dnak_unif_bmi,
  dnak_wunif_bmi,
  gyrB_unif_city,
  gyrB_wunif_city,
  gyrB_unif_bmi,
  gyrB_wunif_bmi,
  nrow = 3, ncol = 4, label_size = 10)
plot_grid(title,bdiv_clr_allmarkers, ncol=1, rel_heights=c(0.1, 1))




### Global shannons (for design variables)
plot_grid(
  ssu_adiv  %>%
    gather(key = metric, value = value, c("Observed")) %>% #, "Shannon")) %>%
    mutate(metric = factor(metric, levels = c("Observed"))) %>% # , "Shannon"))) %>%
    ggplot(aes(x = age_range, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = age_range), height = 0, width = .2) +
    scale_colour_brewer(type="qual", palette="Set1") +
    labs(x = "", y = "") +
    theme(legend.position = "none") +
    stat_compare_means(method = "anova", label.y = 4)+
    ggtitle("Observed"), 
  
  dnak_adiv  %>%
    gather(key = metric, value = value, c("Shannon")) %>% #, "Shannon")) %>%
    mutate(metric = factor(metric, levels = c("Shannon"))) %>% # , "Shannon"))) %>%
    ggplot(aes(x = age_range, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = age_range), height = 0, width = .2) +
    scale_colour_brewer(type="qual", palette="Set1") +
    labs(x = "", y = "") +
    stat_compare_means(method = "anova", label.y = 1)+
    theme(legend.position = "none"),  
  #ggtitle("Shannon dnaK"),
  
  dnak_adiv  %>%
    gather(key = metric, value = value, c("Observed")) %>% #, "Shannon")) %>%
    mutate(metric = factor(metric, levels = c("Observed"))) %>% # , "Shannon"))) %>%
    ggplot(aes(x = age_range, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = age_range), height = 0, width = .2) +
    scale_colour_brewer(type="qual", palette="Set1") +
    labs(x = "", y = "") +
    theme(legend.position = "none") +
    stat_compare_means(method = "anova", label.y = 4),
  #ggtitle("Observed dnaK"),
  
  gyrB_adiv  %>%
    gather(key = metric, value = value, c("Shannon")) %>% #, "Shannon")) %>%
    mutate(metric = factor(metric, levels = c("Shannon"))) %>% # , "Shannon"))) %>%
    ggplot(aes(x = age_range, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = age_range), height = 0, width = .2) +
    scale_colour_brewer(type="qual", palette="Set1") +
    labs(x = "", y = "") +
    stat_compare_means(method = "anova", label.y = 1)+
    theme(legend.position = "none"),
  #ggtitle("Shannon dnaK"),
  
  gyrB_adiv  %>%
    gather(key = metric, value = value, c("Observed")) %>% #, "Shannon")) %>%
    mutate(metric = factor(metric, levels = c("Observed"))) %>% # , "Shannon"))) %>%
    ggplot(aes(x = age_range, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = age_range), height = 0, width = .2) +
    scale_colour_brewer(type="qual", palette="Set1") +
    labs(x = "", y = "") +
    theme(legend.position = "none") +
    stat_compare_means(method = "anova", label.y = 4),
  #ggtitle("Observed gyrB"),
  
  labels = c("ssu", "","dnaK","", "gyrB"),
  
  nrow = 3, ncol = 2, label_size = 10)

### Linear models BMI vs Shannon
#source("lsalazar/proyectos/phylotags/articulo/scripts/phylotags_diversity.R")

### Trees


### Differential abundance


### Linear models counts vs diet

### https://mgimond.github.io/Stats-in-R/Logistic.html