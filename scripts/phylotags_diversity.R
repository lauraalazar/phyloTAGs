#rm(list = ls())
#source("phylotags_preprocessing.R")

idlevel=100 #this matters for esting the results at different levels of identity (phylotags_suppl.R)
marker_rare_100 = marker_rare

#ALPHA DIVERSITY
#Generate a data.frame with alpha div measures
ssu_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(ssu_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ssu_rare, measures = "Shannon"),
  #  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(dnak_rare)))), tree = phyloseq::phy_tree(dnak_rare))[, 1],
  #  "CHS" = phyloseq::sample_data(dnak_rare)$chs_class,
  "BMI" = phyloseq::sample_data(ssu_rare)$bmi,
  "bmi_class" = factor(sample_data(ssu_rare)$bmi_class,levels = c("Lean","Overweight","Obese")),
  "city" = phyloseq::sample_data(ssu_rare)$city,
  "sex" = phyloseq::sample_data(ssu_rare)$sex,
  "age_range" = phyloseq::sample_data(ssu_rare)$age_range)
# head(ssu_adiv)

dnak_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(dnak_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(dnak_rare, measures = "Shannon"),
  #  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(dnak_rare)))), tree = phyloseq::phy_tree(dnak_rare))[, 1],
  #  "CHS" = phyloseq::sample_data(dnak_rare)$chs_class,
  "BMI" = phyloseq::sample_data(dnak_rare)$bmi,
  "bmi_class" = factor(sample_data(dnak_rare)$bmi_class,levels = c("Lean","Overweight","Obese")),
  "city" = phyloseq::sample_data(dnak_rare)$city,
  "sex" = phyloseq::sample_data(dnak_rare)$sex,
  "age_range" = phyloseq::sample_data(dnak_rare)$age_range)
# head(dnak_adiv)

gyrB_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(gyrB_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(gyrB_rare, measures = "Shannon"),
  #  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(dnak_rare)))), tree = phyloseq::phy_tree(dnak_rare))[, 1],
  #  "CHS" = phyloseq::sample_data(dnak_rare)$chs_class,
  "BMI" = phyloseq::sample_data(gyrB_rare)$bmi,
  "bmi_class" = factor(sample_data(gyrB_rare)$bmi_class,levels = c("Lean","Overweight","Obese")),
  "city" = phyloseq::sample_data(gyrB_rare)$city,
  "sex" = phyloseq::sample_data(gyrB_rare)$sex,
  "age_range" = phyloseq::sample_data(gyrB_rare)$age_range)
# head(gyrB_adiv)

## SUMMARY OF SHANNON STATS FOR DIFFERENT TAXA

# loop for each family, get the richness and anova, store them to display which are significant with respect to bmi
### Function to get automated phyloseq with tax level and choice within
# https://github.com/joey711/phyloseq/issues/1471
fn <- function(ps, level, choice) {
  x <- paste0(level,"==","'",choice,"'")
  oldTax <- data.frame(tax_table(ps))
  newTax <- subset(oldTax, eval(parse(text=x)))
  tax_table(ps) <- tax_table(as.matrix(newTax))
  return(ps)
}


#Example
fn(dnak_rare, level="family",choice="f__Acidaminococcaceae")

## FUNCTION MARKER_TAX_STATS
#this function takes marker and level of taxonomy as argument and returns
# 1) dataframe of shannon values for each family
# 2) list of anovas 
# 3) list of glm 

 #include idelvel clustering: 100, 97, 85 (for other than 100 run phylotags_suppl.R)

marker_tax_stats = function(marker=NULL,clustering=NULL,level_tax=NULL){
  #take the rarefied phylose object for the marker
  phyobj=get(paste(marker,paste("rare",clustering,sep="_"),sep="_")) 
  #get a list of the families present
  taxa_list = levels(as.factor(tax_table(phyobj)[,level_tax]))
  #initialize a matrix with equal samples as the physeq object
  shannon_family.df = data.frame(matrix(nrow=ncol(otu_table(phyobj)),ncol=0))
  shannonbinary_family.df = data.frame(matrix(nrow=ncol(otu_table(phyobj)),ncol=0))
  observed_family.df = data.frame(matrix(nrow=ncol(otu_table(phyobj)),ncol=0))
  
  #loop over each taxón and calculate de shannon for each sample, then bind to a dataframe
  family_names = NULL
  for (tx in taxa_list){
    #  print(tx)
    shannon_family.df <- cbind(shannon_family.df,phyloseq::estimate_richness(fn(phyobj, level=level_tax,choice=tx),measures = "Shannon"))
    shannonbinary_family.df <- cbind(shannonbinary_family.df,ifelse(phyloseq::estimate_richness(fn(phyobj, level=level_tax,choice=tx),measures = "Shannon")!=0,1,0))
    observed_family.df <- cbind(observed_family.df,phyloseq::estimate_richness(fn(phyobj, level=level_tax,choice=tx),measures = "Observed"))
    
    # fix names that have stripe or space
    tx[grepl("-",tx)] <- str_replace_all(tx, "-","_")
    tx[grepl(" ", tx)] <- str_replace_all(tx, " ","_")
    tx[grepl("^[[:digit:]]+", tx)] <- paste("f",tx[grepl("^[[:digit:]]+", tx)],sep="")
    family_names <- append(family_names,tx)
  }
  
  # Fix names with non permited characters    
  # fix names that have stripe or space
  # family_names[grepl("-",family_names)] <- str_replace_all(tx, "-","_")
  # spaces removed form names
  # family_names[grepl(" ", family_names)] <- family_names[grepl("", family_names)]
  #R does no deal well with strings starting with number, so here I add a letter to those
  # family_names[grepl("^[[:digit:]]+", family_names)] <- paste("f",family_names[grepl("^[[:digit:]]+", family_names)],sep="")
  
  colnames(shannon_family.df) <- family_names 
  colnames(shannonbinary_family.df) <- family_names
  colnames(observed_family.df) <- family_names
  
  #add bmi as explanatory variable
  shannon_family.df$bmi <- sample_data(phyobj)$bmi
  shannonbinary_family.df$bmi <- sample_data(phyobj)$bmi
  observed_family.df$bmi <- sample_data(phyobj)$bmi
  
  shannon_family.df$bmi_class <- factor(sample_data(phyobj)$bmi_class,levels = c("Lean","Overweight","Obese"))
  shannonbinary_family.df$bmi_class <- factor(sample_data(phyobj)$bmi_class,levels = c("Lean","Overweight","Obese"))
  observed_family.df$bmi_class <- factor(sample_data(phyobj)$bmi_class,levels = c("Lean","Overweight","Obese"))
  
  #add city as explanatory variable
  shannon_family.df$city <- factor(sample_data(phyobj)$city)
  shannonbinary_family.df$city <- factor(sample_data(phyobj)$city)
  observed_family.df$city <- factor(sample_data(phyobj)$city)
  
  diversity.df <- list(shannon=shannon_family.df,shannon_bin=shannonbinary_family.df,observed=observed_family.df)
  return(diversity.df)
}

## ANOVA
#ssu
ssu_shannon.fam <- marker_tax_stats("ssu",idlevel,"family")$shannon
#Check the means for each group to know direction 
#marker_shannon.fam  %>% 
#select(CAG_302,bmi_class) %>% 
#  group_by(bmi_class) %>% 
#  summarize(Mean = mean(CAG_302, na.rm=TRUE))

ssu_shannon.gen <- marker_tax_stats("ssu",idlevel,"genus")$shannon

ssu_shannon.fam.aov <- lapply(setdiff(colnames(ssu_shannon.fam), c("bmi","bmi_class","city")),
                              function(s) {aov(as.formula(paste(as.character(s), " ~ bmi_class")),ssu_shannon.fam)})
ssu_shannon.gen.aov <- lapply(setdiff(colnames(ssu_shannon.gen), c("bmi","bmi_class","city")),
                              function(s) {aov(as.formula(paste(as.character(s), " ~ bmi_class")),ssu_shannon.gen)})

#simple way to do it
#lapply(marker_shannon.fam.aov,function(x) {summary(x)})

#more fancy but gives a list of all statistics for all families
ssu_vars.fam = names(ssu_shannon.fam)[!names(ssu_shannon.fam) %in% c("bmi","bmi_class","city")]
ssu_vars.gen = names(ssu_shannon.gen)[!names(ssu_shannon.gen) %in% c("bmi","bmi_class","city")]

ssu_anova_fam.df = lapply(setNames(ssu_vars.fam, ssu_vars.fam),function(s) {
  aov(as.formula(paste(as.character(s), " ~ bmi_class")),ssu_shannon.fam)})

ssu_anova_gen.df = lapply(setNames(ssu_vars.gen, ssu_vars.gen),function(s) {
  aov(as.formula(paste(as.character(s), " ~ bmi_class")),ssu_shannon.gen)})

#significant table
ssu_anova_fam.sigtab <- map_df(ssu_anova_fam.df, tidy, .id="taxa")[which(map_df(ssu_anova_fam.df, tidy, .id="taxa")[,"p.value"]<0.05$
ssu_anova_gen.sigtab <- map_df(ssu_anova_gen.df, tidy, .id="taxa")[which(map_df(ssu_anova_gen.df, tidy, .id="taxa")[,"p.value"]<0.05$

#dnaK
dnak_shannon.fam.aov <- lapply(setdiff(colnames(dnak_shannon.fam), c("bmi","bmi_class","city")),
                               function(s) {aov(as.formula(paste(as.character(s), " ~ bmi_class")),dnak_shannon.fam)})
dnak_shannon.gen.aov <- lapply(setdiff(colnames(dnak_shannon.gen), c("bmi","bmi_class","city")),
                               function(s) {aov(as.formula(paste(as.character(s), " ~ bmi_class")),dnak_shannon.gen)})

dnak_vars.fam = names(dnak_shannon.fam)[!names(dnak_shannon.fam) %in% c("bmi","bmi_class","city")]
dnak_vars.gen = names(dnak_shannon.gen)[!names(dnak_shannon.gen) %in% c("bmi","bmi_class","city")]

dnak_anova_fam.df = lapply(setNames(dnak_vars.fam, dnak_vars.fam),function(s) {
  aov(as.formula(paste(as.character(s), " ~ bmi_class")),dnak_shannon.fam)})
dnak_anova_gen.df = lapply(setNames(dnak_vars.gen, dnak_vars.gen),function(s) {
  aov(as.formula(paste(as.character(s), " ~ bmi_class")),dnak_shannon.gen)})

#significant table
dnak_anova_fam.sigtab <- map_df(dnak_anova_fam.df, tidy, .id="taxa")[which(map_df(dnak_anova_fam.df, tidy, .id="taxa")[,"p.value"]<0$
dnak_anova_gen.sigtab <- map_df(dnak_anova_gen.df, tidy, .id="taxa")[which(map_df(dnak_anova_gen.df, tidy, .id="taxa")[,"p.value"]<0$

#gyrB
#gyrB_rare_100 <- gyrB_rare
gyrB_shannon.fam <- marker_tax_stats("gyrB",idlevel,"family")$shannon
gyrB_shannon.gen <- marker_tax_stats("gyrB",idlevel,"genus")$shannon
  

gyrB_shannon.fam.aov <- lapply(setdiff(colnames(gyrB_shannon.fam), c("bmi","bmi_class","city")),
                               function(s) {aov(as.formula(paste(as.character(s), " ~ bmi_class")),gyrB_shannon.fam)})
gyrB_shannon.gen.aov <- lapply(setdiff(colnames(gyrB_shannon.gen), c("bmi","bmi_class","city")),
                               function(s) {aov(as.formula(paste(as.character(s), " ~ bmi_class")),gyrB_shannon.gen)})

gyrB_vars.fam = names(gyrB_shannon.fam)[!names(gyrB_shannon.fam) %in% c("bmi","bmi_class","city")]
gyrB_vars.gen = names(gyrB_shannon.gen)[!names(gyrB_shannon.gen) %in% c("bmi","bmi_class","city")]


gyrB_anova_fam.df = lapply(setNames(gyrB_vars.fam, gyrB_vars.fam),function(s) {
  aov(as.formula(paste(as.character(s), " ~ bmi_class")),gyrB_shannon.fam)})
gyrB_anova_gen.df = lapply(setNames(gyrB_vars.gen, gyrB_vars.gen),function(s) {
  aov(as.formula(paste(as.character(s), " ~ bmi_class")),gyrB_shannon.gen)})
  
#significant table
gyrB_anova_fam.sigtab <- map_df(gyrB_anova_fam.df, tidy, .id="taxa")[which(map_df(gyrB_anova_fam.df, tidy, .id="taxa")[,"p.value"]<0$
gyrB_anova_gen.sigtab <- map_df(gyrB_anova_gen.df, tidy, .id="taxa")[which(map_df(gyrB_anova_gen.df, tidy, .id="taxa")[,"p.value"]<0$
                                                                         
# Gamma-hurdle (zero and non-zero modelled "indepentently")
#https://seananderson.ca/2014/05/18/gamma-hurdle/
#two models that give 1) m1: th probability of non-zero 
# 2) given that is non-zero the response of y
# 1) binary model
# adjusted to city variable

## BINARY
#ssu
marker_shannon_bin.fam <- marker_tax_stats("marker",idlevel,"family")$shannon_bin
marker_shannon_bin.gen <- marker_tax_stats("marker",idlevel,"genus")$shannon_bin

marker_vars.fam = names(marker_shannon_bin.fam)[!names(marker_shannon_bin.fam) %in% c("bmi","bmi_class","city")]
marker_vars.gen = names(marker_shannon_bin.gen)[!names(marker_shannon_bin.gen) %in% c("bmi","bmi_class","city")]

marker_glmbin_fam.df <- lapply(setNames(marker_vars.fam, marker_vars.fam),function(s) {
  glm(as.formula(paste(as.character(s), " ~ bmi + city")), 
      data = marker_shannon_bin.fam, family = binomial(link = logit))})
ssu_glmbin_gen.df <- lapply(setNames(marker_vars.gen, marker_vars.gen),function(s) {
  glm(as.formula(paste(as.character(s), " ~ bmi + city")), 
      data = marker_shannon_bin.gen, family = binomial(link = logit))})

marker_glmbin_fam.mapdf <- map_df(marker_glmbin_fam.df, tidy, .id="taxa")
marker_glmbin_gen.mapdf <- map_df(marker_glmbin_gen.df, tidy, .id="taxa")

marker_glmbin_fam.sigtab <- marker_glmbin_fam.mapdf %>%
  filter(term=="bmi" & p.value<0.05)
marker_glmbin_gen.sigtab <- ssu_glmbin_gen.mapdf %>%
  filter(term=="bmi" & p.value<0.05)

## CONTINOUS
# first remove columns where all values are 0 (for numeric values)
marker_shannon_cont.fam <- marker_shannon.fam %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)
marker_shannon_cont.gen <- marker_shannon.gen %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)

marker_vars.fam = names(marker_shannon_cont.fam)[!names(marker_shannon_cont.fam) %in% c("bmi","bmi_class","city")] 
marker_vars.gen = names(marker_shannon_cont.gen)[!names(marker_shannon_cont.gen) %in% c("bmi","bmi_class","city")]

marker_glmcont_fam.df <- lapply(setNames(marker_vars.fam, marker_vars.fam),function(s) {
  glm(as.formula(paste(as.character(s), " ~ bmi ")), 
      data = marker_shannon_cont.fam[marker_shannon_cont.fam[which(colnames(marker_shannon_cont.fam)==s)]!=0,], family = Gamma(link = log))})
marker_glmcont_gen.df <- lapply(setNames(marker_vars.gen, marker_vars.gen),function(s) {
  glm(as.formula(paste(as.character(s), " ~ bmi ")), 
      data = marker_shannon_cont.gen[marker_shannon_cont.gen[which(colnames(marker_shannon_cont.gen)==s)]!=0,], family = Gamma(link = log))})

marker_glmcont_fam.mapdf <- map_df(marker_glmcont_fam.df, tidy, .id="taxa")
marker_glmcont_gen.mapdf <- map_df(marker_glmcont_gen.df, tidy, .id="taxa")

marker_glmcont_fam.sigtab <- marker_glmcont_fam.mapdf %>%
  filter(term=="bmi" & p.value<0.055)
marker_glmcont_gen.sigtab <- marker_glmcont_gen.mapdf %>%
  filter(term=="bmi" & p.value<0.055)

#take significant pvalues for bmi (slope) and run again the model with city
newvars.fam = marker_glmcont_fam.sigtab$taxa
newvars.gen = marker_glmcont_gen.sigtab$taxa

marker_glmcont_city_fam.df <- lapply(setNames(newvars.fam, newvars.fam),function(s) {
  glm(as.formula(paste(as.character(s), " ~ bmi + city")), 
      data = marker_shannon_cont.fam[marker_shannon_cont.fam[which(colnames(marker_shannon_cont.fam)==s)]!=0,], family = Gamma(link = log))})
marker_glmcont_city_fam.mapdf <- map_df(marker_glmcont_city_fam.df, tidy, .id="taxa")
marker_glmcont_city_fam.sigtab <- marker_glmcont_city_fam.mapdf %>%
  filter(term=="bmi" & p.value<0.055)

marker_glmcont_city_gen.df <- lapply(setNames(newvars.gen, newvars.gen),function(s) {
  glm(as.formula(paste(as.character(s), " ~ bmi + city")), 
      data = marker_shannon_cont.gen[marker_shannon_cont.gen[which(colnames(marker_shannon_cont.gen)==s)]!=0,], family = Gamma(link = log))})
marker_glmcont_city_gen.mapdf <- map_df(marker_glmcont_city_gen.df, tidy, .id="taxa")
marker_glmcont_city_gen.sigtab <- marker_glmcont_city_gen.mapdf %>%
  filter(term=="bmi" & p.value<0.055)

# Correct for multiple testing for all genera within each family
# subset genera of the family Lachnospiracea
#gyrB
#cont
genera_gyrB <- levels(as.factor(gyrB_TAX_dada2[which(gyrB_TAX_dada2[,"family"]=="f__Lachnospiraceae"),"genus"]))
gyrB_glmcont_gen.genera_gyrB <- gyrB_glmcont_gen.mapdf%>%
  filter(taxa %in% genera_gyrB & term=="bmi")
gyrB_glmcont_gen.genera_gyrB$p.adjust <- p.adjust(gyrB_glmcont_gen.genera_gyrB$p.value,method = "holm")
# 0.4: Ruminococcus_A
#bin
gyrB_glmbin_gen.genera_gyrB <- gyrB_glmbin_gen.mapdf%>%
  filter(taxa %in% genera_gyrB & term=="bmi")
gyrB_glmbin_gen.genera_gyrB$p.adjust <- p.adjust(gyrB_glmbin_gen.genera_gyrB$p.value,method = "holm")

#dnaK
#cont
oscillo_dnak <- levels(as.factor(dnak_taxmat_dada2[which(dnak_taxmat_dada2[,"family"] %in% "f__Oscillospiraceae"),"genus"]))
oscillo_dnak <- str_replace_all(oscillo_dnak, "-","_")
dnak_glmcont_gen.oscillo_dnak <- dnak_glmcont_gen.mapdf%>%
  filter(taxa %in% oscillo_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
#bin
dnak_glmbin_gen.oscillo_dnak <- dnak_glmbin_gen.mapdf%>%
  filter(taxa %in% oscillo_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
#g__CAG_170 (0.0196)
#g__CAG_83 (0.182)
#g__ER4 (0.236)

acuta_dnak <- levels(as.factor(dnak_taxmat_dada2[which(dnak_taxmat_dada2[,"family"] %in% "f__Acutalibacteraceae"),"genus"]))
acuta_dnak <- str_replace_all(acuta_dnak, "-","_")
#cont
dnak_glmcont_gen.acuta_dnak <- dnak_glmcont_gen.mapdf%>%
  filter(taxa %in% acuta_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
# "g__CAG_177" (0.0817)
#bin
dnak_glmbin_gen.acuta_dnak <- dnak_glmbin_gen.mapdf%>%
  filter(taxa %in% acuta_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))

rumino_dnak <- levels(as.factor(dnak_taxmat_dada2[which(dnak_taxmat_dada2[,"family"] %in% "f__Ruminococcaceae"),"genus"]))
rumino_dnak <- str_replace_all(rumino_dnak, "-","_")
#cont
dnak_glmcont_gen.rumino_dnak <- dnak_glmcont_gen.mapdf%>%
  filter(taxa %in% rumino_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
#bin
dnak_glmbin_gen.rumino_dnak <- dnak_glmbin_gen.mapdf%>%
  filter(taxa %in% rumino_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))

#ssu
lachno_ssu <- levels(as.factor(ssu_taxmat_dada2[which(ssu_taxmat_dada2[,"Family"] %in% "Lachnospiraceae"),"Genus"]))
lachno_ssu <- str_replace_all(lachno_ssu, "-","_")
ssu_glmcont_gen.lachno_ssu <- ssu_glmcont_gen.mapdf%>%
  filter(taxa %in% lachno_ssu & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))

cag138_ssu <- levels(as.factor(ssu_taxmat_dada2[which(ssu_taxmat_dada2[,"Family"] %in% "CAG_138"),"Genus"]))
cag138_ssu <- str_replace_all(cag138_ssu, "-","_")
cag138_glmcont_gen.cag138_ssu <- ssu_glmcont_gen.mapdf%>%
  filter(taxa %in% cag138_ssu & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))

cag74_ssu <- levels(as.factor(ssu_taxmat_dada2[which(ssu_taxmat_dada2[,"Family"] %in% "CAG-74"),"Genus"]))
cag74_ssu <- str_replace_all(cag74_ssu, "-","_")
cag74_glmcont_gen.cag74_ssu <- ssu_glmcont_gen.mapdf%>%
  filter(taxa %in% cag74_ssu & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))

#cont
oscillo_ssu <- levels(as.factor(ssu_taxmat_dada2[which(ssu_taxmat_dada2[,"Family"] %in% "Oscillospiraceae"),"Genus"]))
oscillo_ssu <- str_replace_all(oscillo_ssu, "-","_")
ssu_glmcont_gen.oscillo_ssu <- ssu_glmcont_gen.mapdf%>%
  filter(taxa %in% oscillo_ssu & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
#bin
dnak_glmbin_gen.oscillo_dnak <- dnak_glmbin_gen.mapdf%>%
  filter(taxa %in% oscillo_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
#g__CAG_170 (0.0196)
#g__CAG_83 (0.182)
#g__ER4 (0.236)

acuta_dnak <- levels(as.factor(dnak_taxmat_dada2[which(dnak_taxmat_dada2[,"family"] %in% "f__Acutalibacteraceae"),"genus"]))
acuta_dnak <- str_replace_all(acuta_dnak, "-","_")
#cont
dnak_glmcont_gen.acuta_dnak <- dnak_glmcont_gen.mapdf%>%
  filter(taxa %in% acuta_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
# "g__CAG_177" (0.0817)
#bin
dnak_glmbin_gen.acuta_dnak <- dnak_glmbin_gen.mapdf%>%
  filter(taxa %in% acuta_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))

rumino_dnak <- levels(as.factor(dnak_taxmat_dada2[which(dnak_taxmat_dada2[,"family"] %in% "f__Ruminococcaceae"),"genus"]))
rumino_dnak <- str_replace_all(rumino_dnak, "-","_")
#cont
dnak_glmcont_gen.rumino_dnak <- dnak_glmcont_gen.mapdf%>%
  filter(taxa %in% rumino_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
#bin
dnak_glmbin_gen.rumino_dnak <- dnak_glmbin_gen.mapdf%>%
  filter(taxa %in% rumino_dnak & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))

