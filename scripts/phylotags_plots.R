#source("/lsalazar/proyectos/phylotags/articulo/scripts/phylotags_diversity.R")
#source("G:/BackUp_vidarium/BackUp_2023enero/proyectos/phylotags/articulo/scripts/phylotags_preprocessing.R")
#source("G:/BackUp_vidarium/BackUp_2023enero/proyectos/phylotags/articulo/scripts/phylotags_diversity.R")

library(gridExtra)
#################################
# Figure2: Pairwise comparisons #
#################################

#dnak_16s.rds <- readRDS("/lsalazar/proyectos/phylotags/articulo/figures/figure4A_16s_dnak1.rds")
#dnak_gyrB.rds <- readRDS("/lsalazar/proyectos/phylotags/articulo/figures/figure4B_16s_gyrb.rds")
#gyrB_16s.rds <- readRDS("/lsalazar/proyectos/phylotags/articulo/figures/figure4C_dnak1_gyrb.rds")

#comparison_title <- ggdraw() + draw_label("Pairwise comparison among Clostridia genome assemblies", fontface='bold')
#markers.plot <- plot_grid(dnak_16s.rds, gyrB_16s.rds, dnak_gyrB.rds,   ncol=1, nrow=3)
#comparison.plot <- plot_grid(comparison_title,markers.plot, ncol=1,rel_heights=c(0.1, 1))

#svg("/lsalazar/proyectos/phylotags/articulo/figures/figure2_pairwisecomp.svg")
#comparison.plot
#dev.off()

#Supplemetnary Figures of Pairwise comparison for families
#fig5C <- readRDS("/lsalazar/proyectos/phylotags/articulo/msystems/figures/figure5C_16s_dnak1_fam.rds")
#fig6C<- readRDS("/lsalazar/proyectos/phylotags/articulo/msystems/figures/figure6C_16s_gyrb_fam.rds")
#fig7C<- readRDS("/lsalazar/proyectos/phylotags/articulo/msystems/figures/figure7C_dnak1_gyrb_fam.rds")

#suppl_comparison_title <- ggdraw() + draw_label("Pairwise comparison among Clostridia genome assemblies by families", fontface='bold')
#markers.plot <- plot_grid(fig5C, fig6C, fig7C,   ncol=1, nrow=3)
#suppl_comparison.plot <- plot_grid(suppl_comparison_title,markers.plot, ncol=1,rel_heights=c(0.1, 1))

#svg("/lsalazar/proyectos/phylotags/articulo/msystems/figures/figureSuppl_pairwisecomp.svg")
#suppl_comparison.plot
#dev.off()

#################################
# Figure4: GLMs diversity ~ bmi #
#################################

#colour blind friendly palette
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#SIGNIFICANT FAMILIES: continuous model
# only plot the significant families for the markers gyrB and dnaK, with the corresponding 16S results

## dnaK: Oscillospiraceae
oscillo_dnak_fam.df <- cbind(dnak_shannon_cont.fam %>% select(bmi,shannon_cont = "f__Oscillospiraceae"),
                             dnak_shannon_bin.fam %>% select(shannon_bin = "f__Oscillospiraceae"))
oscillo_dnak_fam.df$labelling = ifelse(oscillo_dnak_fam.df$shannon_bin==1, ">0",oscillo_dnak_fam.df$shannon_bin)

oscillo_dnak_pval <- as.numeric(dnak_glmcont_fam.sigtab[dnak_glmcont_fam.sigtab["taxa"]=="f__Oscillospiraceae","p.value"])

# ssu: Oscillospiraceae 
oscillo_ssu_fam.df <- cbind(ssu_shannon_cont.fam %>% select(bmi,shannon_cont = "Oscillospiraceae"),
                            ssu_shannon_bin.fam %>% select(shannon_bin = "Oscillospiraceae"))
oscillo_ssu_fam.df$labelling = ifelse(oscillo_ssu_fam.df$shannon_bin==1, ">0",oscillo_ssu_fam.df$shannon_bin)

oscillo_ssu_pval <- as.numeric(ssu_glmcont_fam.sigtab[ssu_glmcont_fam.sigtab["taxa"]=="Oscillospiraceae","p.value"])

# Oscillospiraceae dnaK + ssu in one plot
oscillo_plot <- plot_grid(ggplot(oscillo_ssu_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                            geom_point(size=2.5) +
                            theme(legend.title=element_blank()) +
                            theme(legend.position = "none")+
                            scale_colour_manual(name="shannon",values=okabe) +
                            ylim(c(0,4)) + 
                            xlim(c(18,40)) + 
                            geom_smooth( data=subset(oscillo_ssu_fam.df,shannon_bin==1),
                                         aes(bmi,shannon_cont),
                                         method = "glm", method.args = list(family = "Gamma"),
                                         se = FALSE,  linewidth=.75) + 
                            annotate("text", x = 25, y = 4, label = paste("cont.pval",round(oscillo_ssu_pval,digits=4),sep="=")) +
                            ggtitle("16S rRNA") +
                            labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"),
                          
                          ggplot(oscillo_dnak_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                            geom_point(size=2.5) +
                            theme(legend.title=element_blank()) +
                            theme(legend.position = "none")+
                            scale_colour_manual(name="shannon",values=okabe) +
                            ylim(c(0,4)) +
                            xlim(c(18,40)) +
                            geom_smooth(data=subset(oscillo_dnak_fam.df,shannon_bin==1),
                                      aes(bmi,shannon_cont),
                                      method = "glm",method.args = list(family = "Gamma"),
                                      se = FALSE,  linewidth=.75) +  
                            annotate("text", x = 25, y = 4,  label = paste("cont.pval",round(oscillo_dnak_pval,digits=4),sep="=")) + 
                            ggtitle("dnaK") + labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"),
                          nrow = 1, ncol = 2, label_size = 10)


oscillo_title <- ggdraw() + draw_label("Oscillospiraceae") #, fontface='bold')

# significant oscillo families
oscillo_fam.plot <- plot_grid(oscillo_title, oscillo_plot,  ncol=1, rel_heights=c(0.1, 1))

#Acutalibacteraceae ssu NO INCLUIR
acuta_ssu_fam.df <- cbind(ssu_shannon_cont.fam %>% select(bmi,shannon_cont = "Acutalibacteraceae"),
                          ssu_shannon_bin.fam %>% select(shannon_bin = "Acutalibacteraceae"))
acuta_ssu_fam.df$labelling = ifelse(acuta_ssu_fam.df$shannon_bin==1, ">0",acuta_ssu_fam.df$shannon_bin)

acuta_ssu_pval <- round(as.numeric(ssu_glmcont_fam.sigtab[ssu_glmcont_fam.sigtab["taxa"]=="Acutalibacteraceae","p.value"]), 
                        digits=3)              

acuta_plot <- ggplot(acuta_ssu_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) +
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")+
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) + 
  xlim(c(18,40)) + 
  geom_smooth(data=subset(acuta_ssu_fam.df,shannon_bin==1),
              aes(bmi,shannon_cont),
              method = "glm",method.args = list(family = "Gamma"),
              se = FALSE,  linewidth=.75) +  
  annotate("text", x = 25, y = 3.5, label = paste("pval",acuta_ssu_pval,sep="=")) +
  ggtitle("16S rRNA") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n")

acuta_title <- ggdraw() + draw_label("Acutalibacteraceae") #, fontface='bold')

# significant acuta family
acuta_fam.plot <- plot_grid(acuta_title, acuta_plot,  ncol=1, rel_heights=c(0.1, 1))
#plot_grid(acuta_title, acuta_plot,  ncol=1, rel_heights=c(0.1, 1))

#Ruminococceae ssu NO INCLUIR
rumino_ssu_fam.df <- cbind(ssu_shannon_cont.fam %>% select(bmi,shannon_cont = "Ruminococcaceae"),
                          ssu_shannon_bin.fam %>% select(shannon_bin = "Ruminococcaceae"))
rumino_ssu_fam.df$labelling = ifelse(rumino_ssu_fam.df$shannon_bin==1, ">0",rumino_ssu_fam.df$shannon_bin)

rumino_ssu_pval <- round(as.numeric(ssu_glmcont_fam.sigtab[ssu_glmcont_fam.sigtab["taxa"]=="Ruminococcaceae","p.value"]), 
                        digits=3)              

rumino_plot <- ggplot(rumino_ssu_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) +
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")+
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) + 
  xlim(c(18,40)) + 
  geom_smooth(data=subset(rumino_ssu_fam.df,shannon_bin==1),
              aes(bmi,shannon_cont),
              method = "glm",method.args = list(family = "Gamma"),
              se = FALSE,  linewidth=.75) +  
  annotate("text", x = 25, y = 3.5, label = paste("pval",rumino_ssu_pval,sep="=")) +
  ggtitle("16S rRNA") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n")

rumino_title <- ggdraw() + draw_label("Ruminococcaceae") #, fontface='bold')

# significant rumino family
rumino_fam.plot <- plot_grid(rumino_title, rumino_plot,  ncol=1, rel_heights=c(0.1, 1))

#CAG-138 ssu NO INCLUIR
cag138_ssu_fam.df <- cbind(ssu_shannon_cont.fam %>% select(bmi,shannon_cont = "CAG_138"),
                           ssu_shannon_bin.fam %>% select(shannon_bin = "CAG_138"))
cag138_ssu_fam.df$labelling = ifelse(cag138_ssu_fam.df$shannon_bin==1, ">0",cag138_ssu_fam.df$shannon_bin)

cag138_ssu_pval <- round(as.numeric(ssu_glmcont_fam.sigtab[ssu_glmcont_fam.sigtab["taxa"]=="CAG_138","p.value"]), 
                         digits=3)              

cag138_plot <- ggplot(cag138_ssu_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) +
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")+
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) + 
  xlim(c(18,40)) + 
  geom_smooth(data=subset(cag138_ssu_fam.df,shannon_bin==1),
              aes(bmi,shannon_cont),
              method = "glm",method.args = list(family = "Gamma"),
              se = FALSE,  linewidth=.75) +  
  annotate("text", x = 25, y = 3.5, label = paste("pval",cag138_ssu_pval,sep="=")) +
  ggtitle("16S rRNA") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n")

cag138_title <- ggdraw() + draw_label("CAG_138") #, fontface='bold')

# significant rumino family
cag138_fam.plot <- plot_grid(cag138_title, cag138_plot,  ncol=1, rel_heights=c(0.1, 1))



# family plot
#family_title <- ggdraw() + draw_label("Siginificant families", fontface='bold')
#sigfam.plot <- plot_grid(gyrB_fam.plot, oscillo_fam.plot,ncol=1, nrow=2)



## gyrB: Lachnospiraceae
#Lachno gyrB 
lachno_gyrB_fam.df <- cbind(gyrB_shannon_cont.fam %>% select(bmi,shannon_cont = "f__Lachnospiraceae"),
                            gyrB_shannon_bin.fam %>% select(shannon_bin = "f__Lachnospiraceae"))
lachno_gyrB_fam.df$labelling = ifelse(lachno_gyrB_fam.df$shannon_bin==1, ">0",lachno_gyrB_fam.df$shannon_bin)
lachno_gyrB_pval <- as.numeric(gyrB_glmcont_fam.sigtab[gyrB_glmcont_fam.sigtab["taxa"]=="f__Lachnospiraceae","p.value"])

#Lachno ssu
lachno_ssu_fam.df <- cbind(ssu_shannon_cont.fam %>% select(bmi,shannon_cont = "Lachnospiraceae"),
                           ssu_shannon_bin.fam %>% select(shannon_bin = "Lachnospiraceae"))
lachno_ssu_fam.df$labelling = ifelse(lachno_ssu_fam.df$shannon_bin==1, ">0",lachno_ssu_fam.df$shannon_bin)
lachno_ssu_pval <- as.numeric(ssu_glmcont_fam.sigtab[ssu_glmcont_fam.sigtab["taxa"]=="Lachnospiraceae","p.value"])#0.00365 #gyrB_glmcont_city_fam.sigtab$p.value

lachno_title <- ggdraw() + draw_label("Lachnospiraceae") #, fontface='bold')
lachno_plot <- plot_grid(ggplot(lachno_ssu_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                           geom_point(size=2.5) +
                           theme(legend.title=element_blank()) +
                           theme(legend.position = "none") +
                           scale_colour_manual(name="shannon",values=okabe) +
                           ylim(c(0,4)) + 
                           xlim(c(18,40)) + 
                           geom_smooth(data=subset(lachno_ssu_fam.df,shannon_bin==1), 
                                     aes(bmi,shannon_cont), 
                                     method = "glm", method.args = list(family = "Gamma"),
                                     se = FALSE, linewidth =.75) + 
                           annotate("text", x = 25, y = 3.9, label = paste("cont.pval",round(lachno_ssu_pval,digits = 4),sep="=")) +  #"0.023"
                           ggtitle("16S rRNA") +
                           labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"),
                         ggplot(lachno_gyrB_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                           geom_point(size=2.5) +
                           theme(legend.title=element_blank()) +
                           theme(legend.position = "none")+
                           scale_colour_manual(name="shannon",values=okabe) +
                           ylim(c(0,4)) +
                           xlim(c(18,40)) +
                           geom_smooth(data=subset(lachno_gyrB_fam.df,shannon_bin==1),
                                       aes(bmi,shannon_cont), 
                                       method = "glm", method.args = list(family = "Gamma"),
                                       se = FALSE, linewidth =.75) + 
                           annotate("text", x = 25, y = 3.9,  label = paste("cont.pval",round(lachno_gyrB_pval,digits = 4),sep="=")) + 
                           ggtitle("gyrB") + labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"), 
                         nrow = 1, ncol = 2, label_size = 10)

# significant gyrB families
gyrB_fam.plot <- plot_grid(lachno_title, lachno_plot,  ncol=1, rel_heights=c(0.1, 1))


#All families together in panel A.
#sigfam_title <- cowplot::ggdraw() + cowplot::draw_label("A. Significant Families", size = 20,hjust = 1.5)
#family_plot <- cowplot::plot_grid(sigfam_title,gyrB_fam.plot, oscillo_fam.plot, ncol=1, rel_heights=c(0.1,0.4,0.4))

#pdf("/lsalazar/proyectos/phylotags/articulo/msystems/figures/Fig4A_sig_fam.pdf",height=10,width=10)
#family_plot
#dev.off()


#Genera
#Oscillospiraceae
#CAG-83
cag83 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_83") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag83_cont.plot <-ggplot(cag83)+
  geom_point(size=2.5,aes(bmi, shannon_cont , color= labelling)) +
#  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) + 
  xlim(c(18,40)) + 
  geom_smooth(data = cag83, aes(x = bmi, y = shannon_cont),
              method = "lm", #method.args = list(family = "binomial"), 
              se = FALSE, fullrange=TRUE, color = okabe[1]) + 
 # annotate("text", x = 30, y = 0.75, label = "bin.cont=NS") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
 # ggtitle("dnaK: CAG-83 (Oscillospiraceae)") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon \n")

cag83.pval <- as.numeric(dnak_glmbin_gen.sigtab[dnak_glmbin_gen.sigtab["taxa"]=="g__CAG_83","p.value"])
cag83.fdr <- as.numeric(dnak_glmbin_gen.oscillo_dnak[dnak_glmbin_gen.oscillo_dnak["taxa"]=="g__CAG_83","padj"])

cag83_bin.boxplot <-ggplot(cag83)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = 0, y = 20,  label = paste("bin.pval",round(cag83.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag83.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("CAG-83 (Oscillospiraceae) dnaK")


cag83_bin.violin <-ggplot(cag83)+
  geom_dotplot(size=2,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = 39,  label = paste("bin.pval",round(cag83.pval,digits = 4),sep="=")) +
  annotate("text", x = 0.5, y = 37,  label = paste("bin.fdr",round(cag83.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")"))


cag83_model<-plot_grid(cag83_cont.plot,cag83_bin.boxplot,ncol=2,nrow=1)

#cag-83 full plot
cag83_title <- ggdraw() + draw_label("dnaK: CAG-83 (Oscillospiraceae)") #, fontface='bold')
cag83.plot <- plot_grid(cag83_title, cag83_model, nrow=2, ncol=1, rel_heights=c(0.1, 1))

#CAG-170
cag170 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_170") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag170_cont.plot <-ggplot(cag170)+
  geom_point(size=2.5,aes(bmi, shannon_cont , color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) + 
  xlim(c(18,40)) + 
  geom_smooth(data = cag170, aes(x = bmi, y = shannon_cont),
              method = "lm", #method.args = list(family = "binomial"), 
              se = FALSE, fullrange=TRUE, color = okabe[1]) + 
  # annotate("text", x = 30, y = 0.75, label = "bin.cont=NS") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  # ggtitle("dnaK: CAG-170 (Oscillospiraceae)") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon \n")

cag170.pval <- as.numeric(dnak_glmbin_gen.sigtab[dnak_glmbin_gen.sigtab["taxa"]=="g__CAG_170","p.value"])
cag170.fdr <- as.numeric(dnak_glmbin_gen.oscillo_dnak[dnak_glmbin_gen.oscillo_dnak["taxa"]=="g__CAG_170","padj"])

cag170_bin.boxplot <-ggplot(cag170)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = 0, y = 20,  label = paste("bin.pval",round(cag170.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag170.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")"))+
  ggtitle("CAG-170 (Oscillospiraceae) dnaK")

cag170_model<-plot_grid(cag170_cont.plot,cag170_bin.boxplot,ncol=2,nrow=1)

#cag-170 full plot
cag170_title <- ggdraw() + draw_label("dnaK: CAG-170 (Oscillospiraceae)") #, fontface='bold')
#cag170.plot <- plot_grid(cag170_title, cag170_model, nrow=2, ncol=1, rel_heights=c(0.1, 1))

#CAG-177
cag177 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_177") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag177.pval <- as.numeric(dnak_glmcont_gen.sigtab[dnak_glmcont_gen.sigtab["taxa"]=="g__CAG_177","p.value"])
cag177.fdr <- as.numeric(dnak_glmcont_gen.acuta_dnak[dnak_glmcont_gen.acuta_dnak["taxa"]=="g__CAG_177","padj"])

cag177_cont.plot <-ggplot(cag177,  aes(bmi,shannon_cont,color=as.factor(labelling) ))+
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")+
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) +
  xlim(c(18,40)) +
  geom_smooth(data=subset(cag177,shannon_bin==1),
              aes(bmi,shannon_cont), 
              method = "glm", method.args = list(family = "Gamma"),
              se = FALSE, linewidth =.75) + 
  annotate("text", x = 25, y = 3.9,  label = paste("cont.pval",round(cag177.pval,digits = 4),sep="=")) +
  annotate("text", x = 25, y = 3.5, label = paste("cont.fdr",round(cag177.fdr,digits = 4),sep="=")) +
  ggtitle("CAG-177 (Acutalibacteraceae) dnaK") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  # ggtitle("dnaK: CAG-177 (Acutalibacteraceae)") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon \n")


cag177_bin.boxplot <-ggplot(cag177)+
  geom_boxplot(size=2,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = 39,  label = paste("bin.pval","0.015",sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("CAG-177 (Acutalibacteraceae) dnaK")

cag177_model<-plot_grid(cag177_cont.plot,cag177_bin.boxplot,ncol=2,nrow=1)

#All genera together in panel B.
siggen_title <- cowplot::ggdraw() + cowplot::draw_label("B. Significant Genera", size = 20,hjust = 1.5)
oscillo_gen <- plot_grid(cag83_bin.boxplot,cag170_bin.boxplot, ncol=2, nrow = 1, rel_heights=c(0.5,0.5))
genera_plot <- cowplot::plot_grid(siggen_title,oscillo_gen, ncol=2, nrow = 2, rel_heights=c(0.1,0.9))

#CAG-488
cag488 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_488") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag488.pval <- as.numeric(dnak_glmbin_gen.sigtab[dnak_glmbin_gen.sigtab["taxa"]=="g__CAG_488","p.value"])
cag488.fdr <- as.numeric(dnak_glmbin_gen.acuta_dnak[dnak_glmbin_gen.acuta_dnak["taxa"]=="g__CAG_488","padj"])

cag488_bin.boxplot <-ggplot(cag488)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = 0, y = 20,  label = paste("bin.pval",round(cag488.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag488.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("CAG-488 (Acutalibacteraceae) dnaK")

#ER4
er4 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__ER4") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

er4.pval <- as.numeric(dnak_glmbin_gen.sigtab[dnak_glmbin_gen.sigtab["taxa"]=="g__ER4","p.value"])
er4.fdr <- as.numeric(dnak_glmbin_gen.oscillo_dnak[dnak_glmbin_gen.oscillo_dnak["taxa"]=="g__ER4","padj"])

er4_bin.boxplot <-ggplot(er4)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = 0, y = 20,  label = paste("bin.pval",round(er4.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(er4.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("ER4 (Oscillospiraceae) dnaK")


#Fig4 <- plot_grid(family_plot, genera_plot, nrow=1, ncol=2, rel_widths = c(2,1))


#pdf("/lsalazar/proyectos/phylotags/articulo/msystems/figures/Fig4B_sig_gen.pdf",height=10,width=10)
#genera_plot
#dev.off()

#svg("/lsalazar/proyectos/phylotags/articulo/msystems/figures/sig_fam_gen.svg",height=10,width=10)
#Fig4
#dev.off()


## 16S Genera Bin
#Acetatifactor
aceta <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "Acetatifactor") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

aceta.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="Acetatifactor","p.value"])
aceta.fdr <- as.numeric(ssu_glmbin_gen.lachno_ssu[ssu_glmbin_gen.lachno_ssu["taxa"]=="Acetatifactor","padj"])

aceta_bin.boxplot <-ggplot(aceta)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(aceta.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(aceta.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("Acetatifactor (Lachnospiraceae) 16S")

# CAG-41
cag41 <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "CAG_41") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag41.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="CAG_41","p.value"])
cag41.fdr <- as.numeric(UBA1381_glmbin_gen.UBA1381_ssu[UBA1381_glmbin_gen.UBA1381_ssu["taxa"]=="CAG_41","padj"])

cag41_bin.boxplot <-ggplot(cag41)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(cag41.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag41.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("CAG_41 (UBA1381) 16S")

# CAG-83
cag83 <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "CAG_83") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag83.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="CAG_83","p.value"])
cag83.fdr <- as.numeric(ssu_glmbin_gen.oscillo_ssu[ssu_glmbin_gen.oscillo_ssu["taxa"]=="CAG_83","padj"])

cag83_bin.boxplot <-ggplot(cag83)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(cag83.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag83.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("CAG_83 (Oscillospiraceae) 16S")

# Dysosmobacter
dysosmo <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "Dysosmobacter") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

dysosmo.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="Dysosmobacter","p.value"])
dysosmo.fdr <- as.numeric(ssu_glmbin_gen.oscillo_ssu[ssu_glmbin_gen.oscillo_ssu["taxa"]=="Dysosmobacter","padj"])

dysosmo_bin.boxplot <-ggplot(dysosmo)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(cag83.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag83.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("Dysosmobacter (Oscillospiraceae) 16S")

# ER4
er4ssu <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "ER4") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

er4ssu.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="ER4","p.value"])
er4ssu.fdr <- as.numeric(ssu_glmbin_gen.oscillo_ssu[ssu_glmbin_gen.oscillo_ssu["taxa"]=="ER4","padj"])

er4ssu_bin.boxplot <-ggplot(er4ssu)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(er4ssu.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(er4ssu.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("ER4 (Oscillospiraceae) 16S")

# OEMS01
oem <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "OEMS01") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

oem.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="OEMS01","p.value"])
oem.fdr <- as.numeric(cag74_glmbin_gen.cag74_ssu[cag74_glmbin_gen.cag74_ssu["taxa"]=="OEMS01","padj"])

oem_bin.boxplot <-ggplot(oem)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(oem.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(oem.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("OEMS01 (CAG-74) 16S")

# UBA11524
uba11524 <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "UBA11524") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

uba11524.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="UBA11524","p.value"])
uba11524.fdr <- as.numeric(cag74_glmbin_gen.cag74_ssu[cag74_glmbin_gen.cag74_ssu["taxa"]=="UBA11524","padj"])

uba11524_bin.boxplot <-ggplot(oem)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(uba11524.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(uba11524.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("UBA11524 (CAG-74) 16S")

# PHIL1
phil1 <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "Phil1") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

phil1.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="Phil1","p.value"])
phil1.fdr <- as.numeric(cag138_glmbin_gen.cag138_ssu[cag138_glmbin_gen.cag138_ssu["taxa"]=="Phil1","padj"])

phil1_bin.boxplot <-ggplot(phil1)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(phil1.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(phil1.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("Phil1 (CAG-138) 16S")


# UBA3792
uba3792 <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "UBA3792") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

uba3792.pval <- as.numeric(ssu_glmbin_gen.sigtab[ssu_glmbin_gen.sigtab["taxa"]=="UBA3792","p.value"])
uba3792.fdr <- as.numeric(cag138_glmbin_gen.cag138_ssu[cag138_glmbin_gen.cag138_ssu["taxa"]=="UBA3792","padj"])

uba3792_bin.boxplot <-ggplot(uba3792)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = -0.5, y = 20,  label = paste("bin.pval",round(uba3792.pval,digits = 4),sep="=")) +
  annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(uba3792.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("UBA3792 (CAG-138) 16S")


###################################
# Figure5: Differential Abundance #
###################################
source("/lsalazar/proyectos/phylotags/articulo/scripts/phylotags_abundance.R")


#order ASVs according to foldChange 
dnak_gathered_sigASV_nw_ob$asv <- factor(dnak_gathered_sigASV_nw_ob$asv, 
                                         levels=rownames(dnak_sigtab.nw_ob[order(dnak_sigtab.nw_ob$log2FoldChange),]))
gyrB_gathered_sigASV_nw_ob$asv <- factor(gyrB_gathered_sigASV_nw_ob$asv, 
                                         levels=rownames(gyrB_sigtab.nw_ob[order(gyrB_sigtab.nw_ob$log2FoldChange),]))

#change labels of ASVs codes for species 
plot_labels_dnak = sapply(strsplit(as.character(dnak_sigtab.nw_ob[rownames(dnak_sigtab.nw_ob[order(dnak_sigtab.nw_ob$log2FoldChange),]),]$species), 
                                   "\\__"), `[`, 2)

plot_labels_gyrB = sapply(strsplit(as.character(gyrB_sigtab.nw_ob[rownames(gyrB_sigtab.nw_ob[order(gyrB_sigtab.nw_ob$log2FoldChange),]),]$species),
                                   "\\__"), `[`, 2)

#pdf("/lsalazar/proyectos/phylotags/articulo/msystems/figures/diffab_dnak.pdf",height=5.5,width=10)
diffab_dnak <- ggplot(dnak_gathered_sigASV_nw_ob, aes(asv, log10(normalized_counts +1))) +
  geom_jitter(aes(color=bmi_class), position=position_jitter(0.2), size = 2) +
  geom_point(aes(group=bmi_class, colour=bmi_class), stat='summary', fun.y='mean', fun.args = list(mult=1), geom = "pointrange",  size = 0.4) + 
  #  geom_label()+
  geom_line(aes(group=bmi_class, colour=bmi_class), stat='summary', fun.y='mean' , size = 2) +
  scale_color_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07"))+
  #  scale_y_log10() +
  ylab("Log10 Normalized counts") +
  xlab("") +
  scale_x_discrete(labels = plot_labels_dnak) +
  ggtitle("dnaK Siginificant ASVs ") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), 
        plot.margin = margin(1.5,1.5,1.5,1.5,"cm")) +
  guides(color=guide_legend(title="BMI class"))
#dev.off()

#pdf("/lsalazar/proyectos/phylotags/articulo/msystems/figures/diffab_gyrB.pdf",height=5.5,width=10)
diffab_gyrB <- ggplot(gyrB_gathered_sigASV_nw_ob, aes(asv, log10(normalized_counts +1))) +  geom_jitter(aes(color=bmi_class), position=position_jitter(0.2), size = 2) +
  geom_point(aes(group=bmi_class, colour=bmi_class), stat='summary', fun.y='mean', fun.args = list(mult=1), geom = "pointrange",  size = 0.4) +
  geom_line(aes(group=bmi_class, colour=bmi_class), stat='summary', fun.y='mean' , size = 2) +
  scale_color_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07"))+
  #  scale_y_log10() +
  ylab("Log10 Normalized counts") +
  xlab("") +
  scale_x_discrete(labels = plot_labels_gyrB) +
  ggtitle("gyrB Siginificant ASVs ") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), 
        plot.margin = margin(1.5,1.5,1.5,1.5,"cm")) +
  guides(color=guide_legend(title="BMI class"))
#dev.off()


er4 <- dnak_shannon.gen %>% 
  select(bmi,shannon_cont = "g__ER4") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

er4.pval <- as.numeric(dnak_glmbin_gen.sigtab[dnak_glmbin_gen.sigtab["taxa"]=="g__ER4","p.value"])
er4.fdr <- as.numeric(dnak_glmbin_gen.oscillo_dnak[dnak_glmbin_gen.oscillo_dnak["taxa"]=="g__ER4","padj"])

er4_bin.boxplot <-ggplot(er4)+
  geom_boxplot(size=2,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = 39,  label = paste("bin.pval",round(er4.pval,digits = 4),sep="=")) +
  annotate("text", x = 0.5, y = 37,  label = paste("bin.fdr",round(er4.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("ER4 (Oscillospiraceae) dnaK")


#CAG-177
cag177ssu <- ssu_shannon.gen %>% 
  select(bmi,shannon_cont = "CAG_177") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag177ssu.pval <- as.numeric(ssu_glmcont_gen.sigtab[ssu_glmcont_gen.sigtab["taxa"]=="CAG_177","p.value"])
cag177ssu.fdr <- as.numeric(acuta_glmcont_gen.acuta_ssu[acuta_glmcont_gen.acuta_ssu["taxa"]=="CAG_177","padj"])

cag177ssu_cont.plot <-ggplot(cag177ssu,  aes(bmi,shannon_cont,color=as.factor(labelling) ))+
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")+
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) +
  xlim(c(18,40)) +
  geom_smooth(data=subset(cag177ssu,shannon_bin==1),
              aes(bmi,shannon_cont), 
              method = "glm", method.args = list(family = "Gamma"),
              se = FALSE, linewidth =.75) + 
  annotate("text", x = 25, y = 3.9,  label = paste("cont.pval",round(cag177ssu.pval,digits = 4),sep="=")) +
  annotate("text", x = 25, y = 3.5, label = paste("cont.fdr",round(cag177ssu.fdr,digits = 4),sep="=")) +
  ggtitle("CAG-177 (Acutalibacteraceae) 16S") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon \n")


#Hungatella
hung <- gyrB_shannon.gen %>% 
  select(bmi,shannon_cont = "g__Hungatella") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

hung.pval <- as.numeric(gyrB_glmcont_gen.sigtab[gyrB_glmcont_gen.sigtab["taxa"]=="g__Hungatella","p.value"])
hung.fdr <- as.numeric(gyrB_glmcont_gen.genera_gyrB[gyrB_glmcont_gen.genera_gyrB["taxa"]=="g__Hungatella","p.adjust"])

hung_cont.plot <-ggplot(hung,  aes(bmi,shannon_cont,color=as.factor(labelling) ))+
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")+
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) +
  xlim(c(18,40)) +
  geom_smooth(data=subset(cag177ssu,shannon_bin==1),
              aes(bmi,shannon_cont), 
              method = "glm", method.args = list(family = "Gamma"),
              se = FALSE, linewidth =.75) + 
  annotate("text", x = 25, y = 3.9,  label = paste("cont.pval",round(hung.pval,digits = 4),sep="=")) +
  annotate("text", x = 25, y = 3.5, label = paste("cont.fdr",round(hung.fdr,digits = 4),sep="=")) +
  ggtitle("Hungatella (Lachnospiraceae) 16S") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon \n")

#Ruminococcus_A
ruminoA <- gyrB_shannon.gen %>% 
  select(bmi,shannon_cont = "g__Ruminococcus_A") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

ruminoA.pval <- as.numeric(gyrB_glmcont_gen.sigtab[gyrB_glmcont_gen.sigtab["taxa"]=="g__Ruminococcus_A","p.value"])
ruminoA.fdr <- as.numeric(gyrB_glmcont_gen.genera_gyrB[gyrB_glmcont_gen.genera_gyrB["taxa"]=="g__Ruminococcus_A","p.adjust"])

ruminoA_cont.plot <-ggplot(ruminoA,aes(bmi,shannon_cont,color=as.factor(labelling) ))+
  geom_point(size=2.5) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none")+
  scale_colour_manual(name="shannon",values=okabe) +
  ylim(c(0,4)) +
  xlim(c(18,40)) +
  geom_smooth(data=subset(ruminoA,shannon_bin==1),
              aes(bmi,shannon_cont), 
              method = "glm", method.args = list(family = "Gamma"),
              se = FALSE, linewidth =.75) + 
  annotate("text", x = 25, y = 3.9,  label = paste("cont.pval",round(ruminoA.pval,digits = 4),sep="=")) +
  annotate("text", x = 25, y = 3.5, label = paste("cont.fdr",round(ruminoA.fdr,digits = 4),sep="=")) +
  ggtitle("Ruminococcus_A (Lachnospiraceae) 16S") +
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon \n")
