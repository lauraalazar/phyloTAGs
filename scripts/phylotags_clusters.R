#source("/lsalazar/proyectos/phylotags/articulo/scripts/phylotags_preprocessing.R")

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

#ssu_shannon.fam <- marker_tax_stats("ssu",85,"family")$shannon
#ssu_shannon.gen <- marker_tax_stats("ssu",85,"genus")$shannon


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
