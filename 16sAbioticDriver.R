#### orchard 16s analysis from multiple experiments, sequencing by Dr. john chason BYU #########################
## started 1/11/24 by jk beltz
## post qiime2 analysis ########
setwd("/Users/jackbeltz/Documents/PENN/Dissertation/CH516sDRIVER/Orchard_2022_16s/QiimeOutputs")

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("tidyverse")
library("devtools")
library(qiime2R)
library(tidyverse)
library(Biostrings)
library(dplyr)
library(dada2)
library(qiime2R)
library(ggplot2)
library(pheatmap)
library(data.table)

SVAB <-read_qza("table-allbeltz.qza")
SVnwAB <- read_qza("table-allbeltz_noWR_Wphyla.qza")
OTUTABLE_nw <- SVnwAB$data
ncol(OTUTABLE_nw)
SV_nw <- OTUTABLE_nw #####apply(OTUTABLE_nw, 2, function(x) x/sum(x)*100)

##import files##

denoise_stats_nwAB <-read_qza("denoising-stats-allbeltz.qza")
rep_seqs <- read_qza("rep-seqs-allbeltz.qza")
seqs <- rep_seqs$data
sickseq=(seqs[2])
as.data.frame(sickseq) ### sequence for weird strain found in 2021

unrooted_artifact <- read_qza("unrooted-tree-allbeltz.qza")
rooted_artifact <- read_qza("rooted-tree-allbeltz.qza")

##taxonomy
taxonomy <- read_qza("taxonomy-allbeltz.qza")
taxon <- parse_taxonomy(taxonomy$data)
View(taxon)
taxasumgenus <-summarize_taxa(SV_nw, taxon)$Genus
taxasumfamily <-summarize_taxa(SV_nw, taxon)$Family
taxasumorder <-summarize_taxa(SV_nw, taxon)$Order
taxasumclass <-summarize_taxa(SV_nw, taxon)$Class
#View(taxon)

otu.taxon<- merge(taxon, OTUTABLE_nw, by ='row.names', all.y=TRUE)
View(otu.taxon)

#### 09a469f13ec1e2269b3b363405f4e8e4 #### the NA 

##metadata
meta <- read.delim("mapper_allbeltz.txt")
meta <- meta %>% rename("SampleID" = "X.SampleID")

##import metricss##
shannon <- read_qza("metrics-noWR-Wphyla-rd100/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")


observed_features <- read_qza("metrics-noWR-Wphyla-rd100/observed_features_vector.qza")
observed_features <- observed_features$data %>% rownames_to_column("SampleID")
#View(observed_features)

faith_pd <- read_qza("metrics-noWR-Wphyla-rd100/faith_pd_vector.qza")
faith_pd <-faith_pd$data %>% rownames_to_column("SampleID")

#view(faith_pd)

rarefied=read_qza("metrics-noWR-Wphyla-rd100/rarefied_table.qza")

## import PCAs ##
uu_distance_matrix <- read_qza("metrics-noWR-Wphyla-rd100/unweighted_unifrac_distance_matrix.qza")
uu_pcoa <- read_qza("metrics-noWR-Wphyla-rd100/unweighted_unifrac_pcoa_results.qza")
#list(wu_pcoa$data$Vectors$SampleID)

wu_distance_matrix <- read_qza("metrics-noWR-Wphyla-rd100/weighted_unifrac_distance_matrix.qza")
wu_pcoa <- read_qza("metrics-noWR-Wphyla-rd100/weighted_unifrac_pcoa_results.qza")

bray_distance_matrix <- read_qza("metrics-noWR-Wphyla-rd100/bray_curtis_distance_matrix.qza")
bray_curtis_pcoa <- read_qza("metrics-noWR-Wphyla-rd100/bray_curtis_pcoa_results.qza")
#view(bray_curtis_pcoa$data$Vectors)
jaccard_matrix <- read_qza("metrics-noWR-Wphyla-rd100/jaccard_distance_matrix.qza")
jaccard_pcoa <- read_qza("metrics-noWR-Wphyla-rd100/jaccard_pcoa_results.qza")

##adding sample metrics to metadata
metadata=meta %>% 
  full_join(observed_features) %>%
  full_join(faith_pd)%>%
  full_join(shannon)
#View(metadata)

metadata=uu_pcoa$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3, PC4) %>%
  full_join(metadata)%>%
  rename("PC1_uu"="PC1", "PC2_uu"="PC2", "PC3_uu"="PC3","PC4_uu"="PC4")

metadata=wu_pcoa$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3, PC4) %>%
  full_join(metadata)%>%
  rename("PC1_wu"="PC1", "PC2_wu"="PC2", "PC3_wu"="PC3", "PC4_wu"="PC4")

metadata=bray_curtis_pcoa$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3, PC4) %>%
  full_join(metadata)%>%
  rename("PC1_bc"="PC1", "PC2_bc"="PC2", "PC3_bc"="PC3", "PC4_bc"="PC4")

metadata$TreatmentInoculum <- gsub('DecreasedLarvalDensity', 'Decreased Larval Density', metadata$TreatmentInoculum)
metadata$TreatmentInoculum <- gsub('OutdoorMicorbiome', 'OutdoorMicrobiome', metadata$TreatmentInoculum)

##splitting metadata by experiment

metadata2 <- split(metadata, metadata$Experiment)
metadata3 <- split(metadata, metadata$TreatmentInoculum)
pest21meta <- metadata2$"21Pesticide"
lcd22meta <- metadata2$"22LCD"
new22meta <- metadata2$"22NEW"
slur22meta <- metadata2$"22SLURRIES"
swap22meta <- metadata2$"22SWAP"
controlmeta <-metadata3$Control


####splitting taxasum by experiment####

##pest21##
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.taxasumgenus_pest21=t.taxasumgenus[pest21meta$SampleID, ]
taxasum_genus_pest21=as.data.frame(t(t.taxasumgenus_pest21))

t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.taxasumfamily_pest21=t.taxasumfamily[pest21meta$SampleID, ]
taxasum_family_pest21=as.data.frame(t(t.taxasumfamily_pest21))

t.taxasumorder=as.data.frame(t(taxasumorder))
t.taxasumorder_pest21=t.taxasumorder[pest21meta$SampleID, ]
taxasum_order_pest21=as.data.frame(t(t.taxasumorder_pest21))

t.taxasumclass=as.data.frame(t(taxasumclass))
t.taxasumclass_pest21=t.taxasumclass[pest21meta$SampleID, ]
taxasum_class_pest21=as.data.frame(t(t.taxasumclass_pest21))


##LCD22#
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.taxasumgenus_lcd22=t.taxasumgenus[lcd22meta$SampleID, ]
taxasum_genus_lcd22=as.data.frame(t(t.taxasumgenus_lcd22))

t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.taxasumfamily_lcd22=t.taxasumfamily[lcd22meta$SampleID, ]
taxasum_family_lcd22=as.data.frame(t(t.taxasumfamily_lcd22))

t.taxasumorder=as.data.frame(t(taxasumorder))
t.taxasumorder_lcd22=t.taxasumorder[lcd22meta$SampleID, ]
taxasum_order_lcd22=as.data.frame(t(t.taxasumorder_lcd22))

t.taxasumclass=as.data.frame(t(taxasumclass))
t.taxasumclass_lcd22=t.taxasumclass[lcd22meta$SampleID, ]
taxasum_class_lcd22=as.data.frame(t(t.taxasumclass_lcd22))


##NEW22#
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.taxasumgenus_new22=t.taxasumgenus[new22meta$SampleID, ]
taxasum_genus_new22=as.data.frame(t(t.taxasumgenus_new22))

t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.taxasumfamily_new22=t.taxasumfamily[new22meta$SampleID, ]
taxasum_family_new22=as.data.frame(t(t.taxasumfamily_new22))

t.taxasumorder=as.data.frame(t(taxasumorder))
t.taxasumorder_new22=t.taxasumorder[new22meta$SampleID, ]
taxasum_order_new22=as.data.frame(t(t.taxasumorder_new22))

t.taxasumclass=as.data.frame(t(taxasumclass))
t.taxasumclass_new22=t.taxasumclass[new22meta$SampleID, ]
taxasum_class_new22=as.data.frame(t(t.taxasumclass_new22))


##Slurry22#
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.taxasumgenus_slur22=t.taxasumgenus[slur22meta$SampleID, ]
taxasum_genus_slur22=as.data.frame(t(t.taxasumgenus_slur22))

t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.taxasumfamily_slur22=t.taxasumfamily[slur22meta$SampleID, ]
taxasum_family_slur22=as.data.frame(t(t.taxasumfamily_slur22))

t.taxasumorder=as.data.frame(t(taxasumorder))
t.taxasumorder_slur22=t.taxasumorder[slur22meta$SampleID, ]
taxasum_order_slur22=as.data.frame(t(t.taxasumorder_slur22))

t.taxasumclass=as.data.frame(t(taxasumclass))
t.taxasumclass_slur22=t.taxasumclass[slur22meta$SampleID, ]
taxasum_class_slur22=as.data.frame(t(t.taxasumclass_slur22))
view(taxasum_genus_slur22)

##Swap22#
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.taxasumgenus_swap22=t.taxasumgenus[swap22meta$SampleID, ]
taxasum_genus_swap22=as.data.frame(t(t.taxasumgenus_swap22))

t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.taxasumfamily_swap22=t.taxasumfamily[swap22meta$SampleID, ]
taxasum_family_swap22=as.data.frame(t(t.taxasumfamily_swap22))

t.taxasumorder=as.data.frame(t(taxasumorder))
t.taxasumorder_swap22=t.taxasumorder[swap22meta$SampleID, ]
taxasum_order_swap22=as.data.frame(t(t.taxasumorder_swap22))

t.taxasumclass=as.data.frame(t(taxasumclass))
t.taxasumclass_swap22=t.taxasumclass[swap22meta$SampleID, ]
taxasum_class_swap22=as.data.frame(t(t.taxasumclass_swap22))


##control#
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.taxasumgenus_control=t.taxasumgenus[controlmeta$SampleID, ]
taxasum_genus_control=as.data.frame(t(t.taxasumgenus_control))

t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.taxasumfamily_control=t.taxasumfamily[controlmeta$SampleID, ]
taxasum_family_control=as.data.frame(t(t.taxasumfamily_control))

t.taxasumorder=as.data.frame(t(taxasumorder))
t.taxasumorder_control=t.taxasumorder[controlmeta$SampleID, ]
taxasum_order_control=as.data.frame(t(t.taxasumorder_control))

t.taxasumclass=as.data.frame(t(taxasumclass))
t.taxasumclass_control=t.taxasumclass[controlmeta$SampleID, ]
taxasum_class_control=as.data.frame(t(t.taxasumclass_control))

####heatmaps###
{
m=metadata %>% #### meta file 
  filter (Experiment == '22NEW') %>%
  #filter (TreatmentInoculum =='FounderMicrobiome')%>%
  filter (CollectionTreatment == 'Starved') %>%
  filter (Timepoint == '22.7') %>%
  as_data_frame()#%>% 
  #mutate(pop=case_when(str_detect(PopCode,'^E')==T ~ 'Ecage',TRUE ~ 'Indoor'))


taxa=taxasumfamily ##change phylo level as needed
taxa<-taxa[order(-rowSums(taxa)),] 
taxa[taxa==0]<-1                            

a =m %>%
  group_by(m$SampleID) %>%
  dplyr :: summarise(SampleID= SampleID,
                     Timepoint = Timepoint,
                     #SampleType = CollectionTreatment,
                     Treatment = TreatmentInoculum) %>%
                    # Population = pop) %>%
  as_data_frame() %>%
  column_to_rownames('m$SampleID')

select <- a$SampleID

#length(names)
d = taxa%>% #### freq file =
  rownames_to_column('phylo') %>%
  as_tibble() %>%
  select(phylo, any_of(select)) %>% 
  as_data_frame() %>%
  column_to_rownames('phylo') 

#d[, select]

a = subset(a, select = -c(SampleID) )

#View(d)
#view(a)
pheatmap(log10(d[1:10,]), annotation_col = a, cluster_cols = F, cluster_rows = F, labels_col = "") 
}


####PCA Viz####
view(metadata)
metadata %>%
  filter (Experiment == '22NEW') %>%  
  #mutate(pop=case_when(str_detect(PopCode,'^E')==T ~ 'Ecage',TRUE ~ 'Indoor')) %>%
  filter (CollectionTreatment == 'Starved') %>% 
  #filter (!sample_type_detail == 'cage env') %>% 
  #filter (Timepoint == '22.5') %>%   
  #filter(TreatmentInoculum == "None") %>% 
  ggplot(aes(x=PC2_wu, y=PC3_wu, color=TreatmentInoculum)) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() # +

  

