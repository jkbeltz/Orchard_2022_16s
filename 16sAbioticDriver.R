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
library(plotly)
library(vegan)
library(data.table)
library(janitor)

SVAB <-read_qza("table-allbeltz.qza")
SVnwAB <- read_qza("table-allbeltz_noWR_Wphyla.qza")
OTUTABLE_nw <- SVnwAB$data
ncol(OTUTABLE_nw)
SV_nw <- OTUTABLE_nw #####apply(OTUTABLE_nw, 2, function(x) x/sum(x)*100)

##import files##

denoise_stats_nwAB <-read_qza("denoising-stats-allbeltz.qza")
rep_seqs <- read_qza("rep-seqs-allbeltz.qza")
seqs <- rep_seqs$data
#sickseq=(seqs[2])
as.data.frame(sickseq) ### sequence for weird strain found in 2021

unrooted_artifact <- read_qza("unrooted-tree-allbeltz.qza")
rooted_artifact <- read_qza("rooted-tree-allbeltz.qza")

##taxonomy
taxonomy <- read_qza("taxonomy-allbeltz.qza")
taxon <- parse_taxonomy(taxonomy$data)
#View(taxon)
taxasumgenus <-summarize_taxa(SV_nw, taxon)$Genus
taxasumfamily <-summarize_taxa(SV_nw, taxon)$Family
taxasumorder <-summarize_taxa(SV_nw, taxon)$Order
taxasumclass <-summarize_taxa(SV_nw, taxon)$Class
#View(taxon)

otu.taxon<- merge(taxon, OTUTABLE_nw, by ='row.names', all.y=TRUE)
#View(otu.taxon)

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
metadata$TreatmentInoculum <- gsub('Control Pigmentation Selection', 'ControlPigmentationSelection', metadata$TreatmentInoculum)
metadata$TreatmentInoculum <- gsub('Dark Pigmentation Selection', 'DarkPigmentationSelection', metadata$TreatmentInoculum)
metadata$TreatmentInoculum <- gsub('Light Pigmentation Selection', 'LightPigmentationSelection', metadata$TreatmentInoculum)
metadata <- metadata %>%
  mutate(pop=case_when(str_detect(PopCode,'^E')==T ~ 'Ecage',TRUE ~ 'Indoor'))
metadata <- metadata %>%
  mutate(popslur=case_when(str_detect(PopCode,'^E')==T ~ 'Ecage',TRUE ~ 'Indoor'))
  
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

####heatmaps####

## add a mean column group by treatment
{
m=metadata %>% #### meta file 
  filter (Experiment == '21Pesticide') %>%
  #filter (TreatmentInoculum =='None')%>%
  #filter (CollectionTreatment == 'Starved') %>%
  #filter (Timepoint %in% c('22','22.1','22.2','22.3','22.4')) %>%
  #filter (TreatmentInoculum %in% c('Control','ControlPigmentationSelection','Founder')) %>%
  as_data_frame()#%>% 
  #mutate(pop=case_when(str_detect(PopCode,'^E')==T ~ 'Ecage',TRUE ~ 'Indoor'))

taxa=taxasumgenus  ##change phylo level as needed
taxa<-taxa[order(-rowSums(taxa)),] 
taxa[taxa==0]<-1                            

a =m %>%
  group_by(m$TreatmentInoculum, m$Timepoint, m$SampleID) %>% # add m$Timepoint, 
  dplyr :: summarise(
                     SampleID= SampleID,
                     Timepoint = Timepoint) %>%
                     #SampleType = CollectionTreatment) %>%
                      #Population = pop) %>%
  as_data_frame() %>%
  column_to_rownames('m$SampleID')%>%
  rename(Treatment = 'm$TreatmentInoculum')%>%
  relocate(Treatment, .after = Timepoint)

  #rename(Timepoint = 'm$Timepoint')

select <- a$SampleID

#length(names)
d = taxa%>% #### freq file =
  rownames_to_column('phylo') %>%
  as_tibble() %>%
  select(phylo, any_of(select)) %>% 
  as_data_frame() %>%
  column_to_rownames('phylo') 

a = subset(a, select = -c(SampleID))
a = subset(a, select = -c(`m$Timepoint`))
#a = subset(a, select = -c(Cage))
#View(d)
#View(a)
pheatmap(log10(d[1:10,]), annotation_col = a, cluster_cols = F, cluster_rows = F, labels_col = "") 
}





  ####PCA Viz####

##2D
metadata %>%
  drop_na(PC1_bc)%>%
  filter (Experiment == '22LCD') %>%
  filter (Timepoint %in% c('22.6','22.7','22.1','22.2','22.3')) %>%
  #filter (TreatmentInoculum %in% c('Control','Decreased Larval Density')) %>%
  filter (TreatmentInoculum %in% c('ControlPigmentationSelection')) %>%
  #filter (Experiment == '22LCD') %>%  
  #mutate(pop=case_when(str_detect(PopCode,'^E')==T ~ 'Ecage',TRUE ~ 'Indoor')) %>%
  filter (CollectionTreatment == 'Starved') %>% 
  #filter (!sample_type_detail == 'cage env') %>% 
  #filter (pop == 'Ecage') %>%   
  #filter(TreatmentInoculum == "SummerEvolvedMicrobiome") %>% 
  ggplot(aes(x=PC1_bc, y=PC2_bc, color=Timepoint)) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + #alpha controls transparency and helps when points are overlapping
  theme_q2r()# +
  oh #facet_wrap(~Timepoint)

##3D

colfunc<-colorRampPalette(c("yellow","blue","white","red"))
#View(metadata)
metadata %>%
  drop_na(PC1_bc)%>%
  filter (Experiment == '22NEW') %>%
  #filter (Timepoint %in% c('22','22.1','22.2','22.3','22.4')) %>%
  filter (TreatmentInoculum %in% c('Control','Decreased Larval Density')) %>%
  #filter(TreatmentInoculum == "SummerEvolvedMicrobiome") %>%
  #mutate(pop=case_when(str_detect(PopCode,'^E')==T ~'Ecage',TRUE ~ 'Indoor')) %>%
  filter (CollectionTreatment == 'Starved') %>%
  plot_ly(x = ~PC1_bc, y = ~PC2_bc, z = ~PC3_bc, 
          type = 'scatter3d', 
          mode = 'markers', 
          symbol = ~Timepoint, 
          color = ~TreatmentInoculum, 
          colors = c('red','black'),
          #colors = colfunc(8),
          #symbols = c('o', 'circle'), 
          marker = list(size = 6)) 
View(metadata)
####PERMANOVA####
factor(metadata$Timepoint)

View(metadata)

statdata = metadata %>%
  #filter (Experiment == '22NEW') %>%  
  filter (Timepoint %in% c('22.1','22.2','22.6','22.7')) %>%
  mutate(time_cat = case_when(Timepoint < 22.4 ~ "Early",
                              Timepoint > 22.4 ~ "Late")) %>%
  filter (TreatmentInoculum %in% c('ControlPigmentationSelection')) %>%
  #filter (TreatmentInoculum %in% c('Control','Decreased Larval Density')) %>%
  filter (CollectionTreatment == 'food') %>%
  #mutate(pop=case_when(str_detect(PopCode,'^E')==T ~'Ecage',TRUE ~ 'Indoor'), na.rm=TRUE) %>%
  filter(!is.na(PC1_bc))
#View(statdata)
adonis2(
  statdata[ , c("PC1_bc", "PC2_bc", "PC3_bc")] ~ time_cat,
  data = statdata,
  method = "euc"
)

resSE.aov <- aov(shannon_entropy ~ time_cat, data = statdata)
# Summary of the analysis
summary(resSE.aov)

resOF.aov <- aov(observed_features ~ time_cat, data = statdata)
# Summary of the analysis
summary(resOF.aov)

t.test(shannon_entropy ~ time_cat, data = statdata)


##22NEW - starved - timextreatment ALL SIG for bc & wu 
##22NEW - food - timextreatment  ALL SIG for bc & wu / TP very sig interaction slightly for uu 

##22LCD - starved - timextreatment TP very sig bc, wu, uu / Treatment slightly sig bc 
##22LCD - food - timextreatment TP very sig bc, wu / Interaction slightly sig uu

##21pest - unstarved - timextreatment ALL SIG bc, wu /uu(slightly less), 

##are contr



metadata$interaction <- factor(interaction(metadata$Timepoint, metadata$TreatmentInoculum))
levels(metadata$interaction)
    

####line ploys divesity####

metadata %>%
  drop_na(PC1_bc)%>%
  #filter (Experiment == '22NEW') %>%
  #filter (Timepoint %in% c('22.6','22.7','22.1','22.2','22.3')) %>%
  filter (TreatmentInoculum %in% c('Control','Decreased Larval Density', 'ControlPigmentationSelection')) %>%
  #filter (TreatmentInoculum %in% c('ControlPigmentationSelection')) %>%
  #filter (Experiment == '22SWAP') %>%
  filter (CollectionTreatment == 'Starved') %>%
  #filter (TreatmentInoculum == 'OutdoorMicrobiome') %>%
  ggplot(aes(x=as.factor(Timepoint), y=shannon_entropy, 
             group=interaction(as.factor(TreatmentInoculum), as.factor(Timepoint)), 
             fill=as.factor(TreatmentInoculum))) +
  geom_boxplot() +
  #stat_summary(geom="errorbar", fun.data=mean_se, width=.1) +
  #stat_summary(geom="line", fun.data=mean) 
  #stat_summary(geom="point", fun.y=mean, size=6, shape=12) +
  geom_jitter(shape=19, width=0.15, height=0)+
  #geom_signif(
    # which groups should be compared?
    #comparisons = list(c("Adelie", "Gentoo")), 
    #map_signif_level=TRUE)
  #ylim(0,5)+
  xlab("Experimental Control") +
  #facet_wrap(~Timepoint) +
  #scale_colour_manual(name = "Cage Treatment", 
                      #breaks = c("Control","Pesticide Addition"),
                      #labels = c("Control Populations","LB+ Populations", "AT+ Populations"),
                      #values = c("grey","red"))+
  theme_classic()  # try other themes like theme_bw() or theme_classic()

# Compute the analysis of variance
res.aov <- aov(weight ~ group, data = my_data)
# Summary of the analysis
summary(res.aov)


#####taxa abundance plots

{
taxa=taxasumfamily
taxa <-sweep(taxa, 2, colSums(taxa), `/`) * 100
#view(result)

#select_taxa <- c("Lactobacillaceae","Acetobacteraceae","Enterobacteriaceae","Leuconostocacea", "Enterococcaceae") 
select_taxa <- c("Lactobacillaceae","Acetobacteraceae","Enterobacteriaceae","Leuconostocacea", "Enterococcaceae","Corynebacteriaceae", "Streptococcaceae", "Sphingobacteriaceae", "Weeksellaceae", "Pasteurellales; NA") 
#select_taxa <- c("Lactobacillaceae; NA","Acetobacteraceae; Acetobacter","Enterobacteriaceae; NA","Providencia", "Lactobacillaceae; Lactobacillus","Leuconostocaceae; Weissella", "Enterococcus", "Wautersiella", "Corynebacterium", "Serratia") 
taxa_test<- tibble::rownames_to_column(taxa, "names")
result <- filter(taxa_test, grepl(paste(select_taxa, collapse="|"), names))
result = t(result)
result = result %>%
  row_to_names(row_number = 1)


m =metadata %>%
  drop_na(PC1_bc)%>%
  #filter (Experiment == '22LCD') %>%
  filter (Timepoint %in% c('22','22.1','22.2','22.3')) %>%
  filter (TreatmentInoculum %in% c('Control','ControlPigmentationSelection','Founder')) %>%
  #filter(grepl('Control', TreatmentInoculum))%>%
  #filter (TreatmentInoculum == 'Control')%>%
  #filter (!CollectionTreatment == 'food') %>%
  #filter (Timepoint == '22.5') %>%
  as_data_frame() 
#View(m)

a=m%>% 
  group_by(m$TreatmentInoculum, m$SampleID) %>% # add m$Timepoint, 
  dplyr :: summarise(
                     SampleID= SampleID,
                     Timepoint = Timepoint,
                     Experiment = Experiment) %>%
                     #SampleType = CollectionTreatment) %>%
                     #Population = pop) %>%
  as_data_frame() %>%
  column_to_rownames('m$SampleID')%>%
  rename(Treatment = 'm$TreatmentInoculum')#%>%


combined = merge(a, result, by=0)
#View(result)
combined = combined %>% 
  pivot_longer(
  cols = starts_with("Bacteria"),
  names_to = c("Family"),
  values_to = "RelativeAbundance")

combined %>%
  #filter (Experiment == '22SWAP') %>%
  #filter (CollectionTreatment == 'Starved') %>%
  #filter (!Treatment == 'Founder') %>%
  ggplot(aes(x=as.factor(Timepoint), y=as.numeric(RelativeAbundance),
             #group=interaction(as.factor(Timepoint), as.factor(Family)), 
             fill=as.factor(Treatment))) +
  geom_boxplot()+
  #geom_bar(position="stack", stat="identity")+
  #stat_summary(geom="errorbar", fun.data=mean_se, width=.1) +
  #stat_summary(geom="line", fun.data=mean) 
  #stat_summary(geom="point", fun.y=mean, size=6, shape=12) +
  #geom_jitter(shape=19, width=0.15, height=0)+
  #geom_signif(
  # which groups should be compared?
  #comparisons = list(c("Adelie", "Gentoo")), 
  #map_signif_level=TRUE)
  #ylim(0,5)+
  xlab("Time Point") +
  facet_wrap(~Family) +
  #scale_colour_manual(name = "Cage Treatment", 
  #breaks = c("Control","Pesticide Addition"),
  #labels = c("Control Populations","LB+ Populations", "AT+ Populations"),
  #values = c("grey","red"))+
  theme_classic()  # try other themes like theme_bw() or theme_classic()

#geom_signif(comparisons = list(c("A.c", "B.c"),
#c("A.c", "A.d"),
#c("B.c", "B.d"),
#c("A.d", "B.d")),
# test = "t.test", step_increase = 0.075,
#map_signif_level = TRUE, tip_length = 0)
}


