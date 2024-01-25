#### orchard 2020 16s analysis #########################
## started 9/14/21 by jk beltz
## post qiime2 analysis ########
install.packages("tidyverse")
library("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

BiocManager::install("dada2", version = "3.17")

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.17")
packageVersion("dada2")
#!/usr/bin/env Rscript
library(qiime2R)
library(tidyverse)
library(Biostrings)
library(dplyr)
install.packages("dada2")
library(dada2)
library(qiime2R)
library(ggplot2)

### online tutorial https://www.google.com/search?q=qiimer+vs.+qiime2r&oq=qiimer+vs.+qiime2r&aqs=chrome..69i57j0i22i30.11627j0j7&sourceid=chrome&ie=UTF-8
##Reading QZA files#######
##SET WD to HOME/DOCUMENTS/16S_QIIME2
setwd("/Users/jackbeltz/Desktop/16s Drivers/Chaston Files")

meta=beltz_mapper_combined
View(meta)
meta1<-read_tsv("mapper-beltz1to4b.tsv")
meta2<-read_tsv('mapper_beltzrun2.tsv')
meta3<-read.delim("mapper_beltzrun3.txt")
View(meta2)

table1<-read_qza("table-beltz1to4_noWR.qza")
table2<-read_qza("table-beltzrun2_noWR.qza")
table3<-read_qza("table-beltzrun3_noWR.qza")

OTUTABLE1=table1$data
SV1<-apply(OTUTABLE1, 2, function(x) x/sum(x)*100)
OTUTABLE2=table2$data
SV2<-apply(OTUTABLE2, 2, function(x) x/sum(x)*100)
OTUTABLE3=table3$data
SV3<-apply(OTUTABLE3, 2, function(x) x/sum(x)*100)

taxonomy1<-read_qza("taxonomy-beltz1to4.qza")
taxonomy2<-read_qza("taxonomy-beltzrun2.qza")
taxonomy3<-read_qza("taxonomy-beltzrun3.qza")

taxon1<-parse_taxonomy(taxonomy1$data)
taxon2<-parse_taxonomy(taxonomy2$data)
taxon3<-parse_taxonomy(taxonomy3$data)

taxa1sumgenus<-summarize_taxa(SV1, taxon1)$Genus
taxa1sumfamily<-summarize_taxa(SV1, taxon1)$Family
taxa1sumorder<-summarize_taxa(SV1, taxon1)$Order
taxa1sumclass<-summarize_taxa(SV1, taxon1)$Class
taxa1sumking<-summarize_taxa(SV1, taxon1)$Kingdom
View(taxa1sumgenus)

taxa2sumgenus<-summarize_taxa(SV2, taxon2)$Genus
taxa2sumfamily<-summarize_taxa(SV2, taxon2)$Family
taxa2sumorder<-summarize_taxa(SV2, taxon2)$Order
taxa2sumclass<-summarize_taxa(SV2, taxon2)$Class
taxa2sumking<-summarize_taxa(SV2, taxon2)$Kingdom
View(taxa2sumgenus)

otu.taxon1<- merge(taxon1, OTUTABLE1, by ='row.names', all=TRUE)
otu.taxon2<- merge(taxon2, OTUTABLE2, by ='row.names', all=TRUE)
otu.taxon3<- merge(taxon3, OTUTABLE3, by ='row.names', all=TRUE)


SLURRYmeta = subset(meta2, meta2$Experiment =="22SLURRIES")


t.taxa2sumgenus=as.data.frame(t(taxa2sumgenus))
t.slurrytaxasumgenus=t.taxa2sumgenus[SLURRYmeta$`#SampleID`, ]
slurrytaxasumgenus.long=as.data.frame(t.slurrytaxasumgenus)
View(slurrytaxasumgenus.long)
slurrytaxasumgenus=as.data.frame(t(t.slurrytaxasumgenus))
View(slurrytaxasumgenus)

t.taxa2sumfamily=as.data.frame(t(taxa2sumfamily))
t.slurrytaxasumfamily=t.taxa2sumfamily[SLURRYmeta$`#SampleID`, ]
slurrytaxasumfamily=as.data.frame(t(t.slurrytaxasumfamily))
t.taxa2sumorder=as.data.frame(t(taxa2sumorder))
t.slurrytaxasumorder=t.taxa2sumorder[SLURRYmeta$`#SampleID`, ]
slurrytaxasumorder=as.data.frame(t(t.slurrytaxasumorder))
View(slurrytaxasumgenus)

SLURRYmeta$TreatmentInoculum = as.factor(SLURRYmeta$TreatmentInoculum)

slurryall=taxa_heatmap(slurrytaxasumfamily, SLURRYmeta, normalize=TRUE, ntoplot=12)
slurryall + facet_wrap( ~ SLURRYmeta$TreatmentInoculum, scales = "free_x") 

rownames(SLURRYmeta) <- SLURRYmeta$`#SampleID`

slurrygenus<-merge(slurrytaxasumgenus.long, by='row.names')
View(slurrygenus)

mine.long <- pivot_longer(data = slurrytaxasumgenus.long,
                          names_to = "Genus", 
                          values_to = "Abundance")
head(mine.long)


SWAPmeta = subset(meta2, meta2$Experiment =="22SWAP")
View(SWAPmeta)



library(pheatmap)

#load in the OTU or function table (here I am using real counts not relative abundances)
d<-slurrytaxasumgenus

#order the table by descending number of observations (OTUs or functions(e.g. KOs)), 
#so we can just plot the top most number or rows later
d_sorted<-d[order(-rowSums(d)),]

#we are going to log to increase visibility of lower counts so need to get rid of these 0's 
#This is obviously not ok if you have small counts, 
#but on a big heatmap with large counts the difference between 0 or 1 will not be distiguishable
d_sorted[d_sorted==0]<-1
View(d_sorted)
#Load in 'annotations" file. 
#First column must match sample ids from table above. 
#Second column contains grouping information
map<-read.delim("map.tab",row.names=1)
map=SLURRYmeta
View(map)

= SLURRYmeta$TreatmentInoculum
rownames(annotation_column) = SLURRYmeta$`#SampleID`
#draw the heatmap with log scaled intensities
#only top 1000 rows are plotted here. Personally I have used 4000 rows, but I try a few to see if it increases visibility
pheatmap(log10(d_sorted[1:12,]),color = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 
         
pheatmap(log10(d_sorted[1:12,]),colorRampPalette(c('black','blue'))(50),annotation=map)#since labels on legend are log values you can set them manually here AFTER visualizing (or just edit them in GIMP).


pheatmap(log10(d_sorted[1000,]),colorRampPalette(c('black','blue'))(50),show_rownames=FALSE,show_colnames=FALSE,annotation=map,legend_breaks=0:3,legend_labels=c(0,10,100,1000))





test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# Draw heatmaps

