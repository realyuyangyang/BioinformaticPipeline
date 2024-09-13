rm(list = ls())
#options(stringsAsFactors = F)
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if(!dir.exists("03.Venn")){dir.create("03.Venn")}
setwd("03.Venn")
getwd()
#
########################library#################################################
library(dplyr)
library(ggplot2)
library(ggvenn)
library(VennDiagram)
library(ggtext)
###########################Data##################################################
#signif_DEG <- read.csv("/data/nas1/yuyangyang/AZS/ProcessedData/signif_DEG.csv", row.names=1)
DEG_sig <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/01.DEGs/DEG_sig.csv")
DEG <-DEG_sig$X
DEG
MAG <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/03.Venn/MAG.csv")
MAG
MPRGmodule <- read.table("/data/nas1/yuyangyang/Projects/YQKM-10411-2/02.WGCNA/13.module.txt", quote="\"", comment.char="")
#MAG<- read.table("/data/nas1/yuyangyang/AZS/RawData/MitochondrioGnene.csv", 
#                 quote="\"", comment.char="")

#MPRGmodule<- read.table("/data/nas1/yuyangyang/AZS/RawData/PCD_RG.csv", 
#                 quote="\"", comment.char="")
MAG<- MAG$Gene.Symbol
MPRGmodule<- MPRGmodule$V1
#MAG_MPRGmodule_inter <- intersect(MAG, MPRGmodule)
#genelist<-intersect(DEG1,MAG_MPRGmodule_inter)
#genelist
gene1<-intersect(DEG,MAG)
gene2<-intersect(DEG,MPRGmodule)
#genelist<-intersect(gene1,gene2)
Candidate_genelist <-intersect(gene1,gene2)
Candidate_genelist
write.csv(Candidate_genelist,"intersect.gene.csv")
saveRDS(Candidate_genelist,"intersect.gene.rds")
###################################Venn
venn_list<-list(DEG=DEG,
                MAG=MAG,
                MPRGmodule=MPRGmodule)
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
Candidated_genes<-inter[1,]$..values.. $`1`
Candidated_genes<-Candidated_genes[order(Candidated_genes,decreasing = F )]
#colnames(Candidated_genes) <- c("symbol")
#Candidated_genes
write.table(Candidated_genes,
            'Candidated_genes.csv',
            row.names = FALSE,
            col.names="symbol", 
            sep = '\t', 
            quote = FALSE)
##################################Plot##########################################
p1 <- ggvenn(venn_list,
             c('DEG','MAG','MPRGmodule'),
             fill_color =c('#DDA0DD', '#90EE90','#87CEEB'),
             fill_alpha = 0.5,
             stroke_alpha = 1,
             stroke_size = 0.4,
             text_size = 5,
             stroke_color=NA,
             stroke_linetype="solid",
             set_name_color=c('#DDA0DD', '#90EE90','#87CEEB'),
             set_name_size = 6,
             text_color = 'black'
)+
  xlim(-2.5,2.5)+
  ylim(-3,2)
p1
ggsave('Venn.png',width=16,height=12)
ggsave('Venn.jpg',width=16,height=12, dpi = 600)
ggsave('Venn.pdf',width=16,height=12, dpi = 600)
dev.off()
#pdf('02.Candidated_genes.pdf',w=8,h=8,family='Times')

