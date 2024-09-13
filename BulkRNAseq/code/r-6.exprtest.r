#候选基因筛选
rm(list=ls()) 
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if(!dir.exists("06.exprtest")){dir.create("06.exprtest")}
setwd("06.exprtest")
getwd()
######################################################
#####################train##########################
library(rstatix)
library(ggplot2)
library(ggpubr)
library(reshape)
########################################################
#-----------------------DATA---------------------------
Expre <- read.csv("../01.DEGs/SYMBOL_Expre_data.csv", 
                  row.names=1)
Expre <- log2(Expre+1)
range(Expre)
Expre
dif <- read.csv("../01.DEGs/DEG_sig.csv", 
                row.names=1)
train_group <- read.csv("../01.DEGs/colData.csv")
Candidated_genes <-read.csv("../05.ML/05.hub_gene.csv", row.names = NULL)
#Candidated_genes <- c("ALOX5", "GRIA2", "GNMT", "BCAT1") %>% as.data.frame()
#names(Candidated_genes) <- c("symbol")
#Candidated_genes <- Candidated_genes$symbol
dim(Candidated_genes)
#
#
key_eset<-Expre[Candidated_genes$symbol,]%>%na.omit()%>%t()
key_eset<-key_eset[,order(colnames(key_eset),decreasing = F)]
##按照字母进行排序
dat <- key_eset %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,train_group$sample)
dat <- merge(dat, train_group, by = "sample")
dat2 <- tidyr::gather(dat, symbol, expression, -c(sample, group))
dat2
library(rstatix)
stat_res <- dat2 %>% 
  group_by(symbol) %>% 
  wilcox_test(expression ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat_res 
#write.csv(stat_res, "stat_res.csv",quote=F)
hub_dif<-dif[stat_res $symbol,]
hub_dif
write.csv(hub_dif, "hub_dif_train.csv",quote=F)
stat_res$p<-hub_dif$P.Value  #利用差异分析P值
stat_res

stat_res$p.signif[stat_res$p>0.05] <- "ns"   
stat_res$p.signif[stat_res$p<0.05&stat_res$p>0.01] <- "*"   
stat_res$p.signif[stat_res$p<0.01&stat_res$p>0.001] <- "**"  
stat_res$p.signif[stat_res$p<0.001&stat_res$p>0.0001] <- "***"   
stat_res$p.signif[stat_res$p<0.0001] <- "****"
colnames(stat_res)[1] <- 'Gene'
stat_res
##boxplot
boxplot_dat <-data.frame(cbind(key_eset,train_group$group))
boxplot_dat
colnames(boxplot_dat)[6]<-'group'##############################################gene number+1
head(boxplot_dat)
boxplot_dat$id <- rownames(boxplot_dat)
boxplot_dat <- melt(boxplot_dat,id = c("id","group"))
boxplot_dat$value <-as.numeric(boxplot_dat$value)
boxplot_dat$Group<- factor(boxplot_dat$group)
boxplot_dat
colnames(boxplot_dat)<-c('id','Group','Gene','value')
boxplot_dat$Gene<-gsub('[.]','-',boxplot_dat$Gene)
boxplot_dat

##boxplot
##boxplot

###################Data Prepration###################
library(ggsci)
library(ggplot2)
## 调整因子水平
boxplot_dat$Group <- factor(boxplot_dat$Group,levels = c("Control","CRSwNP"))
p<-ggboxplot(
  boxplot_dat,
  x = "Group",
  y = "value",
  fill = "Group", 
  palette =  c( "#D85F06","#1E9C79"),
  # add = "jitter",
  # ylim = c(0, 10)
)+stat_pvalue_manual(stat_res,
                     y.position = c(rep(13, 5)),###columns of Gene,
                     size = 4.5,
                     color = "black",
                     family = "Times",
                     label = "p.signif",
                     #parse = T,
                     face = "bold")+
  ggtitle('GSE136825')+
  theme_bw() + xlab("")+ylab("Expression of genes")+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.x =element_text(size=16,family = "Times", face = "bold",color='black'), #####Group
        axis.title.y =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.y=element_text(size=10,family = "Times", face = "bold",color='black'),
        strip.text = element_text(size = 16,family = "Times", face = "bold",color='black'),####Gene
        plot.title=element_text(size=24,family = "Times", face = "bold",hjust=0.5,color='black'),
        legend.position = "none")+
  facet_wrap(.~Gene,nrow=1
             ,scales = "free_y"
  )+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), legend.title=element_text(size=15) , legend.text=element_text(size=14))
p
#png('04.exprtest_train_gene_boxplot.png',width=12,height=5,units='in',res=600,bg='white')
#pdf('04.exprtest_traingene_boxplot.pdf',width=12,height=5,bg='white')
ggsave('06.exprtest_train_gene_boxplot.png',width=11,height=6)
ggsave('06.exprtest_train_gene_boxplot.jpg',width=11,height=6, dpi = 600)
ggsave('06.exprtest_train_gene_boxplot.pdf',width=11,height=6, dpi = 600)
print(p)
dev.off()
###############################Save#######################################################################
#------------------------------GSE194282---------------------------------------------
###############################validate###########################################################
#候选基因筛选

######################################################
#####################train##########################
library(rstatix)
library(ggplot2)
library(ggpubr)
library(reshape)
########################################################
#-----------------------DATA---------------------------
Expre <- read.csv("../00.rawdata/00.GSE194282/GSE194282_symbol_expr.csv", 
                  row.names=1)
Expre <- log2(Expre+1)
range(Expre)
Expre
dif <- read.csv("../00.rawdata/00.GSE194282/sigDiff_GSE194282.csv", 
                row.names=1)
train_group <- read.csv("../00.rawdata/00.GSE194282/GSE194282_group.csv")

#train_group<-data.frame(sample=rownames(Group_all),Group=Group_all$Group)
########------------------------------------------------------------------------
key_eset<-Expre[Candidated_genes$symbol,]%>%na.omit()%>%t()
key_eset<-key_eset[,order(colnames(key_eset),decreasing = F)]
##按照字母进行排序
dat <- key_eset %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,train_group$sample)
dat <- merge(dat, train_group, by = "sample")
dat2 <- tidyr::gather(dat, symbol, expression, -c(sample, group))
dat2
library(rstatix)
stat_res <- dat2 %>% 
  group_by(symbol) %>% 
  wilcox_test(expression ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat_res 
#write.csv(stat_res, "stat_res.csv",quote=F)
hub_dif<-dif[stat_res $symbol,]
hub_dif
write.csv(hub_dif, "hub_dif_train.csv",quote=F)
stat_res$p<-hub_dif$P.Value  #利用差异分析P值
stat_res

stat_res$p.signif[stat_res$p>0.05] <- "ns"   
stat_res$p.signif[stat_res$p<0.05&stat_res$p>0.01] <- "*"   
stat_res$p.signif[stat_res$p<0.01&stat_res$p>0.001] <- "**"  
stat_res$p.signif[stat_res$p<0.001&stat_res$p>0.0001] <- "***"   
stat_res$p.signif[stat_res$p<0.0001] <- "****"
colnames(stat_res)[1] <- 'Gene'
stat_res
##boxplot
boxplot_dat <-data.frame(cbind(key_eset,train_group$group))
boxplot_dat
colnames(boxplot_dat)[6]<-'group'##############################################gene number+1
head(boxplot_dat)
boxplot_dat$id <- rownames(boxplot_dat)
boxplot_dat <- melt(boxplot_dat,id = c("id","group"))
boxplot_dat$value <-as.numeric(boxplot_dat$value)
boxplot_dat$Group<- factor(boxplot_dat$group)
boxplot_dat
colnames(boxplot_dat)<-c('id','Group','Gene','value')
boxplot_dat$Gene<-gsub('[.]','-',boxplot_dat$Gene)
boxplot_dat

##boxplot
##boxplot

###################Data Prepration###################
library(ggsci)
library(ggplot2)
boxplot_dat$Group <- factor(boxplot_dat$Group,levels = c("Control","CRSwNP"))
p<-ggboxplot(
  boxplot_dat,
  x = "Group",
  y = "value",
  fill = "Group", 
  palette =  c( "#D85F06","#1E9C79"),
  # add = "jitter",
  # ylim = c(0, 10)
)+stat_pvalue_manual(stat_res,
                     y.position = c(rep(3.3, 5)),###columns of Gene,
                     size = 4.5,
                     color = "black",
                     family = "Times",
                     label = "p.signif",
                     #parse = T,
                     face = "bold")+
  ggtitle('GSE194282')+
  theme_bw() + xlab("")+ylab("Expression of genes")+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.x =element_text(size=16,family = "Times", face = "bold",color='black'), #####Group
        axis.title.y =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.y=element_text(size=10,family = "Times", face = "bold",color='black'),
        strip.text = element_text(size = 16,family = "Times", face = "bold",color='black'),####Gene
        plot.title=element_text(size=24,family = "Times", face = "bold",hjust=0.5,color='black'),
        legend.position = "none")+
  facet_wrap(.~Gene,nrow=1
             ,scales = "free_y"
  )+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), legend.title=element_text(size=15) , legend.text=element_text(size=14))
p
ggsave('06.GSE194282.png',width=11,height=6)
ggsave('06.GSE194282.jpg',width=11,height=6, dpi = 600)
ggsave('06.GSE194282.pdf',width=11,height=6, dpi = 600)
dev.off()
#png('04.exprtest_train_gene_boxplot.png',width=12,height=5,units='in',res=600,bg='white')
#pdf('04.exprtest_traingene_boxplot.pdf',width=12,height=5,bg='white')
###############################Save#######################################################################
#------------------------------GSE36830---------------------------------------------
###############################validate###########################################################
#候选基因筛选

######################################################
#####################train##########################
library(rstatix)
library(ggplot2)
library(ggpubr)
library(reshape)
########################################################
#-----------------------DATA---------------------------
Expre <- read.csv("../00.rawdata/GSE36830.csv", 
                  row.names=1)
Expre <- log2(Expre+1)
range(Expre)
Expre
dif <- read.csv("../00.rawdata/DEG_GSE36830.csv", 
                row.names=1)
train_group <- read.csv("../00.rawdata/validate_group.csv")

#train_group<-data.frame(sample=rownames(Group_all),Group=Group_all$Group)
########------------------------------------------------------------------------
key_eset<-Expre[Candidated_genes$symbol,]%>%na.omit()%>%t()
key_eset<-key_eset[,order(colnames(key_eset),decreasing = F)]
##按照字母进行排序
dat <- key_eset %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,train_group$sample)
dat <- merge(dat, train_group, by = "sample")
dat2 <- tidyr::gather(dat, symbol, expression, -c(sample, group))
dat2
library(rstatix)
stat_res <- dat2 %>% 
  group_by(symbol) %>% 
  wilcox_test(expression ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat_res 
#write.csv(stat_res, "stat_res.csv",quote=F)
hub_dif<-dif[stat_res $symbol,]
hub_dif
write.csv(hub_dif, "hub_dif_train.csv",quote=F)
stat_res$p<-hub_dif$P.Value  #利用差异分析P值
stat_res

stat_res$p.signif[stat_res$p>0.05] <- "ns"   
stat_res$p.signif[stat_res$p<0.05&stat_res$p>0.01] <- "*"   
stat_res$p.signif[stat_res$p<0.01&stat_res$p>0.001] <- "**"  
stat_res$p.signif[stat_res$p<0.001&stat_res$p>0.0001] <- "***"   
stat_res$p.signif[stat_res$p<0.0001] <- "****"
colnames(stat_res)[1] <- 'Gene'
stat_res
##boxplot
boxplot_dat <-data.frame(cbind(key_eset,train_group$group))
boxplot_dat
colnames(boxplot_dat)[6]<-'group'##############################################gene number+1
head(boxplot_dat)
boxplot_dat$id <- rownames(boxplot_dat)
boxplot_dat <- melt(boxplot_dat,id = c("id","group"))
boxplot_dat$value <-as.numeric(boxplot_dat$value)
boxplot_dat$Group<- factor(boxplot_dat$group)
boxplot_dat
colnames(boxplot_dat)<-c('id','Group','Gene','value')
boxplot_dat$Gene<-gsub('[.]','-',boxplot_dat$Gene)
boxplot_dat

##boxplot
##boxplot

###################Data Prepration###################
library(ggsci)
library(ggplot2)
boxplot_dat$Group <- factor(boxplot_dat$Group,levels = c("Control","CRSwNP"))
p<-ggboxplot(
  boxplot_dat,
  x = "Group",
  y = "value",
  fill = "Group", 
  palette =  c( "#D85F06","#1E9C79"),
  # add = "jitter",
  # ylim = c(0, 10)
)+stat_pvalue_manual(stat_res,
                     y.position = c(rep(3.3, 5)),###columns of Gene,
                     size = 4.5,
                     color = "black",
                     family = "Times",
                     label = "p.signif",
                     #parse = T,
                     face = "bold")+
  ggtitle('GSE36830')+
  theme_bw() + xlab("")+ylab("Expression of genes")+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.x =element_text(size=16,family = "Times", face = "bold",color='black'), #####Group
        axis.title.y =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.y=element_text(size=10,family = "Times", face = "bold",color='black'),
        strip.text = element_text(size = 16,family = "Times", face = "bold",color='black'),####Gene
        plot.title=element_text(size=24,family = "Times", face = "bold",hjust=0.5,color='black'),
        legend.position = "none")+
  facet_wrap(.~Gene,nrow=1
             ,scales = "free_y"
  )+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), legend.title=element_text(size=15) , legend.text=element_text(size=14))
p
ggsave('06.GSE36830.png',width=11,height=6)
ggsave('06.GSE36830.jpg',width=11,height=6, dpi = 600)
ggsave('06.GSE36830.pdf',width=11,height=6, dpi = 600)
dev.off()
