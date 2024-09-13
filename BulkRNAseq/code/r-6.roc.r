rm(list = ls())

setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if (!dir.exists("06.roc")){dir.create("06.roc")}
setwd("06.roc")
library(tidyverse)
library(lance)
#library(ROCR)
library(ggplot2)
library(pROC)

###读取和整理训练集train和验证集validate的ROC输入数据
dat_train <- read.csv('../01.DEGs/SYMBOL_Expre_data.csv',row.names = 1,header = T)|>t()|>data.frame()
dat_validate <- read.csv('../00.rawdata/00.GSE194282/GSE194282_symbol_expr.csv',row.names = 1,header = T)|>t()|>data.frame()
group_train <- read.csv('../01.DEGs/colData.csv',header = T)
group_validate <- read.csv('../00.rawdata/00.GSE194282/GSE194282_group.csv',header = T)
table(group_train$group)
# control involved 
# 64       58 
table(group_validate$group)
# control involved 
# 21       33 
condition_train <- group_train[match(rownames(dat_train),group_train$sample),2]
condition_validate <- group_validate[match(rownames(dat_validate),group_validate$sample),2]
#
#-----------------------------------hub gene------------------------------------
hub_gene <- read.csv("../05.ML/05.hub_gene.csv")
KeyGene <- hub_gene$symbol
#lasso_gene <- read.csv('../08_Machine_Learning/01.lasso_gene_adj.csv',skip=1,header =F)
#svm_gene <- read.csv('../08_Machine_Learning/02.svm_gene_adj.csv',skip=1,header = F)
#lasso_interset_svm <- intersect(lasso_gene$V1,svm_gene$V1)
#"LCE3D" "PI3"  
dat_train <- dat_train[,colnames(dat_train) %in% KeyGene]
dat_validate <- dat_validate[,colnames(dat_validate) %in% KeyGene]
dat_train$condition <- condition_train
dat_validate$condition <- condition_validate


# lasso_union_svm <- union(lasso_gene$V1,svm_gene$V1)
# dat_train_union <- dat_train[,colnames(dat_train) %in% lasso_union_svm]
# dat_validate_union <- dat_validate[,colnames(dat_validate) %in% lasso_union_svm]
# dat_train_intersect$condition <- condition_train
# dat_validate_union$condition <- condition_validate


###绘制训练集trainROC曲线-----------
project <- 'train'
for (i in c(1:length(KeyGene))) {
  roc<-roc(dat_train$condition,dat_train[,i],levels=c("CRSwNP", "Control"),direction='>')
  png(paste0(i, ".", colnames(dat_train)[i],'-',project,".png"),width = 4,height = 4,units = 'in',res = 600)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       auc.polygon=T,##将area显示为多边形
       #auc.polygon.con="#fff7f7",
       grid=c(0.2,0.2),#添加背景网格，横轴以0.2为间隔在x轴上添加网格线，以0.2为间隔在y轴上添加网格线
       grid.col=c("black","black"),
       print.thres=T,##添加戳点cutoff和95%CI
       main=paste0(colnames(dat_train)[i]),
       col="#FF2E63",
       legacy.axes=T)##更爱x轴格式)
  dev.off()
  pdf(paste0(i, ".", colnames(dat_train)[i],'-',project,".pdf"),width = 4,height = 4)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       auc.polygon=T,##将area显示为多边形
       #auc.polygon.con="#fff7f7",
       
       grid=c(0.2,0.2),#添加背景网格，横轴以0.2为间隔在x轴上添加网格线，以0.2为间隔在y轴上添加网格线
       grid.col=c("black","black"),
       print.thres=T,##添加戳点cutoff和95%CI
       main=paste0(colnames(dat_train)[i]),
       col="#FF2E63",
       legacy.axes=T)##更爱x轴格式)
  dev.off()
}


###绘制验证集validateROC曲线-----------
project <- 'validate'
for (i in c(1:length(KeyGene))) {
  roc<-roc(dat_validate$condition,dat_validate[,i],levels=c("CRSwNP", "Control"),direction='>')
  png(paste0(i, ".", colnames(dat_validate)[i],'-',project,".png"),width = 4,height = 4,units = 'in',res = 600)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       auc.polygon=T,##将area显示为多边形
       #auc.polygon.con="blue",
       
       grid=c(0.2,0.2),#添加背景网格，横轴以0.2为间隔在x轴上添加网格线，以0.2为间隔在y轴上添加网格线
       grid.col=c("black","black"),
       print.thres=T,##添加戳点cutoff和95%CI
       main=paste0(colnames(dat_validate)[i]),
       col="#FF2E63",
       legacy.axes=T)##更爱x轴格式)
  dev.off()
  
  pdf(paste0(i, ".", colnames(dat_validate)[i],'-',project,".pdf"),width = 4,height = 4)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       auc.polygon=T,##将area显示为多边形
       #auc.polygon.con="blue",
       
       grid=c(0.2,0.2),#添加背景网格，横轴以0.2为间隔在x轴上添加网格线，以0.2为间隔在y轴上添加网格线
       grid.col=c("black","black"),
       print.thres=T,##添加戳点cutoff和95%CI
       main=paste0(colnames(dat_validate)[i]),
       col="#FF2E63",
       legacy.axes=T)##更爱x轴格式)
  dev.off()
}

write.csv(KeyGene,file='biomarker.csv',quote = F,row.names = F)  
