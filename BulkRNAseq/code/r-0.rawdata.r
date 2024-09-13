#DEGs
rm(list=ls())
getwd()
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/") 
if (!dir.exists("00.rawdata")) {dir.create("00.rawdata")}
setwd("00.rawdata")
library(dplyr)
library(tidyverse)
library('Biobase')
library(org.Hs.eg.db)#人类
library(clusterProfiler)
#-----------------------------Counts-> Tpm--------------------------------------------------
#DEGs
expr_df <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/00.rawdata/GSE136825_genecounts_length.csv", 
                                        row.names=NULL)

names(expr_df)
counts <-  expr_df[,3:ncol(expr_df)]
dim(counts)
expr_df <- expr_df[which(rowSums(counts) > 0),]
expr_df
dim(expr_df)
counts <-  expr_df[,3:ncol(expr_df)]
dim(counts)
#counts <- counts[rowSums(counts)>0,]
#基因长度，目标基因的外显子长度之和除以1000，单位是Kb，不是bp
kb <- expr_df$Length / 1000
rpk <- counts / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
dim(tpm)
rownames(tpm) <- expr_df$Geneid
write.csv(tpm,file = "Ensemble_Tpm.csv")

################################################################################
#-----------------------------------ID Convert----------------------------------
Ensemble_Tpm <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/00.rawdata/Ensemble_Tpm.csv", 
                                 row.names=1)
data <- Ensemble_Tpm
data$Ensembl_ID <- substr(rownames(data),1,15)

gene.df <- bitr(data$Ensembl_ID , fromType = "ENSEMBL",
                toType = c("SYMBOL"),
                OrgDb = org.Hs.eg.db)
SYMBOL_lst <- gene.df$SYMBOL
SYMBOL_lst

data <- data[match(gene.df$ENSEMBL,data$Ensembl_ID),]
dim(data)
identical(data$Ensembl_ID,gene.df$ENSEMBL)
data$Ensembl_ID <- gene.df$SYMBOL
data
names(data)[names(data)=="Ensembl_ID"] <- "symbol"
unique_SYMBOL_data <- data[!duplicated(data$symbol), ]
unique_SYMBOL_data
dim(unique_SYMBOL_data)


getwd()
write.csv(unique_SYMBOL_data,file = "Symbol_Tpm.csv")
##########################################################################################
#-----------------------------validate data-----------------------------------------------#
#GEO芯片数据下载和探针ID转换: https://blog.csdn.net/qq_43138237/article/details/128008548

options(stringAsFactors = F)	# # 避免引入factor
#Sys.setenv("VROOM_CONNECTION_SIZE"=131072*600)	# 设置内存，避免不够用
library(GEOquery)
library(limma)
library(affy)
gset <- getGEO('GSE36830', destdir=".",
                               AnnotGPL = T,     ## 注释文件
                               getGPL = T)       ## 平台文件
gset[[1]]
exp<-exprs(gset[[1]])
cli<-pData(gset[[1]])	## 获取临床信息
GPL <- fData(gset[[1]])	## 获取平台信息 
gpl<-GPL[,c(1,3)]
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,
                                     function(x)unlist(strsplit(x,"///"))[1]),
                              stringsAsFactors=F)[,1]
exp<-as.data.frame(exp)
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol<-na.omit(exp_symbol)
table(duplicated(exp_symbol$`Gene symbol`))
exp_unique <-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`Gene symbol`)
validate_data <- exp_unique[,-(7:12)]
dim(validate_data)
write.csv(validate_data,"GSE36830.csv")

dat <- log2(exp_unique+1)
data <- log(data+1)
boxplot(data.frame(data), col = "red")
######################################################################################
#-------------------------------------validate deg-----------------------------------
#library
##################################################################################
##################################差异分析########################################
##
#write.csv(nrDEG,file = "C:/Users/10149/Desktop/azs_20240325/dataset/ProcessedData/Symbol_Azs_Ctrl_DEG.csv") 

#Symbol_Azs_Ctrl_Expre <- read.csv("C:/Users/10149/Desktop/azs_20240325/dataset/ProcessedData/Symbol_Azs_Ctrl_Expre.csv", row.names=1)
rm(list=ls())
getwd()
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/") 
if (!dir.exists("00.rawdata")) {dir.create("00.rawdata")}
setwd("00.rawdata")
data1 <- read.csv("../00.rawdata/GSE36830.csv", row.names=1)
head(data1)
#list <- c(rep("AZS", 5), rep("Control",6)) %>% factor(., levels = c("Control", "AZS"), ordered = F)
group_list <- c(rep("Control", 6), rep("CRSwNP",12))
group_list <- factor(group_list,levels = c("Control","CRSwNP"),ordered = F)
group_list
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
design
## 比较信息
contrast.matrix<-makeContrasts("CRSwNP-Control",
                               levels = design)
contrast.matrix##查看比较矩阵的信息，这里我们设置的是Treat Vs. control
## 拟合模型
eset_dat <- data1
fit <- lmFit(eset_dat,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
head(DEG)
write.csv(DEG,file = "DEG_GSE36830.csv") 
##################################################################################################################
allDiff <- DEG
allDiff
allDiff$type <- ifelse(allDiff$logFC > 0.5 & allDiff$P.Value < 0.05, "up",
                       ifelse(allDiff$logFC < -0.5 & allDiff$P.Value < 0.05, "down", "not-sig")
)
allDiff$gene <- rownames(allDiff) 
allDiff
write.csv(allDiff,file = "sigDiff_GSE36830.csv")






