rm(list=ls()) 
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/00.rawdata/")
if(!dir.exists("00.GSE194282")){dir.create("00.GSE194282")}
setwd("00.GSE194282")
getwd()
#
library(data.table)
b=fread('GPL17692-38203.txt',data.table = F)[,c(2,8)] #D读取下载到本地的soft文件
head(b)
library(stringr)
b$gene=str_split(b$gene_assignment,' // ',simplify = T)[,2] #将基因信息列按照 // 符号分离，抽取symble编号
ids=b[c(1,3)]
###########################
library(GEOquery)
library(limma)
library(affy)
gset <- getGEO('GSE194282', destdir=".",
               AnnotGPL = T,     ## 注释文件
               getGPL = F)       ## 平台文件
gset[[1]]
exp<-exprs(gset[[1]])
#cli<-pData(gset[[1]])	## 获取临床信息
#GPL <- fData(gset[[1]])	## 获取平台信息 
#gpl<-GPL[,c(1,3)]
#gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,
#                                     function(x)unlist(strsplit(x,"///"))[1]),
#                              stringsAsFactors=F)[,1]
################################################################################
#--------------------------------ids--------------------------------------------
library(data.table)
b=fread('GPL17692-38203.txt',data.table = F)[,c(2,8)] #D读取下载到本地的soft文件
head(b)
library(stringr)
b$gene=str_split(b$gene_assignment,' // ',simplify = T)[,2] #将基因信息列按照 // 符号分离，抽取symble编号
ids=b[c(1,3)]
################################################################################
#-----------------------------gene expression matrix-----------------------------
exp<-as.data.frame(exp)
exp$probeset_id <- rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息
exp_symbol<-merge(exp,ids,by="probeset_id")
exp_symbol<-na.omit(exp_symbol)

table(duplicated(exp_symbol$`gene`))
exp_unique <-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`gene`)
GSE194282_symbol_expr <- exp_unique
dim(GSE194282_symbol_expr)
write.csv(GSE194282_symbol_expr,"GSE194282_symbol_expr.csv")
################################################################################
#------------------------deg-----------------------------------------------------
rm(list=ls()) 
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/00.rawdata/")
if(!dir.exists("00.GSE194282")){dir.create("00.GSE194282")}
setwd("00.GSE194282")
getwd()
data1 <- read.csv("../00.GSE194282/GSE194282_symbol_expr.csv", row.names=1)
head(data1)
#list <- c(rep("AZS", 5), rep("Control",6)) %>% factor(., levels = c("Control", "AZS"), ordered = F)
group_list <- c(rep("Control", 7), rep("CRSwNP",7))
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
write.csv(DEG,file = "DEG_GSE194282.csv") 
##################################################################################################################
allDiff <- DEG
allDiff
allDiff$type <- ifelse(allDiff$logFC >= 0.5 & allDiff$P.Value < 0.05, "up",
                       ifelse(allDiff$logFC <= -0.5 & allDiff$P.Value < 0.05, "down", "not-sig")
)
allDiff$gene <- rownames(allDiff) 
allDiff
write.csv(allDiff,file = "sigDiff_GSE194282.csv")



