#DEGs
rm(list=ls())
getwd()
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/") 
if (!dir.exists("02.WGCNA")) {dir.create("02.WGCNA")}
setwd("02.WGCNA")

#------------------------------DATA---------------------------------------------
expr <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/00.rawdata/Symbol_Tpm.csv",row.names = 1)
dim(expr)
range(expr)
expr <- log2((expr)+1)
range(expr)
MPRG <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/02.WGCNA/MP-RGs.csv")
group = read.csv('/data/nas1/yuyangyang/Projects/YQKM-10411-2/02.WGCNA/group.csv')
#------------------------------packages----------------------------------------
library(dplyr)
library(ggplot2)
library(limma)
library(WGCNA)
library(stringr)
library(GSVA)
library(magrittr)
library(tidyr)
library(ggpubr)
library(tibble)
library(ggsci)
# ssgsea ------------------------------------------------------------------
gene.set.list<-data.frame(MPRG$gene)
ssgsea_res <-
  gsva(
    as.matrix(expr),
    gene.set.list,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = T
  )
ssgsea_score<-data.frame(ssgsea_res)%>%t(.)
write.csv(ssgsea_score,'01.ssgsea_score.csv',quote=F)

  
dat <- ssgsea_score  %>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
identical(dat$sample,group$sample)
dat <- merge(dat, group, by = "sample")
dat2 <- tidyr::gather(dat, Gene, Score, -c(sample, group))


library(rstatix)
stat_res <- dat2 %>%
  group_by(Gene) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat_res
write.table(stat_res, file = "02.ssgsea_wilcoxon_res.xls", sep = "\t", row.names = F, quote = F)

library(reshape2)
boxplot_dat <-data.frame(cbind(ssgsea_score,group$group))
colnames(boxplot_dat)[2]<-'group'
head(boxplot_dat)
boxplot_dat$id <- rownames(boxplot_dat)
boxplot_dat <- melt(boxplot_dat,id = c("id","group"))
boxplot_dat$value <-as.numeric(boxplot_dat$value)
boxplot_dat$group<- factor(boxplot_dat$group)
colnames(boxplot_dat)<-c('id','group','Gene','value')

library(ggsci)
library(ggplot2)
library(ggpubr)

p<-ggviolin(
  boxplot_dat,
  x = "group",
  y = "value",
  fill = "group", 
  palette = c('#104680','#B72230'),
  add = "boxplot",
  add.params = list(fill = "white")
)+
  stat_pvalue_manual(stat_res,
                     # x = "Gene",
                     label = "{p}{p.signif}", 
                     #tip.length = 0, 
                     #hide.ns = TRUE,
                     y.position = 5, #bar position
                     size = 5,
                     face = "bold",
                     color='red',
                     family = "Times"
  )+
  theme_bw() + xlab("")+ylab("ssgseaScore of MP-RGs")+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.x =element_text(size=16,family = "Times", face = "bold",color='black'),
        axis.title.y =element_text(size=20,family = "Times", face = "bold",color='black'),
        axis.text.y=element_text(size=16,family = "Times", face = "bold",color='black'),
        strip.text = element_text(size = 14,family = "Times", face = "bold",color='black'),
        legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), legend.title=element_text(size=15) , legend.text=element_text(size=14))
p
ggsave(filename = "03.ssgsea_score_violin.pdf", width = 8, height = 6, p)
ggsave(filename = "03.ssgsea_score_violin.png", width = 8, height = 6, dpi = 600, p)

#-----------------------------------------------------------------------
expro<-na.omit(expr)
datExprOri=as.data.frame(t(expro))
nGenes0 = ncol(datExprOri)
nSamples0 = nrow(datExprOri)
nGenes0 ;nSamples0 
dim(datExprOri)  ##332个sample   10000个基因
#m.mad <- apply(eset,1,mad)
#dataExprVar <-eset[which(m.mad >max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

#转换为样品在行，基因在列的矩阵
# expro<-eset
# datExprOri=as.data.frame(t(expro))
nGenes0 = ncol(datExprOri)
nSamples0 = nrow(datExprOri)
nGenes0 ;nSamples0 
## 过滤缺失基因数目多于10%的样本
library(WGCNA)
##检验数据质量
gsg = goodSamplesGenes(datExprOri , verbose = 3)   #检测缺失值
gsg$allOK   ##为TRUE不用删除基因
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(datExprOri)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:",
                     paste(rownames(datExprOri)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  datExprOri = datExprOri[gsg$goodSamples, gsg$goodGenes]
  
}      ##输出的datExprOri已经过滤基因
datExpr1<-datExprOri   
nGenes1<-ncol(datExpr1)      #基因数目
nSamples1 <- nrow(datExpr1) #样品名称变量
nGenes1 ;nSamples1  

##根据层次聚类绘制样本
# 观察是否有离群样品需要剔除
tree=hclust(dist(datExpr1),method ='mcquitty')   ##观察是否有离群样本  c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
pdf(file='04.sampleClustering.pdf',w=12,h=8)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif') 
plot(tree,xlab="", sub="", main="Sample Clustering",
     # labels=FALSE,
     cex=1.0,
     font=2,
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
abline(h =290, col = 'red')
dev.off()

png(file='04.sampleClustering.png',w=12,h=8,units='in',res=600,bg='white')
par(pin = c(4,4), mar = c(6,6,6,1),family='serif') 
plot(tree,xlab="", sub="", main="Sample Clustering",
     # labels=FALSE,
     cex=1.0,
     font=2,
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
abline(h =290, col = 'red')
dev.off()

## 从上图可以看出，疑似异样，考虑去除。
clust = cutreeStatic(tree, cutHeight =290, minSize = 10) 
table(clust)  ###无离群样本
## 去除异常样品后的表达矩阵
datExpr = datExpr1#[clust == 1, ]
nGenes = ncol(datExpr)      #基因数目
nSample =nrow(datExpr) #样品名称变量
nGenes ;nSample   
SampleName<-rownames(datExpr )

## 表型数据
datTraits<-data.frame(row.names=rownames(ssgsea_score),ssgsea_score=ssgsea_score[,1],sample=rownames(ssgsea_score))
datTraits<-datTraits[SampleName,]
datTraits<-data.frame(row.names=datTraits$sample,ssGSEA_score=datTraits$ssgsea_score)
identical(rownames(datExpr),rownames(datTraits))
datExpr[1:4,1:4]
datTraits[1:4,]

traitColors = numbers2colors(datTraits, signed = FALSE)
tree2=hclust(dist(datExpr),method ='mcquitty')   ##观察是否有离群样本  c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")

pdf(file='05.sampleClustering2.pdf',w=10,h=7)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif') 
plotDendroAndColors(tree2,
                    traitColors, 
                    #dendroLabels = FALSE, 
                    hang = 0.03,
                    cex.dendroLabels = 0.8, 
                    # addGuide = TRUE,   ##添加网格线
                    font=2,
                    guideHang = 0.05,
                    cex.axis=1.4,  ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.lab=1.6,   ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.main=1.6,   ##标题的缩放倍数。
                    font.axis = 2, 
                    font.lab = 2, 
                    font.main = 2, 
                    font.sub =2,
                    cex.colorLabels=1.0,
                    groupLabels = colnames(datTraits), 
                    main = "Sample Clustering and trait heatmap")
dev.off()

png(file='05.sampleClustering2.png',w=10,h=7,units='in',res=600,bg='white')
par(pin = c(4,4), mar = c(6,6,6,1),family='serif') 
plotDendroAndColors(tree2,
                    traitColors, 
                    # dendroLabels = FALSE, 
                    hang = 0.03,
                    cex.dendroLabels = 0.8, 
                    # addGuide = TRUE,   ##添加网格线
                    font=2,
                    guideHang = 0.05,
                    cex.axis=1.4,  ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.lab=1.6,   ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.main=1.6,   ##标题的缩放倍数。
                    font.axis = 2, 
                    font.lab = 2, 
                    font.main = 2, 
                    font.sub =2,
                    cex.colorLabels=1.0,
                    groupLabels = colnames(datTraits), 
                    main = "Sample Clustering and trait heatmap")
dev.off()


# 一步构建共表达网络和划分模块 
##设置网络构建参数选择范围，计算beta值（为了计算power）
powers <- c(seq(1, 10, by=1), seq(12, 20, by=2)) #设置beta值的取值范围
#rsquare.cut=0.85  #设置r2阈值为0.85，一般选择0.85或0.9

#开启多线程，加快运行速度
enableWGCNAThreads()
#disbleWGCNAThreads()     #如果运行pickSoftThreshold报call()相关的错误，需要运行禁用多线程
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
#软阈值选择绘图
#左图为每个候选beta对应的R2值,
#绘制阈值线0.9（红线）和0.85（蓝线），如没有超过蓝线，则没找到软阈值，sft$powerEstimate为空
#右图为每个候选beta对应的平均连接度
pdf('06.softThreshold.pdf',w=12,h=8)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif') 
cex1=1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"),
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

png('06.softThreshold.png',w=12,h=8,units='in',res=600,bg='white')
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif') 
cex1=1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"),
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，
# 数值越高，网络越符合无标度特征 (non-scale)。
##自动构建WGCNA模型
#表达矩阵转换成邻接矩阵，然后再将邻接矩阵转换成拓扑矩阵,识别模块
cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power =sft$powerEstimate,
  minModuleSize =100,                    #每个模块最少的基因数
  deepSplit = 4,                         #剪切树参数，一般设置为2，取值0-4（最灵敏）
  mergeCutHeight = 0.4,                 #模块合并参数，越大模块越少，默认0.25  越小，模块越多
  numericLabels = TRUE,                  # T返回数字，F返回颜色
  networkType  = "signed",               #网络类型，一般不需要修改
  maxBlockSize = ncol(datExpr),     
  pamRespectsDendro = FALSE, 
  saveTOMs = TRUE,
  saveTOMFileBase = "TPM-TOM",
  loadTOMs = TRUE,
  verbose = 3
)
cor<-stats::cor
table(net$colors)     ## 18个模块

##模块重新排序
#获取模块颜色
moduleColors = labels2colors(net$colors)
#获取每个模块的特征值
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); #重新排序，相似颜色的挨在一起
useMEs=subset(MEs, select = -c(MEgrey))  #去掉未分组的基因(灰色模块）

#useMEs=subset(MEs)  #去掉未分组的基因(灰色模块）
## 输出每个基因所在的模块，以及与该模块的KME值
if(file.exists('All_Gene_KME.txt')) file.remove('All_Gene_KME.txt')
for(module in substring(colnames(useMEs),3)){
  ME=as.data.frame(useMEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
}

##模块绘图
mergedColors = labels2colors(net$colors)
mergedColors[net$blockGenes[[1]]]
pdf("07.wgcna.dendroColors.pdf",height = 7,width = 9)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif') 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.main=1.8,   ##标题的缩放倍数。
                    font.axis = 2, 
                    font.lab = 2, 
                    font.main = 2, 
                    font.sub =2)

dev.off()

png("07.wgcna.dendroColors.png",height = 7,width = 9,units='in',res=600)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1),family='serif') 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.main=1.8,   ##标题的缩放倍数。
                    font.axis = 2, 
                    font.lab = 2, 
                    font.main = 2, 
                    font.sub =2)

dev.off()

# 将结果保存成cytoscape的输入文件格式
#读入构建时保存的TOM矩阵
load('TPM-TOM-block.1.RData')
TOM=as.matrix(TOM)
#save(file="WGCNA.RData",datExpr,SampleName,useMEs)


##模块与表型的关联分析
#load('WGCNA.RData')
#表型
moduleTraitCor = cor(useMEs, datTraits, use = "p") # use=p代表去掉缺失计算
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(SampleName))#
#整理要显示在图中的数字,第一行为相关系数，第二行为p值
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

#绘制关联热图
pdf("08.wgcna.Module-trait.heatmap.pdf", width = 8, height =12)
par(mar = c(14, 12, 3, 1.5),family='serif')  ##下左上右   ##par(pin = c(4,4), mar = c(6,6,6,1),family='serif') 
labeledHeatmap(Matrix = moduleTraitCor,   #相关系数
               xLabels = colnames(datTraits), #x轴为表型
               yLabels = names(useMEs),#y轴为模块
               ySymbols = names(useMEs),
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, #每个单元格的内容
               setStdMargins = FALSE,
               cex.text =1.1,cex.lab.y = 1.5,cex.legendLabel =1.5,cex.lab.x=1.5,
               zlim = c(-1,1),
               main = "Module-trait relationships",
               cex.lab = 1.2, 
               font.lab.x = 2, font.lab.y = 2 )
dev.off()

png("08.wgcna.Module-trait.heatmap.png", width = 8, height =12,unit='in',res=600,bg='white')
par(mar = c(14, 12, 3, 1.5),family='serif')  ##下左上右   ##par(pin = c(4,4), mar = c(6,6,6,1),family='serif') 
labeledHeatmap(Matrix = moduleTraitCor,   #相关系数
               xLabels = colnames(datTraits), #x轴为表型
               yLabels = names(useMEs),#y轴为模块
               ySymbols = names(useMEs),
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, #每个单元格的内容
               setStdMargins = FALSE,
               cex.text =1.1,cex.lab.y = 1.5,cex.legendLabel =1.5,cex.lab.x=1.5,
               zlim = c(-1,1),
               main = "Module-trait relationships",
               cex.lab = 1.2, 
               font.lab.x = 2, font.lab.y = 2 )
dev.off()

#基因共表达网络热图
dissTOM = 1-TOM
kME=signedKME(datExpr, useMEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
write.table(kME,"kME.txt",quote = F,sep = "\t",col.names = T)
if (dim(datExpr)[2]>=1500) nSelect=1500 else nSelect=dim(datExpr)[2]
set.seed(1)   #For reproducibility, we set the random seed
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
# 重新画聚类图
selectTree = hclust(as.dist(selectTOM),method ='complete')
selectColors = moduleColors[select]

#不同模块的基因显著性图
geneTraitSignificance  = as.data.frame(cor(datExpr, datTraits, use = "p"))
write.table(geneTraitSignificance,"GS.txt",quote = F,sep = "\t",row.names = T)
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 70))
names(geneTraitSignificance) = paste("GS.", colnames(datTraits), sep="")
names(GSPvalue) = paste("GS.", colnames(sample), sep="")
modNames = substring(names(MEs), 3)  ##选出两种颜色
save.image('temp.RData')


# select genes ------------------------------------------------------------
rm(list=ls())
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if (!dir.exists("02.WGCNA")){dir.create("02.WGCNA")}
setwd("02.WGCNA")
load('temp.RData')
modNames1<-c("red","brown","green","black")
###################################################################################################
gs <-0.5
# 
mm <-0.8
# 
# 
# for (module in modNames1){
#   if(module== "grey"){ next }
#   column = match(module, modNames); # col number of interesting modules
#   png(paste("09.GS_MM.", module, ".png", sep=""),height = 6,width = 7,family='Times',units='in',res=600)
#   moduleGenes = moduleColors==module;
#   par(mfrow = c(1,1))
#   verboseScatterplot(abs(kME[moduleGenes, column]),
#                      abs(geneTraitSignificance[moduleGenes, 1]),
#                      xlab = paste("Module Membership in", module, "module"),
#                      ylab = "Gene significance",
#                      main = paste("Module membership vs. gene significance
# 
# "),cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module)
#   abline(v=mm,lwd=2,col="red")
#   abline(h=gs,lwd=2,col="red")
#   dev.off()
# }
# 
# for (module in modNames1){
#   if(module== "grey"){ next }
#   column = match(module, modNames); # col number of interesting modules
#   pdf(paste("09.GS_MM.", module, ".pdf", sep=""),height = 6,width = 7,family='Times')
#   moduleGenes = moduleColors==module;
#   par(mfrow = c(1,1))
#   verboseScatterplot(abs(kME[moduleGenes, column]),
#                      abs(geneTraitSignificance[moduleGenes, 1]),
#                      xlab = paste("Module Membership in", module, "module"),
#                      ylab = "Gene significance",
#                      main = paste("Module membership vs. gene significance
# 
# "),cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module)
#   abline(v=mm,lwd=2,col="red")
#   abline(h=gs,lwd=2,col="red")
#   dev.off()
# }

##salmon
column = match("red", modNames);
moduleGenes = moduleColors=="red";
red_keep = abs(kME[moduleGenes, column])>mm & abs(geneTraitSignificance[moduleGenes, 1]) >gs  
red_kme<-kME[moduleGenes, ]
red_gene<-rownames(red_kme[red_keep,])
write.table(rownames(red_kme),file = "09.red_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
#write.table(salmon_gene,file = "12.salmon_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)

##brown
column = match("brown", modNames);
moduleGenes = moduleColors=="brown";
brown_keep = abs(kME[moduleGenes, column])>mm&abs(geneTraitSignificance[moduleGenes, 1]) >gs  
brown_kme<-kME[moduleGenes, ]
brown_gene<-rownames(brown_kme[brown_keep,])
write.table(rownames(brown_kme),file = "10.brown_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
#write.table(brown_gene,file = "12.brown_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)

##magenta
column = match("green", modNames);
moduleGenes = moduleColors=="green";
green_keep = abs(kME[moduleGenes, column])>mm&abs(geneTraitSignificance[moduleGenes, 1]) >gs  
green_kme<-kME[moduleGenes, ]
green_gene<-rownames(green_kme[green_keep,])
write.table(rownames(green_kme),file = "11.green_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
#write.table(magenta_gene,file = "11.magenta_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)


##black
column = match("black", modNames);
moduleGenes = moduleColors=="black";
black_keep = abs(kME[moduleGenes, column])>mm&abs(geneTraitSignificance[moduleGenes, 1]) >gs  
black_kme<-kME[moduleGenes, ]
black_gene<-rownames(black_kme[black_keep,])
write.table(rownames(black_kme),file = "12.black_all.txt",sep = "\t",quote = F,row.names = F,col.names=F)
#write.table(black_gene,file = "11.black_hub.txt",sep = "\t",quote = F,row.names = F,col.names=F)

module<-c(rownames(black_kme),rownames(red_kme),rownames(brown_kme))
write.table(module,file = "13.module.txt",sep = "\t",quote = F,row.names = F,col.names=F)
