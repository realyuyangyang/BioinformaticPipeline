rm(list = ls())
#options(stringsAsFactors = F)
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if(!dir.exists("05.ML")){dir.create("05.ML")}
setwd("05.ML")
getwd()
#
########################library#################################################
library(tidyverse)
library(dplyr)
#library(ggplot2)
#library(ggvenn)
#library(VennDiagram)
#library(ggtext)
library(magrittr)
library(tibble)
# 加载库和包
library(mlbench)
library(caret)
library(lance)
###########################Data##################################################
intersect <- read.csv('../04.PPI/Gene_PPI_Select.csv', header = T) %>% as.data.frame()
intersect <- data.frame(x=intersect$gene)
dat_expr <- read.csv('../00.rawdata/Symbol_Tpm.csv',row.names = 1,check.names = F) %>% lc.tableToNum()
dat_expr <- log2(dat_expr+1)
gene_exp <- dat_expr[intersect$x,]
gene_exp <- t(gene_exp)%>%data.frame()
group <- read.csv('../01.DEGs/colData.csv')
#
gene_exp$type <- group$group[match(rownames(gene_exp),group$sample)]
table(gene_exp$type)##统计分组情况
#control involved 
#64       58 
gene_exp$type <- factor(gene_exp$type, levels = c("CRSwNP",'Control'))##分组列转为因子类型
dat <- gene_exp
#SVM-REF
result <- data.frame()
i <- 1
for (i in c(1:10)) {
  #i <- 2
  print(i)
  set.seed(i)
  dat1 <- dat
  num <- ncol(dat1)-1
  control <- rfeControl(functions = caretFuncs,method = "cv", number = 5)
  # 执行SVM-RFE算法
  results <- rfe(dat1[,1:num],            #数据框格式
                 as.factor(dat1[,num+1]),   #因子格式
                 sizes = c(1:num),
                 rfeControl = control,
                 method = "svmRadial"
  )
  # # 结果分析
  print(results)
  svmrfe_result <- predictors(results)
  svmrfe_result
  #print(svmrfe_result)
  #if (length(table(svmrfe_result))>3) {break}
  ifelse (length(table(svmrfe_result)) == 0,
          result[1,i] <- NA,
          result[c(1:length(table(svmrfe_result))),i] <- svmrfe_result)
  colnames(result)[i] <- paste('seedd', i)
  re <- unique(t(result)) %>% as.data.frame()
  write.csv(re,file = 'svm20.csv',quote = F,row.names = T)
  
}
re2 <- read.csv('./svm20.csv',row.names = 1)
set.seed(1)
dat1 <- dat
num <- ncol(dat1)-1
control <- rfeControl(functions = caretFuncs,method = "cv", number = 5)
# 执行SVM-RFE算法
results <- rfe(dat1[,1:num],            #数据框格式
               as.factor(dat1[,num+1]),   #因子格式
               sizes = c(1:num),
               rfeControl = control,
               method = "svmRadial"
)
# # 结果分析
print(results)
svmrfe_result <- predictors(results)
svmrfe_result
#"JMJD6" "IDI1"  "SNIP1" "DNMBP"
length(svmrfe_result)
#4
svmrfe_result.write <- data.frame(symbol=svmrfe_result)
write.csv(svmrfe_result.write,"01.svmrfe_gene4.csv",row.names = F,quote = F)
# saveRDS(svmrfe_result,"../00_rawdat/00_temp/svmrfe_GSE1result.rds")
# 绘制结果
# save.image(file = '02.svmrfe.RData')
# load(file = './02.svmrfe.RData')
pdf(file = paste0("01.SVM_RFE_Accuracy.pdf"),width = 5,height = 4,family='Times')
a <- dev.cur()   #记录pdf设备
png(file = paste0("01.SVM_RFE_Accuracy.png"),width = 5, height=4, units="in", res=600,family='Times')
dev.control("enable")
par(mar = c(2,2,2,2));#下、左、上、右
plot(results, type=c("o"),
     xgap.axis = 1
)
dev.copy(which = a)  #复制来自png设备的图片到pdf
dev.off()



###################################################################################
#--------------------------------random forest--------------------------------------
rm(list = ls())
#options(stringsAsFactors = F)
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if(!dir.exists("05.ML")){dir.create("05.ML")}
setwd("05.ML")
getwd()


library(tidyverse)
library(leaps)
library(magrittr)
library(ggplot2)
library(randomForest)

#最大阈值筛选
dat_expr <- read.csv('../00.rawdata/Symbol_Tpm.csv',row.names = 1,check.names = F) %>% lc.tableToNum()
train.dat <- log2(dat_expr+1)
train.dat <- na.omit(train.dat)
hub_gene <- read.csv('../04.PPI/Gene_PPI_Select.csv', header = T) %>% as.data.frame()
group <- read.csv('../01.DEGs/colData.csv')

#control involved 


dat<-train.dat[hub_gene$gene,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('CRSwNP','Control'))

####

# colnames(dat)<-gsub('-','_',colnames(dat),fixed = T)
#dat <- gsub(dat,colnames(dat))
# 随机森林
set.seed(1)  #14
res.rf <- randomForest(group ~ ., data = dat,ntree=1000,mtry = 2,importance = T,fixed = T )
#ntree指定随机森林所包含的决策树数目
#mtry指定节点中用于二叉树的变量个数，默认情况下数据集变量个数的二次方根（分类模型）或三分之一（预测模型）
plot(res.rf, main = NULL)
png("01.RF.ntree_adj.png", 
    width = 8, height = 4, 
    bg = "white", 
    units = "in", 
    res = 300, 
    family = "Times")
plot(res.rf, main = NULL)
dev.off()
pdf("01.RF.ntree_adj.pdf", 
    width = 8, 
    height = 4, 
    family = "Times")
plot(res.rf, main = NULL)
dev.off()

##混淆矩阵----
ntree = which.min(res.rf$err.rate[,1])
res.rf = randomForest::randomForest(group ~ ., data = dat, ntree = ntree)
df.rf = data.frame(sample = rownames(res.rf$votes), vote = res.rf$votes[,2])
df.rf$group = group$group[match(df.rf$sample, group$sample)]
data <- res.rf$confusion[2:1,2:1] %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "true_label")
colnames(data)[2:3] <- c("pred_patient", "pred_normal")
data <- tidyr::gather(data, pred_label, value, -1)

#变量重要性评分
importance(res.rf,type = 1)
#type可以是1，也可以是2，用于判别计算变量重要性的方法，1表示使用精度平均较少值作为度量标准；2表示采用节点不纯度的平均减少值最为度量标准。值越大说明变量的重要性越强；
#重要性绘图
pdf(file = '02.importance_adj.pdf',w=5,h=6)
varImpPlot(res.rf,main = 'importance')
dev.off()
png(file = '02.importance_adj.png',w=400,h=500)
varImpPlot(res.rf,main = 'importance')
dev.off()

result <- data.frame(res.rf$importance)
write.csv(result,file = '03.rf_importance_adj.csv',row.names = T,quote = F)

top <- result%>%rownames_to_column(var = 'symbol')
top<-data.frame(top[order(-top$MeanDecreaseGini),])
top1<-top[1:5,]
write.csv(top1,file = '04.rf_gene_adj.csv',row.names = F,quote = F)
top1
#symbol MeanDecreaseGini
#27   HMOX1         3.844266
#17 PIK3C2G         2.924283
#44    GNMT         2.704289
#5     CYBB         2.267485
#50 SLC22A3         2.178192
#############################################################################
#---------------------------基因交集---------------------------------------
# # 
# # 
# # 
# #交集-----------
rm(list = ls())
#options(stringsAsFactors = F)
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if(!dir.exists("05.ML")){dir.create("05.ML")}
setwd("05.ML")
getwd()
# 
library(ggvenn)
# XGBoost <- read.csv('03.XGBoost.gene.csv')
# `SVM-RFE` <- read.csv('02.svmrfe_gene.csv')
# DEERG <- data.frame(symbol=intersect(XGBoost$symbol,`SVM-RFE`$symbol))
# #44
# write.table(DEERG,file = '03.hub_gene7.csv',row.names = F,quote = F)
# 
LASSO <- read.csv("01.lasso.gene.csv")
RF <- read.csv("04.rf_gene_adj.csv")
LASSO_EF_Interset <- data.frame(symbol=intersect(LASSO$symbol,RF$symbol))
# #44
write.table(LASSO_EF_Interset,file = '05.hub_gene.csv',row.names = F,quote = F)
mydata<-list('LASSO'=LASSO$symbol,'RF'=RF$symbol)
png('05.LASSO_EF_Interset_venn.png',w=5,h=5, units = 'in',res = 600,family='Times')
ggvenn(mydata,c('LASSO','RF'),
        fill_color = c("#efbe3e","#7952ab"),
        show_percentage = T,
        stroke_alpha = 0.5,
        stroke_size = 1,
        text_size = 4,
        stroke_color="#554225",
        stroke_linetype="dotdash",
        set_name_color=c("#efbe3e","#7952ab"),
        text_color = 'black')
 dev.off()
# 
# pdf('03.venn.pdf',w=5,h=5,family='Times')
# ggvenn(mydata,c('LASS0','SVM-RFE'),
#        fill_color = c("#efbe3e","#7952ab"),
#        show_percentage = T,
#        stroke_alpha = 0.5,
#        stroke_size = 1,
#        text_size = 4,
#        stroke_color="#554225",
#        stroke_linetype="dotdash",
#        set_name_color=c("#efbe3e","#7952ab"),
#        text_color = 'black')
# dev.off()


