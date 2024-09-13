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
set.seed(10)  #14
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
top1<-top[1:11,]
top1
write.csv(top1,file = '04.rf_gene_adj.csv',row.names = F,quote = F)
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