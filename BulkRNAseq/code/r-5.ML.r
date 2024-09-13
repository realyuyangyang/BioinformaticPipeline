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
library(ggplot2)
library(ggvenn)
library(VennDiagram)
library(ggtext)
library(dplyr)
library(ggplot2)
library(magrittr)
library(ggplot2)
library(lance)
library(magrittr)
library(tibble)
# 加载库和包
library(mlbench)
library(caret)
library(lance)
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
# #LASSO----------
library(glmnet)
result <- data.frame()
for (j in c(1:2000)) {
   set.seed(j)
   res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$type, family = "binomial",
                          type.measure = "auc")
   l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
   l.coef
   coef.min = coef(res.lasso, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
   res.lasso$lambda.min
   active.min = which(coef.min@i != 0)
   lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] %>% as.data.frame()
   print(j)
   ifelse (nrow(lasso_geneids) == 0,
           result[1,j] <- NA,
           result[c(1:nrow(lasso_geneids)),j] <- lasso_geneids)
   colnames(result)[j] <- paste('seed', j)
 }
re <- unique(t(result)) %>% as.data.frame()
write.csv(re,file = 'lasso2000.csv',quote = F,row.names = T)
##################################################################################
re <- read.csv('./lasso2000.csv',row.names = 1)
set.seed(5)
library(glmnet)
res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$type, family = "binomial",
                        type.measure = "deviance")
# 
plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
# 
pdf('01.lasso.CV.pdf', w=8,h=6,family='Times')
plot(res.lasso)
dev.off()
png('01.lasso.CV.png', w=8,h=6,family='Times',units='in',res=600)
plot(res.lasso)
dev.off()
# 
# 
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
coef.min = coef(res.lasso, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
res.lasso$lambda.min
# #0.09792358
# # 找出那些回归系数没有被惩罚为0的
# 
pdf('01.lasso.Coef.pdf', w=8,h=6,family='Times')
plot(res.lasso$glmnet.fit, xvar = 'lambda')
abline(v = log(res.lasso$lambda.min), lty = 3,lwd = 1,col = "black")
dev.off()
png('01.lasso.Coef.png', w=8,h=6,family='Times',units='in',res=600)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
abline(v = log(res.lasso$lambda.min), lty = 3,lwd = 1,col = "black")
dev.off()
# 
# 
active.min = which(coef.min@i != 0)
# # 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
# #  "C15orf48" "ANKRD9"   "GLCE"     "NMNAT1"   "HSD17B4"
lasso_geneids.write <- data.frame(symbol=lasso_geneids)
write.csv(lasso_geneids.write,file = '01.lasso.gene.csv',row.names = F,quote = F)
#svm-rfe----------
#intersect <- read.csv('../04_Venn/01.Intersection_gene5.csv') %>% as.data.frame()
#intersect <- data.frame(x=intersect$symbol)
#dat_expr<-read.csv('../00_rawdata/01.exprlog.csv',row.names = 1,check.names = F) %>% lc.tableToNum
#group <- read.csv('../00_rawdata/01.group.csv')
#colnames(group) <- c('sample', 'type')
#dat<-dat_expr[intersect$x,group$sample]%>%t%>%as.data.frame()
#dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
#table(group$type)
#dat$type<-factor(dat$type,levels = c('control','DN'))
result <- data.frame()
i <- 1
for (i in c(1:100)) {
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
dev.off()

# # 
# # 
# # 
# #交集-----------
# rm(list = ls())
# setwd("/data/nas1/liky/project/33_YQ805-1/")
# if (! dir.exists('./05_svmrfe')){
#   dir.create('./05_svmrfe')
# }
# setwd('./05_svmrfe')
# 
# library(ggvenn)
# XGBoost <- read.csv('03.XGBoost.gene.csv')
# `SVM-RFE` <- read.csv('02.svmrfe_gene.csv')
# DEERG <- data.frame(symbol=intersect(XGBoost$symbol,`SVM-RFE`$symbol))
# #44
# write.table(DEERG,file = '03.hub_gene7.csv',row.names = F,quote = F)
# 
# 
# mydata<-list('LASS0'=LASS0$symbol,'SVM-RFE'=`SVM-RFE`$symbol)
# png('03.venn.png',w=5,h=5, units = 'in',res = 600,family='Times')
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




# # # Boruta ------------------------------------------------------------------
# rm(list = ls())
# setwd("/data/nas1/liky/project/33_YQ805-1/")
# if (! dir.exists('./05_svmrfe')){
#   dir.create('./05_svmrfe')
# }
# setwd('./05_svmrfe')
# 
# intersect <- read.csv('../02_Venn_GO_KEGG/01.Intersection_gene10.csv') %>% as.data.frame()
# #intersect <- read.csv('../05_MR/01.resultsPPI0.4.csv',row.names = 1) %>% as.data.frame()
# #intersect <- data.frame(x=intersect$SYMBOL)
# dat_expr <- read.csv('../00_rawdata/01.exprlog.csv',row.names = 1,check.names = F) %>% lc.tableToNum
# dat_group <- read.csv('../00_rawdata/01.group.csv')
# 
# 
# dat<-dat_expr[intersect$symbol, dat_group$sample]%>%t%>%as.data.frame()
# dat<-merge(dat,dat_group,by.x='row.names',by.y='sample')#%>%tibble::column_to_rownames(var = 'Row.names')
# rownames(dat) <- dat$Row.names
# dat<-dat[,-1]
# 
# df.exp <- dat
# df.exp$type<-factor(df.exp$type)
# 
# library("Boruta")
# set.seed(44242424)
# res.Boruta<-Boruta(x=df.exp[,1: ncol(df.exp)-1], y=as.factor(df.exp[,ncol(df.exp)]), pValue=0.01, mcAdj=T,
#                    maxRuns=300)
# Boruta<-attStats(res.Boruta) #给出Boruta算法的结果
# table(res.Boruta$finalDecision)
# boruta_geneids<-Boruta[Boruta$decision=='Confirmed',]%>%rownames(.)
# # Tentative Confirmed  Rejected
# # 1         6         3
# boruta_geneids
# #"GSTM4"  "ABHD10" "PPP2CB" "FDXR"   "NUP160" "BGN"
# write.csv(boruta_geneids,"03.boruta_gene.csv",row.names = F,quote = F)
# 
# ##定义一个函数提取每个变量对应的重要性值。
# 
# library(dplyr)
# boruta.imp <- function(x){
#   imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
#   colnames(imp) <- c("Variable","Importance")
#   imp <- imp[is.finite(imp$Importance),]
#   
#   variableGrp <- data.frame(Variable=names(x$finalDecision),
#                             finalDecision=x$finalDecision)
#   
#   showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
#                         finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
#   
#   variableGrp <- rbind(variableGrp, showGrp)
#   
#   boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
#   
#   sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
#     summarise(median=median(Importance)) %>% arrange(median)
#   sortedVariable <- as.vector(sortedVariable$Variable)
#   
#   
#   boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
#   
#   invisible(boruta.variable.imp)
# }
# 
# boruta.variable.imp <- boruta.imp(res.Boruta)
# 
# head(boruta.variable.imp)
# #devtools::install_github("Tong-Chen/YSX")
# library(YSX)
# library(ImageGP)
# 
# sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
#            legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
#            xtics_angle = 90)
# 
# pdf("03.Boruta.pdf",w = 7, h = 6,family='Times')
# sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
#            legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
#            xtics_angle = 90)
# dev.off()
# 
# png("03.Boruta.png",w = 7, h = 6,family='Times',units='in',res=600)
# sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
#            legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
#            xtics_angle = 90)
# dev.off()
# # # 
# # # 
# # # 
# xgboost -----------------------------------------------------------------
# rm(list = ls())
# setwd("/data/nas1/liky/project/33_YQ805-1/")
# if (! dir.exists('./05_svmrfe')){
#   dir.create('./05_svmrfe')
# }
# setwd('./05_svmrfe')
# 
# intersect <- read.csv('../02_Venn_GO_KEGG/01.Intersection_gene5.csv') %>% as.data.frame()
# #intersect <- read.csv('../05_MR/01.resultsPPI0.4.csv',row.names = 1) %>% as.data.frame()
# #intersect <- data.frame(x=intersect$SYMBOL)
# dat_expr <- read.csv('../00_rawdata/01.exprlog.csv',row.names = 1,check.names = F) %>% lc.tableToNum
# dat_group <- read.csv('../00_rawdata/01.group.csv')
# colnames(dat_group) <- c('sample', 'type')
# 
# dat<-dat_expr[intersect$symbol, dat_group$sample]%>%t%>%as.data.frame()
# dat<-merge(dat,dat_group,by.x='row.names',by.y='sample')#%>%tibble::column_to_rownames(var = 'Row.names')
# rownames(dat) <- dat$Row.names
# dat<-dat[,-1]
# 
# df.exp <- dat
# df.exp$type<-factor(df.exp$type)
# library('xgboost')
# library("Matrix")
# library('PRROC')
# library('ggplot2')
# set.seed(1765727)
# ### XGB模型
# colnames(df.exp) <- gsub("-", "_", colnames(df.exp))
# train_matrix <- sparse.model.matrix(type ~.-1, data = df.exp)
# # 训练集的数据预处理
# # 将trainset的1-8列（自变量）转换为矩阵
# traindata1<- data.matrix(df.exp[,-ncol(df.exp)])
# # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
# traindata2 <- Matrix(traindata1,sparse = T)
# # 将因变量转换为numeric类型，-1是为了从0开始计数
# train_y <- as.numeric(df.exp[,ncol(df.exp)])-1
# # 将自变量和因变量拼接为list
# traindata <- list(data=traindata2,label=train_y)
# dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label)
# set.seed(989655)
# res.xgb <- xgboost(data = dtrain,max_depth=2, eta=0.3,
#                    objective='binary:logistic', nround=25)
# # eta：range: [0,1]，参数值越大，越可能无法收敛。把学习率 eta 设置的小一些，小学习率可以使得后面的学习更加仔细。
# # max_depth：每颗树的最大深度，树高越深，越容易过拟合
# # nround：迭代次数
# # objective：定义最小化损失函数类型，常用参数binary:logistic,multi:softmax,multi:softprob
# 
# xgb_importance <- xgb.importance(train_matrix@Dimnames[[2]], model = res.xgb)    ##特征重要度
# xgb_importance $Feature<-gsub('.','-',xgb_importance $Feature,fixed = T )
# 
# ggplot(xgb_importance, aes(x= reorder( Feature,Gain), y=Gain,fill='blue')) +
#   geom_bar(stat="identity") +
#   theme_classic() +
#   guides(fill=FALSE)+
#   #theme(legend.position = )+
#   #geom_hline(yintercept = mean(xgb_importance$Gain),lty = 4,col = "darkred",lwd = 0.8) +
#   coord_flip()+
#   theme_bw()+
#   ggtitle('XGBoost')+
#   theme(plot.title = element_text(size=24,color='black', face = "bold",family='Times'),
#         axis.title.x =element_text(size=18,color='black', face = "bold",family='Times'),
#         axis.text.x =element_text(size=16, color='black', face = "bold",family='Times'),
#         axis.title.y =element_blank(),
#         axis.text.y=element_text(size=16,   color='black',face = "bold",family='Times'),
#         legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
#         legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
#         title=element_text(size=20, color='black', face = "bold",family='Times'),
#         strip.text = element_text(size = 14,family = "Times", face = "bold"))+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
#   labs(x="gene",y="Gain",fill="")
# ggsave('03.XGBoost_importance.pdf',w=8,h=8)
# ggsave('03.XGBoost_importance.png',w=8,h=8)
# xgb.write <- data.frame(symbol=xgb_importance$Feature)
# xgb.write$symbol
# #[1] "TRAF4"    "SLC2A3"   "IDI1"     "TP53INP2" "MAFF"     "JOSD1"    "DNMBP"    "RBM38"    "SNIP1"   
# write.csv(xgb.write,'03.XGBoost.gene.csv',quote=F,row.names=F)
# 
# #交集-------------
# rm(list = ls())
# setwd("/data/nas1/liky/project/28_TY-0119-4/")
# if (! dir.exists('./05_svmrfe')){
#   dir.create('./05_svmrfe')
# }
# setwd('./05_svmrfe')
# 
# Lasso <- read.csv('./01.lasso.gene.csv')$symbol
# `SVM-RFE` <- read.csv('./02.svmrfe_gene.csv')$symbol
# xgb <- read.csv('./03.XGBoost.gene.csv')$symbol
# 
# DEERG <- data.frame(symbol=intersect(Lasso,intersect(`SVM-RFE`,xgb)))
# DEERG$symbol
# #"ZC3H12A" "NUAK1"   "LENG8"
# # modlegene <- read.csv('../02_WGCNA/06.modlegene22.csv')
# # intersect <- data.frame(symbol = intersect(DEERG$symbol, modlegene$.))
# 
# write.table(DEERG,file = '04.hubgene.csv',row.names = F,quote = F)
# 
# 
# library(ggvenn)
# mydata<-list('Lasso'=Lasso,'SVM-RFE'=`SVM-RFE`,'XGBoost'=xgb)
# pdf('04.hubgene_venn.pdf',w=5,h=5,family='Times')
# ggvenn(mydata,c('Lasso','SVM-RFE','XGBoost'),
#        fill_color = c("#7aa4dd","#ceadd0",'#b2e7cb'),
#        show_percentage = T,
#        stroke_alpha = 0.5,
#        stroke_size = 0.7,
#        text_size = 4,
#        stroke_color="white",
#        stroke_linetype="dotdash",
#        set_name_color=c("#7aa4dd","#ceadd0",'#b2e7cb'),
#        text_color = 'black')
# dev.off()
# png('04.hubgene_venn.png',w=5,h=5, units = 'in',res = 600,family='Times')
# ggvenn(mydata,c('Lasso','SVM-RFE','XGBoost'),
#        fill_color = c("#7aa4dd","#ceadd0",'#b2e7cb'),
#        show_percentage = T,
#        stroke_alpha = 0.5,
#        stroke_size = 0.7,
#        text_size = 4,
#        stroke_color="white",
#        stroke_linetype="dotdash",
#        set_name_color=c("#7aa4dd","#ceadd0",'#b2e7cb'),
#        text_color = 'black')
# dev.off()


# # # 
# library(ggVennDiagram)
# library(VennDiagram)
# venn_list<-list(`XGBoost`=xgb,`Lasso`=Lasso,`SVM-RFE`=`SVM-RFE`)
# #venn_list<-list(`Boruta`=boruta_geneids,`Xgboost`=xgb_importance$Feature)
# venn.plot <-venn.diagram(venn_list,  filename = NULL ,
#                          force.unique = T,
#                          print.mode = c("percent", "raw"),
#                          fill = c('skyblue','sandybrown','palegreen'), alpha = 0.8,
#                          col = 'grey',
#                          cex = 1.4,
#                          cat.fontface = "bold",
#                          cat.col = c('skyblue','sandybrown','palegreen'),
#                          #cat.pos = c(0,0),
#                          cat.cex = 1.5,
#                          margin = 0.1,
#                          scaled =FALSE
# )
# inter <- get.venn.partitions(venn_list)
# #将venn.plot通过grid.draw画到pdf文件中
# pdf("04.hubgene_venn.pdf")
# grid.draw(venn.plot)
# dev.off()
# 
# png("04.hubgene_venn.png")
# grid.draw(venn.plot)
# dev.off()
# 
# venn.plot <-venn.diagram(venn_list,  filename = NULL ,
#                          force.unique = T,
#                          print.mode = c("percent", "raw"),
#                          fill = c('skyblue','sandybrown','palegreen'), alpha = 0.8,
#                          col = 'grey',
#                          cex = 1.4,
#                          cat.fontface = "bold",
#                          cat.col = c('skyblue','sandybrown','palegreen'),
#                          #cat.pos = c(0,0),
#                          cat.cex = 1.5,
#                          margin = 0.1,
#                          scaled =FALSE
# )
# inter <- get.venn.partitions(venn_list)
# for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
# model_gene<-inter[1,]$..values..
# model_gene[["1"]]
# model_gene.write <- data.frame(symbol=model_gene[["1"]])
# write.csv(model_gene.write,'04.hubgene_venn.csv', row.names = F, quote = F)
# # 
# # 
# # rm(list = ls())
# # setwd('/data/nas1/liky/project/11_YQKM-30101-9/')
# # if (! dir.exists('./06_lasso_svmrfe_boruta')){
# #   dir.create('./06_lasso_svmrfe_boruta')
# # }
# # setwd('./06_lasso_svmrfe_boruta')
# # 
# # library(tidyverse)
# # library(leaps)
# # 
# # # 1.导入数据--BSV----
# # train.dat <- read.csv('../00_rawdata/01.count_vsd_GSE89408.csv',check.names = F,row.names = 1)
# # hub_gene <- read.csv('../05_MR/03.MRgene.csv',check.names = F,row.names = 1)
# # hub_gene <- read.csv('../05_MR/00.results.csv',row.names = 1) %>% as.data.frame()
# # hub_gene <- data.frame(x=hub_gene$SYMBOL)
# # 
# # group <- read.csv('../00_rawdata/01.group_GSE89408.csv',check.names = F)
# # colnames(group) <- c('sample','group')
# # #lasso_gene <- read.csv('../04_LASSO_GSE103552/03.lasso.gene.csv',check.names = F,row.names = 1)
# # dat<-train.dat[hub_gene$x,group$sample]%>%t%>%as.data.frame()
# # dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
# # dat <- dat[order(dat$group),]
# # table(group$group)
# # # Normal  Tumor 
# # # 32     18
# # 
# # da.scale <- dat %>%
# #   mutate(across(.cols = !contains("group"),.fns = scale), # 对数值数据进行标准化
# #          group = factor(group, levels = c("control","OA")) # 将因子变量设置为实验设计顺序。
# #          #若要将分类变量用作线性回归分析，一定要进行此设置，将对照组排在第一位，分析时其会作为哑变量，即，它在回归方程中的回归系数会为0。
# #   ) %>% 
# #   data.frame()
# # 
# # # 2.最优子集法进行特征选择------
# # sub.fit <- regsubsets(group ~ ., data = da.scale)
# # sub.fit %>%
# #   summary()
# # 
# # best_summary <- summary(sub.fit)
# # names(best_summary) # 查看列表数据中包含的对象
# # which.min(best_summary$bic)
# # 
# # # 使用赤池信息量准则，马洛斯的Cp，贝叶斯信息量准则和修正R方四种方法进行特征选择
# # ## leaps不输出AIC值，但AIC与Cp成正比。
# # pdf("01.model_select.pdf",width=unit(8,"cm"),height=unit(4,"cm"))
# # par(mfrow=c(1,2), oma = c(0,0,0,0))
# # plot(best_summary$adjr2,ylab = "adjr2" ,xlab = "number of features")
# # plot(sub.fit,scale = "adjr2")
# # dev.off()
# # 
# # png("01.model_select.png",w=8,h=4,units = "in", res = 600)
# # par(mfrow=c(1,2), oma = c(0,0,0,0))
# # plot(best_summary$bic,ylab = "bic" ,xlab = "number of features")
# # plot(sub.fit,scale = "bic")
# # dev.off()
# # 
# # #当自变量为4时，bic最小
# # BSR_gene <- data.frame(symbol=c('OLR1','PCNA','KPNA2','CDKN3'))
# # write.csv(BSR_gene,file = '02.hub_gene.csv',row.names = T)
# # 
# # 
# # 
# # 
# #非自限性--------------
# rm(list = ls())
# setwd("/data/nas1/liky/project/33_YQ805-1/")
# if (! dir.exists('./05_svmrfe')){
#   dir.create('./05_svmrfe')
# }
# setwd('./05_svmrfe')
# 
# train.dat <- read.csv('../00_rawdata/01.exprlog.csv',check.names = F,row.names = 1)
# hub_gene <- read.csv('../05_MR/03.MRgene.csv',check.names = F,row.names = 1)
# hub_gene <- read.csv('../02_Venn_GO_KEGG/01.Intersection_gene15.csv',row.names = 1) %>% as.data.frame()
# hub_gene <- data.frame(x=hub_gene$symbol)
# 
# group <- read.csv('../00_rawdata/01.group.csv',check.names = F)
# colnames(group) <- c('sample','Group')
# #lasso_gene <- read.csv('../04_LASSO_GSE103552/03.lasso.gene.csv',check.names = F,row.names = 1)
# data<-train.dat[hub_gene$x,group$sample]%>%t%>%as.data.frame()
# data<-merge(data,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
# data <- data[order(data$Group),]
# data$Group = ifelse(data$Group=='DN',1,0)
# table(group$Group)
# 
# library(tidyverse)
# library(caret)
# library(DALEX)
# 
# ##  https://zhuanlan.zhihu.com/p/354108603
# #set.seed(44)
# set.seed(1002342)
# ozone_rf<-train(Group~.,data = data,
#                 method = "rf",
#                 ntree = 400
# )
# ozone_glm<-train(Group~.,data = data,
#                  method = "glm")
# ozone_gbm<-train(Group~.,data =  data,
#                  method = "gbm"
# )
# 
# #模型解释
# explainer_rf<-explain(ozone_rf,label = "rf",
#                       data = data,
#                       y = data$Group)
# 
# explainer_glm<-explain(ozone_glm,label = "glm",
#                        data = data,
#                        y = data$Group)
# explainer_gbm<-explain(ozone_gbm,label = "gbm",
#                        data = data,
#                        y = data$Group)
# per_rf<-model_performance(explainer_rf)
# per_glm<-model_performance(explainer_glm)
# per_gbm<-model_performance(explainer_gbm)
# 
# #累积残差分布图
# png("01.cumulative.png", width = 6.5, height = 6, bg = "white", units = "in", res = 600,family='Times')
# par(pin = c(4,4), mar = c(6,6,6,1))
# plot(per_rf)
# dev.off()
# 
# pdf("01.cumulative.pdf", width = 6.5, height = 6,family='Times')
# par(pin = c(4,4), mar = c(6,6,6,1))
# plot(per_rf,per_glm,per_gbm)
# dev.off()
# 
# #累积残差分布图
# png("02.boxplot.png", width = 6.5, height = 6, bg = "white", units = "in", res = 600,family='Times')
# plot(per_rf,per_glm,geom = "boxplot")
# dev.off()
# 
# pdf("02.boxplot.pdf", width = 6.5, height = 6,family='Times')
# plot(per_rf,per_glm,per_gbm,geom = "boxplot")
# dev.off()
# 
# 
# #变量重要性分析
# # 分析在不同的模型中，不同变量对于模型预测的相对重要性程度。
# #
# # 此处损失函数为均方根误差，解释为缺了该变量会对响应变量的预测值带来多大程度的影响。
# importance_rf<-variable_importance(
#   explainer_rf,
#   loss_function = loss_root_mean_square
# )
# importance_glm<-variable_importance(
#   explainer_glm,
#   loss_function = loss_root_mean_square
# )
# importance_gbm<-variable_importance(
#   explainer_gbm,
#   loss_function = loss_root_mean_square
# )
# 
# plot(importance_glm,importance_rf)
# 
# pdf('03.Feature_Importance.pdf',w=6.5,h=7,family='Times')
# plot(importance_glm,importance_rf,importance_gbm)
# dev.off()
# 
# png('03.Feature_Importance.png',w=6.5,h=7,bg = "white", units = "in", res = 600,family='Times')
# plot(importance_glm,importance_rf,importance_gbm)
# dev.off()
# 
# importance_rf1 <- data.frame(importance_rf)
# gene<-c('TNF','SERPINA3','F12','PLAU','CSF3')
# write.table(gene,'06.gene.txt',sep='\t',quote=F,row.names=F,col.names=F)
# 
