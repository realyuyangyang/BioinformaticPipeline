rm(list = ls())
setwd('/data/nas1/yuyangyang/Projects/YQKM-10411-2/')
if (! dir.exists('07.ROC_DCA')){
  dir.create('07.ROC_DCA')
}
setwd('07.ROC_DCA')
##################################################################
expr <- read.csv(file = "../01.DEGs/SYMBOL_Expre_data.csv", check.names = F, row.names = 1)
expr <- log2(expr+1)
hubgene <- read.csv('../05.ML/05.hub_gene.csv')
group <- read.csv('../01.DEGs/colData.csv')
group$y<-ifelse(group$group=='CRSwNP',1,0)

hub.expe <- t(expr[hubgene$symbol,]) %>% as.data.frame()

d <- merge(hub.expe, group, by.x = "row.names", by.y = "sample")


# https://rpubs.com/clayford/nomogram
library(rms)
ddist <- datadist(d)
options(datadist='ddist')
hubgene
f <- lrm(y~ ALOX5+HMOX1+PLA2G7,
         #y~ HAAO+MDK+SAMD11+CERS1+TNFRSF4,
         data = d, x = TRUE, y = TRUE)
f
summary(f)
#write.csv(f$coefficients, file = "01.lrm_coefficients.csv")

#列线图-----------------
nom.cox <- nomogram(f,  ##最佳模型
                    fun=plogis, #进行logit转换
                    funlabel="Risk of CRSwNP",
                    lp=F , ##是否显示线性预测值
                    fun.at=c(0.1,seq(0.2,0.8,by=0.1),0.9)  ##概率坐标轴的范围和刻度
)
plot(nom.cox, cex.axis  = 1.2, cex.var = 1.7)

png(filename = "01.nomogram_line_points.png", height = 5, width = 10,units='in',res=600, family = "Times")
par(family = "Times",font=2)
plot(nom.cox, cex.axis  = 0.8 ,cex.var = 1.5)
dev.off()
pdf(file = "01.nomogram_line_points.pdf", height = 5, width = 10, family = "Times")
par(family = "Times",font=2)
plot(nom.cox, cex.axis  = 0.8, cex.var = 1.5)
dev.off()

library(nomogramEx)
nomogramEx(nomo=nom.cox,np=1,digit=9)
# CA3 <- 8
# CA3.points = 11.111111111 * CA3 + -22.222222222
# 
# PPP1R15B <- 5.5
# PPP1R15B.points =  -10.371319642 * PPP1R15B + 98.527536601
# 
# MAPT <- 7
# MAPT.points =  8.713578276 * MAPT + -17.427156551
# 
# MMP9 <- 4
# MMP9.points =  7.864167562 * MMP9 + -15.728335125
# 
# ECT2 <- 6
# ECT2.points =  4.78946376 * ECT2 + -14.368391279
# 
# total points = CA3.points + PPP1R15B.points + MAPT.points + MMP9.points + ECT2.points
# Risk of COPD = 0.01111025 * total points + -1.237351876

# https://blog.csdn.net/sinat_35779431/article/details/105663759
# 不能直接用lrm函数结果
# m <- glm(f$y ~ f$x, family = binomial)
# hl <- ResourceSelection::hoslem.test(m$y, fitted(m), g=10)
# hl
###校准曲线-----------------------
library(ResourceSelection)
# https://zhuanlan.zhihu.com/p/417677244
hl1 <- ResourceSelection::hoslem.test(d$y,predict(f,d),g=10)
hl1

hl2 <- stats::resid(f,"gof")
hl2
pval <- signif(hl2[5], 3)
pval
#0.514 
set.seed(3)
cal1 <- calibrate(f, cmethod='KM', method='boot', B=1000)

pdf("02.calibrate.pdf",width=9,height=6, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.4, cex.axis=1.4, cex.main=1.8, cex.sub=1.4, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability", 
     ylab="Actual Probability", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), cex=1.2,
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
# text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "),cex=1.0)
text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))),cex=1.0)
dev.off()
png("02.calibrate.png", width=9,height=6, units = "in", res = 600, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.4, cex.axis=1.4, cex.main=1.8, cex.sub=1.4, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability", 
     ylab="Actual Probability", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), cex=1.2,
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
# text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "),cex=1.0)
text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))),cex=1.0)
dev.off()

#ROC曲线----------------------------
predicted<-predict(f,newdata = d)
lrm_predict<-ifelse(predicted >0.5, "CRSwNP", "Control") %>% as.factor()
library(PRROC)
roc_curve = roc.curve(lrm_predict, weights.class0 =  d$type == "CRSwNP", curve = T)
roc_curve$auc <- round(roc_curve$auc, digits = 2)
png(file='02.roc.png', width = 8, height = 8, res = 600, units = "in", bg = "white",family='Times')
par(pin = c(4,4), mar = c(6,6,6,1)) 
plot(roc_curve, auc.main = T, legend = F, color = 'darkblue', xlab = "1-Specificity", asp = 1 ,cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=2.0,   ##标题的缩放倍数
     main='',
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
abline(0,1); dev.off()

pdf(file='02.roc.pdf', width = 8, height = 8, family='Times')
par(pin = c(4,4), mar = c(6,6,6,1)) 
plot(roc_curve, auc.main = T, legend = F, color = 'darkblue', xlab = "1-Specificity", asp = 1 ,cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=2.0,   ##标题的缩放倍数
     main='',
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
abline(0,1); dev.off()






# dev<-d
# dev$lrm_predict<-lrm_predict
# dev$predict.y<-ifelse(dev$lrm_predict=='COPD',1,0)
# dev$`COPD prediction nomogram`<-dev$y
# ##DCA决策曲线----------------
# # source("dca.R")   #自定义函数
# # dca(data=dev, outcome="predict.y", predictors=c("DCM prediction nomogram"),smooth=F, probability=c("TRUE")) 
# 
# ##https://www.jianshu.com/p/5f3dda814cb8
# ## 列线图模型dca曲线
# library(ggDCA)
# data  <-ggDCA:: dca(f, 
#                     model.names =c('COPD prediction nomogram'))

## https://www.jianshu.com/p/a120f3f9ad78  

# library(ggprism)
# pdf(file='04.dca.pdf', width =8, height = 7, family='Times')
# par(pin = c(4,4), mar = c(6,6,6,1)) 
# ggplot(data,linetype =T,lwd = 1.0)+ 
#   theme(legend.position="top")+
#   # scale_y_continuous(
#   #   limits = c(-0.01, 0.2),
#   #   guide = "prism_minor"
#   #   )+
#   scale_colour_prism(         
#     palette = "floral",
#     # name = "Cylinders",
#     #label = c("模型1", "ALL", "None")
#   )+
#   theme(axis.title.x = element_text(size = 22, face = "bold", family = "Times"),
#         axis.title.y = element_text(size = 22, face = "bold", family = "Times"),
#         axis.text.x = element_text(size = 17, color='black',face = "bold", family = "Times",  vjust = 1, hjust = 1),
#         axis.text.y = element_text(size = 17, face = "bold", family = "Times"),
#         legend.text = element_text(size = 17, family = "Times",face = "bold"),
#         legend.title = element_blank(),
#         #plot.margin = ggplot2::margin(t=.3,b=0,l=2,r=.5, unit = "cm"),
#         text = element_text(family = "Times"),
#         panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# dev.off()
# 
# png(file='04.dca.png', width = 8, height = 7, family='Times',units='in',res=600)
# par(pin = c(4,4), mar = c(6,6,6,1)) 
# ggplot(data,linetype =T,lwd = 1.0 )+ 
#   theme(legend.position="top")+
#   # scale_y_continuous(
#   #   limits = c(-0.01, 0.2),
#   #   guide = "prism_minor"
#   #   )+
#   scale_colour_prism(         
#     palette = "floral",
#     # name = "Cylinders",
#     #label = c("模型1", "ALL", "None")
#   )+
#   theme(axis.title.x = element_text(size = 22, face = "bold", family = "Times"),
#         axis.title.y = element_text(size = 22, face = "bold", family = "Times"),
#         axis.text.x = element_text(size = 17, face = "bold", family = "Times",  vjust = 1, hjust = 1),
#         axis.text.y = element_text(size =17, face = "bold", family = "Times"),
#         legend.text = element_text(size = 17, family = "Times",face = "bold"),
#         legend.title = element_blank(),
#         #plot.margin = ggplot2::margin(t=.3,b=0,l=2,r=.5, unit = "cm"),
#         text = element_text(family = "Times"),
#         panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# dev.off()

## DCA美化--------------
library(rmda)
Nomogram<-decision_curve(y~ ALOX5+HMOX1+PLA2G7,data= d,
                         family = binomial(link ='logit'),
                         thresholds= seq(0,1, by = 0.01),
                         confidence.intervals = 0.95,
                         study.design = 'case-control',
                         population.prevalence = 0.3)


List<- list(Nomogram)
pdf('02.dca.pdf',w=10,h=9,family='Times')
par(pin = c(4,4), mar = c(6,6,6,4),font.lab=2,cex.lab=2,cex.axis=2) 
plot_decision_curve(List,
                    curve.names=c('Nomogram'),
                    cost.benefit.axis =FALSE,
                    #col= c('red','blue'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

png('02.dca.png',w=10,h=9,family='Times',units='in',res=600)
par(pin = c(4,4), mar = c(6,6,6,4),font.lab=2,cex.lab=2,cex.axis=2) 
plot_decision_curve(List,
                    curve.names=c('Nomogram'),
                    cost.benefit.axis =FALSE,
                    #col= c('red','blue'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
