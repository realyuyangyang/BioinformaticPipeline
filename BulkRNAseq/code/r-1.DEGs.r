#DEGs
rm(list=ls())
getwd()
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/") 
if (!dir.exists("01.DEGs")) {dir.create("01.DEGs")}
setwd("01.DEGs")
#######################################
library(dplyr)
library(tidyverse)
library('Biobase')
library('limma')
library('ggplot2')
library(ggplot2)
library(ggrepel)
library(DESeq2)
################################################################################
#-----------------------------------ID Convert----------------------------------
GSE136825_genecounts <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/Dataset/RawData/GSE136825_genecounts.csv", 
                                 row.names=1)
countData <- GSE136825_genecounts
library(org.Hs.eg.db)#人类
library(clusterProfiler)
data <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/Dataset/RawData/GSE136825_genecounts.csv")
data$Ensembl_ID <- substr(data$Geneid,1,15)
data$Ensembl_ID 
data
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
getwd()
#write.csv(data,"Ensembl_SYMBOL_data.csv",row.names = F)
#########################################################################################################
# 数据预处理
SYMBOL_data <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/01.DEGs/SYMBOL_data.csv", row.names=NULL)#counts
SYMBOL_data
unique_SYMBOL_data <- SYMBOL_data[!duplicated(SYMBOL_data$SYMBOL), ]
unique_SYMBOL_data
write.csv(unique_SYMBOL_data,"SYMBOL_Expre_data.csv",row.names = F)
##############################################################################
#DEGs https://yangfangs.github.io/2016/04/21/RNAseq-DEseq-analysis/
rm(list=ls())
getwd()
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/") 
if (!dir.exists("01.DEGs")) {dir.create("01.DEGs")}
setwd("01.DEGs")
####################################################################
SYMBOL_Expre_data <- read.csv("SYMBOL_Expre_data.csv", 
                              row.names=1)
countData <- SYMBOL_Expre_data
colData <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/01.DEGs/colData.csv", 
                    row.names=NULL)
#sampleNames <- c("ATNZ_1", "ATNZ_2", "ATNZ_3", "NS_1", "NS_2", "NS_3", "NS_4")
#colData = data.frame(
#  samples = c("ATNZ_1", "ATNZ_2", "ATNZ_3", "NS_1", "NS_2", "NS_3", "NS_4"),
#  group = c("AZS","AZS","AZS","Control","Control","Control","Control"))
#print(colData) # 查看 table 数据
#database <- data.frame(name=sampleNames, condition=c("AZS", "AZS", "AZS", "Control", "Control", "Control", "Control"))
#正式构建dds矩阵
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~group)
head(dds)  #查看一下构建好的矩阵
dds <- dds[rownames(counts(dds)) > 1, ] # 过滤掉在任何样本中计数都小于或等于1的基因。
dds <- estimateSizeFactors(dds) # 估计每个样本的大小因子
head(dds)

#rld <- rlog(dds)
#plotPCA(rld, intgroup=c("name","condition"))
#normalize
dds <- DESeq(dds)
#res <- results(dds, name="Group_MMF_vs_WTF") 
#res <- results(dds, contrast=c("group", "CRSwNP", "Control")) #后面的是对照 
res <- results(dds, contrast = c("group","CRSwNP","Control"))#后面的是对照 
summary(res)
head(res)##一定要查看分组
DEG <- as.data.frame(res) 
#给差异分析结果添加上下调标签
DEG$change <- as.factor(
  ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > 0.5,
         ifelse(DEG$log2FoldChange > 0.5,'Up','Down'),'Not')) 
table(DEG$change)
DEG_write <- cbind(symbol = rownames(DEG), DEG)
DEG_write <- na.omit(DEG_write)
#筛选出表达具有差异的基因 （筛选条件：|log2FoldChange| >= 0.5 & pvalue < 0.05）
sig_diff <- subset(DEG, DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > 0.5)
sig_diff_write <- cbind(symbol = rownames(sig_diff), sig_diff)
#Save
sig_diff_write <- na.omit(sig_diff_write)
write.csv(DEG_write, "DEG_all.csv", row.names = F)
write.csv(sig_diff_write, "DEG_sig.csv", row.names = F)
##############################################################################################
#-------------------------------Volcano------------------------------------------------------
rm(list=ls())
getwd()
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/") 
if (!dir.exists("01.DEGs")) {dir.create("01.DEGs")}
setwd("01.DEGs")
setwd("/data/nas1/tancj/Project/87.YQNN-10405-12/") 
if (!dir.exists("01_Degs")) {dir.create("01_Degs")}
setwd("01_Degs")
Res<-read.csv('/data/nas1/tancj/Project/87.YQNN-10405-12/01_Degs/01.DESeq2_out.csv',row.names=1)
log2FoldChange<-1
Res$sig[Res$padj>=0.05  | abs(Res$log2FoldChange) <log2FoldChange] <- "Not"    
Res$sig[Res$padj<0.05  & Res$log2FoldChange>log2FoldChange] <- "Up"     
Res$sig[Res$padj <0.05  & Res$log2FoldChange<= -log2FoldChange] <- "Down"
write.csv(Res, "01.DESeq2_out.csv",quote=F)
table(Res$sig)    #Down   Not    Up 1949 14106  3480    
## dif
dif<-Res[Res$sig!='Not',]
write.csv(dif, "02.dif_out.csv",quote=F)
## dif gene
DEGs<-rownames(dif)
write.table(DEGs,'03.DEGs.txt',sep='\t',quote=F,row.names=F,col.names=F)
DEG<-Res
DEG$Symbols <- rownames(DEG)
up_DEG <- DEG[which(DEG$sig == "Up"),]
up_DEG <- up_DEG[order(-up_DEG$log2FoldChange, up_DEG$padj),]
down_DEG <- DEG[which(DEG$sig == "Down"),]
down_DEG <- down_DEG[order(down_DEG$log2FoldChange, down_DEG$padj),]
data_repel <- rbind(up_DEG[1:10,], down_DEG[1:10,])

library(ggplot2)
library(ggrepel)
DEG$sig <- factor(DEG$sig, levels =  c("Up", "Not", "Down"))
max_lfc <- max(DEG$log2FoldChange)
min_lfc <- min(DEG$log2FoldChange)
max_padj <- round(max(-log10(DEG$padj)))
bk <- seq(0, (max_padj), (max_padj %/% 5))
bk <- append(bk, -log10(0.05), 1)

library(ggplot2)
library(ggrepel)

DEG$sig <- factor(DEG$sig, levels =  c("Up", "Not", "Down"))
max_lfc <- max(DEG$log2FoldChange)
min_lfc <- min(DEG$log2FoldChange)
max_padj   <- round(max(-log10(DEG$padj  )))
bk <- seq(0, (max_padj  ), (max_padj   %/% 5))
bk <- append(bk, -log10(0.05), 1)

ggplot(data = DEG,
       aes(x = log2FoldChange,
           y = -log10(padj  ), 
           colour = sig)) +
  scale_color_manual(values = c("#B72230", "#7b7c7d","#104680")) +
  geom_point(size = 2.5, alpha = 0.7, na.rm=T) +
  scale_x_continuous(limits = c(min_lfc*1.5, max_lfc*1.5)) +
  # scale_x_continuous(limits = c(-15, 15)) +
  scale_y_continuous(breaks = bk) +
  ylim(0,20) +
  geom_label_repel(
    data = data_repel[which(data_repel$sig=="Up"),],
    aes(label = (data_repel[which(data_repel$sig=="Up"),])$Symbols),
    #  fontface = "italic",
    size =4,
    family='Times',
    color = "#B72230",
    segment.color = "grey",
    segment.size =0.5,
    segment.alpha = 0.5,
    show.legend = F,
    direction = "y",
    hjust = 0,
    max.overlaps = 50,
    xlim = c(max_lfc + max_lfc/4, max_lfc + max_lfc/2),
    ylim = c(2, max_padj   * 0.99)
  ) +
  geom_label_repel(
    data = data_repel[which(data_repel$sig=="Down"),],
    aes(label = (data_repel[which(data_repel$sig=="Down"),])$Symbols),
    # fontface = "italic",
    size =4,
    family='Times',
    color = "#104680",
    segment.color = "grey",
    segment.size = 0.5,
    segment.alpha = 0.5,
    show.legend = F,
    direction = "y",
    hjust = 0,
    max.overlaps =50,
    xlim = c(min_lfc + min_lfc/6, min_lfc + min_lfc/2),
    ylim = c(0, max_padj   * 1.1)
  ) +
  geom_hline(yintercept = -log10(0.05),lty = 4,col = "darkgray",lwd = 0.6) +
  geom_vline(xintercept = c(-log2FoldChange, log2FoldChange),lty = 4,col = "darkgray",lwd = 0.6) +
  theme_bw(base_size = 12) +
  theme(legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=13,family='Times', face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold",color = "black",size = 24,family='Times'),
        axis.text.x = element_text(color = "black",size = 19,family='Times'),
        axis.text.y = element_text(color = "black",size = 19,family='Times'),
        axis.title.x = element_text(face = "bold",color = "black",size = 22,family='Times'),
        axis.title.y = element_text(face = "bold",color = "black",size = 22,family='Times'),
        plot.subtitle = element_text(hjust = 0.6,family = "Times", size = 15, face = "italic", colour = "black")) +
  labs(x = "log2FoldChange",y = "-log10(padj)",
       title = "LUAD vs Normal",
       subtitle = paste(sprintf('padj: %.2f;', 0.05),
                        sprintf('log2FoldChange: %.2f;', log2FoldChange),
                        sprintf('Up: %1.0f; Down: %1.0f;', dim(up_DEG)[1], dim(down_DEG)[1]),
                        sprintf('Total: %1.0f', dim(dif)[1])))
ggsave('04.volcano.png',width=7,height=7)
ggsave('04.volcano.pdf',width=7,height=7)


condition<-Group_all
expr<-subset(fpkms_all,select=rownames(condition))
gene<-data_repel
select.expr<-expr[rownames(gene),]

library(data.table)
mat<- t(scale(t(select.expr)))#归一化
mat[mat < -1] <- -1
mat[mat > 1] <- 1
annotation_col=data.frame(row.names=rownames(condition), Condition=condition$Condition)
annotation_row<-data.frame(row.names=rownames(data_repel),sig=data_repel$sig)

library(ComplexHeatmap)
##  https://www.jianshu.com/p/671ea93108e4
pdf('05.heat.pdf',w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Condition= annotation_col$Condition, col = list(Condition = c("LUAD" = "#B72230", "Normal" = "#104680"))) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),show_column_names = F,name = "mat", height = unit(6, "cm"),col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()

png('05.heat.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Condition= annotation_col$Condition, col = list(Condition = c("LUAD" = "#B72230", "Normal" = "#104680"))) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),show_column_names = F,name = "mat", height = unit(6, "cm"),col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()
##################################################################
#火山图-----
library(dplyr)
library(ggfun)
library(grid)
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(scales)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(tidyverse)
library(lance)
#------------------------------------DATA------------------------------------------------
DEG <- read.csv("DEG_all.csv", row.names = 1)
DEG<- na.omit(DEG)
sig_diff <- read.csv("DEG_sig.csv", row.names = 1)
sig_diff <- na.omit(sig_diff)
logFC_cutoff <- 0.5
DEG$symbol <- rownames(DEG)
#sig_diff <- subset(DEG,DEG$adj.P.Val < 0.01 & abs(DEG$logFC) > 2 )
dat_rep <-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2FoldChange,decreasing = F),],10))),]
dat_rep
deg2 <-rbind(deg,dat_rep) 
#指定显示上调前10名的基因，和下调前10名的基因
volcano_plot <- ggplot(data = DEG) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), 
                 color = log2FoldChange,
                 size = -log10(padj))) + 
  geom_text_repel(data =  deg2 %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(change != "NOT") %>%
                    dplyr::arrange(desc(log2FoldChange)) %>%
                    # dplyr::slice(1:10) %>%
                    dplyr::filter(change == "UP"),
                  aes(x = log2FoldChange, y = -log10(padj), label = symbol),
                  box.padding = 0.5,
                  nudge_x = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 10,
                  direction = "y",
                  hjust = "left",
                  max.overlaps = 200
  )+
  geom_text_repel(data =  deg2 %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(change != "NOT") %>%
                    dplyr::filter(change != "UP") %>%
                    dplyr::arrange(desc(-log2FoldChange)) %>%
                    # dplyr::slice(1:10) %>%
                    dplyr::filter(change == "DOWN"),
                  aes(x = log2FoldChange, y = -log10(padj), label = symbol),
                  box.padding = 0.5,
                  nudge_x = -0.2,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "left",
                  max.overlaps = 200
  ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  scale_size(range = c(1,7)) + 
  ggtitle(label = "Volcano Plot") + 
  #xlim(c(-15, 11)) + 
  ylim(c(-1, 191)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5)
  ) + 
  annotate(geom = "text", x = 4, y = 0, label = "adj.p = 0.05", size = 5) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")
    ), 
    xmin = (-logFC_cutoff)-5, 
    xmax = -logFC_cutoff,
    ymin = 190,
    ymax = 190
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "DOWN",
      gp = grid::gpar(col = "#74add1")
    ),
    xmin = (-logFC_cutoff)-5, 
    xmax = -logFC_cutoff,
    ymin = 190,
    ymax = 190
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#d73027")
    ), 
    xmin = logFC_cutoff+5, 
    xmax = logFC_cutoff,
    ymin = 190,
    ymax = 190
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "UP",
      gp = grid::gpar(col = "#d73027")
    ),
    xmin = logFC_cutoff+5, 
    xmax = logFC_cutoff,
    ymin = 190,
    ymax = 190
  ) 
volcano_plot
ggsave('03.volcano_DEG.png', volcano_plot,width = 8, height = 6)
ggsave('03.volcano_DEG.pdf', volcano_plot,width = 8, height = 6,dpi = 600)


library(ComplexHeatmap)
library(circlize)

dat_rep <-sig_diff[rownames(sig_diff)%in%
                    rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange,decreasing = T),],10),
                                   head(sig_diff[order(sig_diff$log2FoldChange,decreasing = F),],10))),]
deg <- sig_diff[rownames(dat_rep),]
dat_rep <- rbind(deg,dat_rep)
# dat_rep <- deg2
rt <- count
# group <- colData
rt <- rt[,rownames(group)]
rt <- rt[rownames(dat_rep),]
#group <- group
#x <- rt
x<-log2(rt+1)
mat <- t(scale(t(x)))#归一化
df1 <- as.data.frame(mat)
mat[mat < (-1)] <- (-1)
mat[mat > 1] <- 1

# mat1 <- as.matrix(mat)
# cir1<-t(scale(t(mat1)))
# head(cir1)
# mat1<-cir1

# pdf('04.heatmap_DEG.pdf',  w=6,h=6,family='Times')
# col_fun1 = colorRamp2(c(-2, 0, 2), c("darkgreen", "white", "orange"))
# circos.par(gap.after = c(80))#这里代码就是空出一段用于添加label
# circos.heatmap(mat1, col = col_fun1,dend.side = "inside",rownames.side = "outside",track.height = 0.4,)
# lgd = Legend(title = "Log2FC", col_fun = col_fun1)
# grid.draw(lgd)
# circos.clear()
# dev.off()
# 
# png('04.heatmap_DEG.png',w=6,h=6,units='in',res=600,family='Times')
# col_fun1 = colorRamp2(c(-2, 0, 2), c("darkgreen", "white", "orange"))
# circos.par(gap.after = c(40))#这里代码就是空出一段用于添加label
# circos.heatmap(mat1, col = col_fun1,dend.side = "inside",rownames.side = "outside",track.height = 0.4)
# lgd = Legend(title = "Log2FC", col_fun = col_fun1)
# grid.draw(lgd)
# circos.clear()
# dev.off()
library(ComplexHeatmap)
library(circlize)
# df = read.csv('../00_rawdata/11.dat.GSE58294.csv',row.names = 1)
# df <- log2(df+30)
# rt <- df
# rt <- rt[rownames(dat_rep),]
# rt <- rt[rownames(sig_diff),]
#r<-na.omit(rt)

# group <- df.group
# group <- group[order(group$group),]#分组排序
# rt<- rt[,group$sample]#样本按照分组进行排序
# r<-r[,group$sample]


# x<-rt
# x <- log2(x+5)
# mat <- t(scale(t(x)))#归一化
# df1 <- as.data.frame(mat)
# mat[mat < (-2)] <- (-2)
# mat[mat > 2] <- 2

pdf('04.heatmap.pdf',  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$Type, col = list(Group = c("Tumor" = "orange", "Normal" = "darkgreen"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("darkgreen", "white","orange"))(100))
dev.off()

png('04.heatmap.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$Type, col = list(Group = c("Tumor" = "orange", "Normal" = "darkgreen"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("darkgreen", "white","orange"))(100))
dev.off()

###############################################################################################
#########################expretest boxplot#####################################################
#######################Heatmap###############
##提取显著差异基因
## 提取FC前100
up_100<-DEG %>% as_tibble() %>% 
  mutate(genename=rownames(DEG)) %>% 
  dplyr::arrange(desc(log2FoldChange)) %>% 
  .$genename %>% .[1:100] ## 管道符中的提取
## FC低前100
down_100<-DEG %>% as_tibble() %>% 
  mutate(genename=rownames(DEG)) %>% 
  dplyr::arrange(log2FoldChange) %>% 
  .$genename %>% .[1:100] ## 管道符中的提取
index<-c(up_100,down_100)  

## 开始绘图-最简单的图
library(pheatmap)
#pheatmap(eset_dat[index,],show_colnames =F,show_rownames = F)
eset_dat <- GSE9257_expr
## 稍微调整细节
index_matrix<-t(scale(t(eset_dat[index,])))##归一化
index_matrix[index_matrix>1]=1
index_matrix[index_matrix<-1]=-1
head(index_matrix)
## 添加注释
#anno=data.frame(group=group_list)
#rownames(anno)=colnames(index_matrix)
#anno##注释信息的数据框
pheatmap(index_matrix,
         scale = "row",
         color=colorRampPalette(c("navy", "white", "red"))(50),
         show_colnames =F,
         show_rownames = F,
         cluster_cols = F)
#         annotation_col=anno)

#################################################################################
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2FoldChange,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = log2FoldChange,
                          y = -log10(pvalue), 
                          color =change)) +
  scale_color_manual(values = c( "green","darkgray","firebrick2")) +
  # scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-0.5,0.5),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log2(Fold Change)",
       y = "-log10 (pvalue)")
volcano_plot
ggsave('01.volcano.png', volcano_plot,width = 8, height = 6, dpi = 600)
ggsave('01.volcano.pdf', volcano_plot,width = 8, height = 6)

#密度热图
library(ComplexHeatmap)
library(circlize)
df = read.csv('SYMBOL_Expre_data.csv',row.names = 1)
df <- log2(df+1)
rt <- df
rt <- rt[rownames(dat_rep),]
# rt <- rt[rownames(sig_diff),]
#r<-na.omit(rt)
rt <- t(rt)%>%data.frame()
group <- read.csv("colData.csv")
rt$group <- group$group[match(rownames(rt),group$samples)]


