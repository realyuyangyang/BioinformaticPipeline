#DEGs
rm(list=ls())
getwd()
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/") 
if (!dir.exists("01_DegAnalysis")) {dir.create("01_DegAnalysis")}
setwd("01_DegAnalysis")
#####################Library##########################
library(ggplot2)
library(ggrepel)
library(data.table)
library(ComplexHeatmap)
#####################DATA###########################
Expr <- read.csv("../01.DEGs/SYMBOL_Expre_data.csv", row.names = 1)
group <- read.csv("../01.DEGs/colData.csv", row.names = 1)
Res<-read.csv('../01.DEGs/DEG_all.csv',row.names=1)
############################################################
#------------------------Volcano----------------------------
log2FoldChange<-1
Res$sig[Res$padj>0.05  | abs(Res$log2FoldChange) <log2FoldChange] <- "Not"    
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


DEG$sig <- factor(DEG$sig, levels =  c("Up", "Not", "Down"))
max_lfc <- max(DEG$log2FoldChange)
min_lfc <- min(DEG$log2FoldChange)
max_padj <- round(max(-log10(DEG$padj)))
bk <- seq(0, (max_padj), (max_padj %/% 5))
bk <- append(bk, -log10(0.05), 1)



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
       title = "CRSwNP vs Control",
       subtitle = paste(sprintf('padj: %.2f;', 0.05),
                        sprintf('log2FoldChange: %.2f;', log2FoldChange),
                        sprintf('Up: %1.0f; Down: %1.0f;', dim(up_DEG)[1], dim(down_DEG)[1]),
                        sprintf('Total: %1.0f', dim(dif)[1])))
ggsave('01.volcano.png',width=8,height=6, dpi = 600)
ggsave('01.volcano.pdf',width=7,height=7)

#######################################################################
#----------------------Heatmap---------------------------------------
condition<-group
expr<-subset(Expr,select=rownames(condition))
data_repel <- rbind(up_DEG[1:20,], down_DEG[1:20,])
gene<-data_repel
select.expr<-expr[rownames(gene),]


mat<- t(scale(t(select.expr)))#归一化
mat[mat < -1] <- -1
mat[mat > 1] <- 1
annotation_col=data.frame(row.names=rownames(condition), Condition=condition$group)
annotation_row<-data.frame(row.names=rownames(data_repel),sig=data_repel$sig)


##  https://www.jianshu.com/p/671ea93108e4
pdf('05.heat.pdf',w=8,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Condition= annotation_col$Condition, col = list(Condition = c("CRSwNP" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 5),show_column_names = F,name = "mat", height = unit(6, "cm"),col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()

png('05.heat.png',w=8,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Condition= annotation_col$Condition, col = list(Condition = c("CRSwNP" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 5),show_column_names = F,name = "mat", height = unit(6, "cm"),col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()
