####-DAY10 VS CTL-###
#Note: 
#1. To change countData and colData.
#2. To rename figures at line 161, 162, 173, 178(optional).
rm(list=ls())
getwd()
setwd("../")
if (!dir.exists("01.DAY10vsCTL")) {dir.create("01.DAY10vsCTL")}
setwd("01.DAY10vsCTL")
getwd()
#library
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(ComplexHeatmap)
#DATA
countData_SYMBOL <- read.csv("../01.DEGs/countData_SYMBOL.csv",row.names = 1)
countData <- countData_SYMBOL[,c(1:4, 9:12)]
countData
##group#
colData = data.frame(
  samples = c("CTL_r1", "CTL_r2", "CTL_r3", "CTL_r4", "DAY10_r1" ,"DAY10_r2" ,"DAY10_r3", "DAY10_r4"),
  group = c("Control","Control","Control","Control", "Treat","Treat","Treat","Treat"))
print(colData) # 查看 table 数据
#正式构建dds矩阵
dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~group)
head(dds)  #查看一下构建好的矩阵
dds <- dds[rownames(counts(dds)) > 1, ] # 过滤掉在任何样本中计数都小于或等于1的基因。
dds <- estimateSizeFactors(dds) # 估计每个样本的大小因子
head(dds)
#
#rld <- rlog(dds)
#plotPCA(rld, intgroup=c("name","condition"))
#normalize
dds <- DESeq(dds)
#res <- results(dds, name="Group_MMF_vs_WTF") 
#res <- results(dds, contrast=c("group", "CRSwNP", "Control")) #后面的是对照 
res <- results(dds, contrast = c("group","Treat","Control"))#后面的是对照 
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
#write.csv(DEG_write, "DAY1vsCTL_DEG_all.csv", row.names = F)
#write.csv(sig_diff_write, "DAY1vsCTL_DEG_sig.csv", row.names = F)
################Plot#####################################
#-------------------------------Volcano------------------------------------------------------
#Res<-read.csv('../DAY1vsCTL_DEG_all.csv',row.names=1)
Res <- DEG_write
log2FoldChange<-1
Res$sig[Res$padj>=0.05  | abs(Res$log2FoldChange) <log2FoldChange] <- "Not"    
Res$sig[Res$padj<0.05  & Res$log2FoldChange>log2FoldChange] <- "Up"     
Res$sig[Res$padj <0.05  & Res$log2FoldChange<= -log2FoldChange] <- "Down"
#write.csv(Res, "01.DESeq2_DAY1vsCTL_out.csv",quote=F)
table(Res$sig)    #Down   Not    Up 1949 14106  3480    
## dif
dif<-Res[Res$sig!='Not',]
#write.csv(dif, "01.DESeq2_DAY1vsCTL_out.csv",quote=F)
## dif gene
DEGs<-rownames(dif)
#write.table(DEGs,'01.DAY1vsCTL_DEGs.txt',sep='\t',quote=F,row.names=F,col.names=F)
DEG<-Res
DEG$Symbols <- rownames(DEG)
up_DEG <- DEG[which(DEG$sig == "Up"),]
up_DEG <- up_DEG[order(-up_DEG$log2FoldChange, up_DEG$padj),]
down_DEG <- DEG[which(DEG$sig == "Down"),]
down_DEG <- down_DEG[order(down_DEG$log2FoldChange, down_DEG$padj),]
data_repel <- rbind(up_DEG[1:10,], down_DEG[1:10,])
#
DEG$sig <- factor(DEG$sig, levels =  c("Up", "Not", "Down"))
max_lfc <- max(DEG$log2FoldChange)
min_lfc <- min(DEG$log2FoldChange)
max_padj   <- round(max(-log10(DEG$padj+1  )))
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
       title = "Treat vs Control",
       subtitle = paste(sprintf('padj: %.2f;', 0.05),
                        sprintf('log2FoldChange: %.2f;', log2FoldChange),
                        sprintf('Up: %1.0f; Down: %1.0f;', dim(up_DEG)[1], dim(down_DEG)[1]),
                        sprintf('Total: %1.0f', dim(dif)[1])))
ggsave('Volcano.png',width=6,height=6)
ggsave('V.volcano.pdf',width=6,height=6)
##############-Heatmap-######
gene<-data_repel
select.expr<-countData[rownames(gene),]
library(data.table)
mat<- t(scale(t(select.expr)))#归一化
mat[mat < -1] <- -1
mat[mat > 1] <- 1
annotation_col=data.frame(row.names=rownames(colData), Condition=colData$group)
annotation_row<-data.frame(row.names=rownames(data_repel),sig=data_repel$sig)
#
pdf('Heatmap.pdf',w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Condition= annotation_col$Condition, col = list(Condition = c("Treat" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),show_column_names = F,name = "mat", height = unit(6, "cm"),col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()
png('Heatmap.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Condition= annotation_col$Condition, col = list(Condition = c("Treat" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9),show_column_names = F,name = "mat", height = unit(6, "cm"),col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()
getwd()
#
####################-Enrichment-#######################
#DATA: up_DEG, down_DEG
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtext)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
#########################-UP-##########################################
GO_database <- 'org.Hs.eg.db'
#data
gene <- bitr(up_DEG$symbol,
             fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb = GO_database)
eG <- enrichGO(gene$ENTREZID, #需要分析的基因的EntrezID
               OrgDb = org.Hs.eg.db, #人基因数据库
               pvalueCutoff =0.5, #设置pvalue界值
               qvalueCutoff = 0.5, #设置qvalue界值(FDR校正后的p值）
               ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
               readable =T)
png('up_GO.png',w=10,h=8,units='in',res=600)
dotplot(eG,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =10)
dev.off()
#KEGG
KEGG_database <- 'hsa' 
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism =KEGG_database,
                 pvalueCutoff =0.5,
                 qvalueCutoff = 0.5,
                 keyType = "kegg")
#write.csv(KEGG,file="07.KEGG.csv")
png('up_KEGG.png',w=10,h=8,units='in',res=600)
barplot(kk,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =10)
dev.off()
#########################-Down-##########################################
GO_database <- 'org.Hs.eg.db'
#data
gene <- bitr(down_DEG$symbol,
             fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb = GO_database)
eG <- enrichGO(gene$ENTREZID, #需要分析的基因的EntrezID
               OrgDb = org.Hs.eg.db, #人基因数据库
               pvalueCutoff =0.5, #设置pvalue界值
               qvalueCutoff = 0.5, #设置qvalue界值(FDR校正后的p值）
               ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
               readable =T)
png('down_GO.png',w=10,h=8,units='in',res=600)
dotplot(eG,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =10)
dev.off()
#KEGG
KEGG_database <- 'hsa' 
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism =KEGG_database,
                 pvalueCutoff =0.5,
                 qvalueCutoff = 0.5,
                 keyType = "kegg")
#write.csv(KEGG,file="07.KEGG.csv")
png('down_KEGG.png',w=10,h=8,units='in',res=600)
barplot(kk,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =10)
dev.off()
############################-GSEA-##########################
getwd()
# 加载所需的R包
library(fgsea)          # GSEA分析主程序
library(data.table)     # 数据处理
library(ggplot2)        # 画图处理
library(dplyr)          # 数据处理
library(msigdb)         # 包含基因集合，通常和GSEA分析共同使用
library(GSEABase)       # 可以提供GSEA基础结构和函数,也会被其他包调用

DEG_DESeq2 <- Res
# 读取差异表达的数据，做排序和权重（一般是根据log2FC来作为排序依据）
head(DEG_DESeq2)
# 得到一个排序权重变量，一维有名称变量，后续分析要用

# 以文本格式保存结果
#fwrite(fgseaRes, file = "fgseaRes.txt", sep = "\t", sep2 = c("", " ", ""))
#GSEA Analysis
# 加载所需的R包
library(org.Hs.eg.db) # human的OrgDB
library(clusterProfiler)
# ID转化
gene_entrezid <- bitr(geneID = rownames(DEG_DESeq2), 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", # 转成ENTREZID
                      OrgDb = "org.Hs.eg.db"
)

head(gene_entrezid)
gene_entrezid$logFC <- DEG_DESeq2$log2FoldChange[match(gene_entrezid$SYMBOL,                                                        rownames(DEG_DESeq2))]
genelist = gene_entrezid$logFC
genelist = sort(genelist, decreasing = TRUE)
names(genelist) = gene_entrezid$ENTREZID 
###################################################################
library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH"
)
library(enrichplot)
library(ggplot2)
gsearank(gsea_res,
         geneSetID = 1 # 要展示的基因集
)

#plot
#把entrezid变为symbol
gsea_res_symbol <- setReadable(gsea_res, "org.Hs.eg.db", "ENTREZID")

#multiple pathways
tmp <- as.data.frame(gsea_res_symbol)
png('REACTOME_CELLULAR_SENESCENCE.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "REACTOME_CELLULAR_SENESCENCE",
          geneSetID = "REACTOME_CELLULAR_SENESCENCE")
dev.off()
#
png('REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS",
          geneSetID = "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS")
dev.off()
#
png('REACTOME_APOPTOSIS.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "REACTOME_APOPTOSIS",
          geneSetID = "REACTOME_APOPTOSIS"
)
dev.off()
#
png('REACTOME_CELLULAR_SENESCENCE.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "REACTOME_CELLULAR_SENESCENCE",
          geneSetID = "REACTOME_CELLULAR_SENESCENCE"
)
dev.off()
#
png('REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
          geneSetID = "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP"
)
dev.off()
#
png('REACTOME_INFLAMMASOMES.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "REACTOME_INFLAMMASOMES",
          geneSetID = "REACTOME_INFLAMMASOMES"
)
dev.off()
#
png('GOBP_CELL_CYCLE.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "GOBP_CELL_CYCLE",
          geneSetID = "GOBP_CELL_CYCLE"
)
dev.off()
#
png('TOOKER_GEMCITABINE_RESISTANCE_UP.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "TOOKER_GEMCITABINE_RESISTANCE_UP",
          geneSetID = "TOOKER_GEMCITABINE_RESISTANCE_UP"
)
dev.off()
#
png('REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS",
          geneSetID = "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS"
)
dev.off()
#
png('KEGG_TGF_BETA_SIGNALING_PATHWAY.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "KEGG_TGF_BETA_SIGNALING_PATHWAY",
          geneSetID = "KEGG_TGF_BETA_SIGNALING_PATHWAY"
)
dev.off()
#
png('KEGG_WNT_SIGNALING_PATHWAY.png',w=10,h=8,units='in',res=600)
gseaplot2(gsea_res,
          title = "KEGG_WNT_SIGNALING_PATHWAY",
          geneSetID = "KEGG_WNT_SIGNALING_PATHWAY"
)
dev.off()
#END
