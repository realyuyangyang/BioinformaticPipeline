rm(list = ls())
#options(stringsAsFactors = F)
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if(!dir.exists("06.Enrich")){dir.create("06.Enrich")}
setwd("06.Enrich")
getwd()
#-----------------library------------------------------
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
#----------------------data------------------------------
Candidated_genes <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/03.Venn/Candidated_genes.csv")
# enrichment --------------------------------------------------------------
GO_database <- 'org.Hs.eg.db'
search_kegg_organism(str, by = "scientific_name", ignore.case = FALSE)
KEGG_database <- 'hsa' 
###GO富集分析
gene <- bitr(Candidated_genes$symbol,
             fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb = GO_database)
go_result<-enrichGO(gene$ENTREZID,
                    OrgDb = GO_database,
                    keyType = "ENTREZID",             
                    ont = "ALL",  
                    pAdjustMethod = "BH",
                    readable = T)
#
BP <- go_result[which(go_result$ONTOLOGY=='BP'),]
CC <- go_result[which(go_result$ONTOLOGY=='CC'),]
MF <- go_result[which(go_result$ONTOLOGY=='MF'),]
display_number <- 10
go_result_BP = as.data.frame(BP)[1:display_number, ]
go_result_CC = as.data.frame(CC)[1:display_number, ]
go_result_MF = as.data.frame(MF)[1:display_number, ]
go_result_plot <- rbind(go_result_BP, go_result_CC, go_result_MF)
colnames(go_result_plot)
go_result_plot$'-log10(p.adjust)' <- -log10(go_result_plot$p.adjust)
go_result_plot<-na.omit(go_result_plot)
#
#
#GO analyse
eG <- enrichGO(gene$ENTREZID, #需要分析的基因的EntrezID
               OrgDb = org.Hs.eg.db, #人基因数据库
               pvalueCutoff =0.05, #设置pvalue界值
               qvalueCutoff = 0.05, #设置qvalue界值(FDR校正后的p值）
               ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
               readable =T)
#eG <- enrichGO(gene = ENTREZID_lst,
#                OrgDb=org.Hs.eg.db,
#                pvalueCutoff = 0.05,
#                qvalueCutoff = 0.05,
#                readable = TRUE)
write.table(eG,file="eG.txt", sep="\t", quote=F, row.names = F)
# ont = "MF",默认
# ont: One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
dotplot(eG,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =5,#只显示前5
        split="ONTOLOGY") + #以ONTOLOGY类型分开
  facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分屏绘图
library(tidyverse)
eGo = read.table("eG.txt",header=TRUE,sep="\t",quote = "")  #读取第1部分enrichGO分析输出的文件eG。
eGo <- separate(data=eGo, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor 
eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)
library(ggplot2)
p <- ggplot(eGo10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor)))+ 
  geom_point(aes(size=Count,color=qvalue,shape=ONTOLOGY)) +
  scale_color_gradient(low="red",high = "green") + 
  labs(color="pvalue",size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
p + facet_wrap( ~ ONTOLOGY)
ggsave("facet_wrap_go.png", w= 8, h=6, dpi = 600)
dev.off()



library(GOplot)
GOterms = data.frame(category = eGo10$ONTOLOGY,
                     ID = eGo10$ID,
                     term = eGo10$Description, 
                     genes = gsub("/", ",", eGo10$geneID), 
                     adj_pval = eGo10$p.adjust)
signif_dat <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/01.DEGs/DEG_sig.csv")
genelist <- data.frame(ID = signif_dat$symbol, logFC = signif_dat$log2FoldChange) #从已有“数据”中提取genelist，1列ID，1列logFC。
circ <- circle_dat(GOterms, genelist)

GOBubble(circ, labels = 5, # 标注的界值：-log(adjusted p-value) (默认5)
         table.legend =T, #是否显示右侧legend，默认是
         ID=T, # T标注term名称，F标注ID
         display="single") #是否分屏
GOBubble(circ, labels = 5, # 标注的界值：-log(adjusted p-value) (默认5)
         table.legend =F, #不显示右侧legend
         ID=F, # 标注term名称
         display='single') # 不分屏
chord <- chord_dat(circ, #前面的circ数据
                   genelist[1:100,], #选择需要的基因
                   GOterms$term[1:10]) #选择要显示的term
GOChord(chord, space = 0, #弦的间隔
        gene.order = 'logFC', #基因排序方式
        gene.space = 0.25, #基因名称和图的间隔
        gene.size = 5, #基因名的大小
        nlfc =1, #是否有logFC
        border.size= NULL, #彩虹边框大小
        lfc.col=c('red','black','cyan')) #自定义logFC 颜色

GOCircle(circ)
GOCircle(circ,rad1=2, #内环半径
         rad2=3, #外环半径
         label.size= 5, #标签大小
         label.fontface = 'bold', #标签字体
         nsub=10, #显示的terms数，前10个。（画图前需要先把数据整理好，想要哪些term）
         zsc.col = c('red', 'yellow', "blue"), # z-score颜色
         lfc.col = c('red', 'blue')) # 基因up或down的颜色
#class(ego)
#head(ego)





ggplot(go_result_plot,aes(x=Count,y=`-log10(p.adjust)`))+
  geom_point(aes(size=Count,color=Description),show.legend = FALSE,alpha=0.6)+
  scale_size(range=c(5,10))+
  # scale_fill_gradient(low = "DeepSkyBlue",high = "Tomato")+
  ylim(2,15)+
  xlim(0,10)+
  theme_bw()+
  theme_bw() + xlab("Count")+ylab("-log10(p.adjust)")+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
        axis.text.x =element_text(size=16,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=16,family = "Times", face = "bold"),
        strip.text = element_text(size = 16,family = "Times", face = "bold"),
        plot.title=element_text(size=22,family = "Times", face = "bold"),
        legend.title=element_text(size=16,family = "Times", face = "bold"),
        legend.text=element_text(size=14,family = "Times", face = "bold"),
        #legend.position = "bottom",
        
  )+
  #theme(legend.position = c("none"))+
  geom_text_repel(fontface = "italic",
                  size = 4,
                  color = "black",
                  family = "Times", size=3,face = "bold",
                  data =go_result_plot,
                  aes(label = Description),
                  show.legend = FALSE )+
  facet_grid(.~ONTOLOGY)

ggsave('06.cor_point.pdf',w=14,h=8)
ggsave('06.cor_point.png',w=14,h=8)







# GO ----------------------------------------------------------------------
gg<-enrichGO(gene$ENTREZID,
             OrgDb = GO_database,
             keyType = "ENTREZID",             
             ont = "ALL",  
             pAdjustMethod = "BH",
             readable = T)
GO<-data.frame(gg)  ##富集结果
table(GO$ONTOLOGY)  ##  BP  CC  MF 130  47   3
write.csv(GO,file="06.GO.csv")

GO_BP_top10<-GO[GO$ONTOLOGY=='BP',][c(1:20),]
#
library(treemap)
pdf('06.GO_BP.pdf',w=10,h=8)
treemap(GO_BP_top10,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='GO:BP',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()

png('06.GO_BP.png',w=10,h=8,units='in',res=600)
treemap(GO_BP_top10,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='GO:BP',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()

GO_CC_top10<-GO[GO$ONTOLOGY=='CC',][c(1:20),]

library(treemap)
pdf('06.GO_CC.pdf',w=10,h=8)
treemap(GO_CC_top10,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='GO:CC',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()

png('06.GO_CC.png',w=10,h=8,units='in',res=600)
treemap(GO_CC_top10,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='GO:CC',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()

GO_MF_top10<-GO[GO$ONTOLOGY=='MF',][c(1:3),]
library(treemap)
pdf('06.GO_MF.pdf',w=10,h=8)
treemap(GO_MF_top10,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='GO:MF',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()

png('06.GO_MF.png',w=10,h=8,units='in',res=600)
treemap(GO_MF_top10,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='GO:MF',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()
#------------------------------------kegg--------------------------------------
#################################################################################
kk <- enrichKEGG(gene = gene$ENTREZID,organism =KEGG_database, keyType = "kegg")
enrichKK<-DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
KEGG<-data.frame(enrichKK)  ##13个通路
write.csv(KEGG,file="07.KEGG.csv")

library(treemap)
pdf('08.KEGG.pdf',w=10,h=8)
treemap(KEGG,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='KEGG pathways',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()

png('08.KEGG.png',w=10,h=8,units='in',res=600)
treemap(KEGG,
        fontsize.title = 14,
        fontsize.labels = 13,
        fontsize.legend = 17,
        title='KEGG pathways',
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        #   palette = "Spectral", #自定义颜色画板
        type="value", #指定颜色填充数据的类型
        format.legend = list(scientific = T, big.mark = " ")
)
dev.off()

# Load ggplot2
library(ggplot2)

# Create data
data <- data.frame(
  name=c("A","B","C","D","E") ,  
  value=c(3,12,5,18,45)
)

# Barplot
p<- ggplot(data, aes(x=name, y=value)) + 
  geom_bar(stat = "identity")
p
ggsave("test.png", w= 8, h=6)
dev.off()