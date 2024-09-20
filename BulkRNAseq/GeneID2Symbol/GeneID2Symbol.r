library(dplyr)
#drop_duplicates
data <- read.csv("~/Downloads/HCC/data/GSE143233_resistantHCC_genes_fpkm.csv",
                                                      stringsAsFactors=TRUE)
countData <- data %>% group_by(symbol) %>% filter (! duplicated(symbol))
setwd("../data/")
write.csv(countData, "countData.csv", row.names = F)
annot <- read.delim("~/Downloads/HCC/data/Human.GRCh38.p13.annot.tsv", 
                    row.names=NULL, stringsAsFactors=TRUE)
raw_counts <- read.delim("~/Downloads/HCC/data/GSE143233_raw_counts_GRCh38.p13_NCBI.tsv", 
                                         row.names=NULL, stringsAsFactors=TRUE)
#
#IDREF -> Symbol
fdata_2 <- na.omit(annot)
ids <- fdata_2[, c("GeneID", "Symbol")] 
#
# install.packages("dplyr")
colnames(ids) = c("GeneID" ,"Symbol")
rownames(raw_counts) = raw_counts$GeneID
exp <- na.omit(raw_counts)
exp2= merge(exp,ids,by.x="GeneID", by.y="GeneID") # 合并数据********
exp2=exp2[!duplicated(exp2$Symbol),]# 按照symbol列去重
rownames(exp2)=exp2$Symbol # 数据框probe_exp的行名变成symbol
gene_exp_matrix <- na.omit(exp2 ) # 去空值
countData_SYMBOL <- gene_exp_matrix[,c(8, 2:7)]
#输出文件
setwd("../data/")
write.csv(countData_SYMBOL,file = "countData_SYMBOL.csv",row.names=F,col.names = T)