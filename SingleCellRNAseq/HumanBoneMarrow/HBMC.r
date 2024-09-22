rm(list = ls())
getwd()
setwd("../")
setwd("~/Downloads/HumanBoneMarrow20240913/")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
###-1.Packages-###
library(Seurat)
library(harmony)
library(celldex)
library(SingleR)
ls("package:celldex")
library(BiocParallel)
library(tinyarray)
if(!require("multtest"))BiocManager::install('multtest')
if(!require("metap"))install.packages('metap')
library(ggplot2)
library(ggrepel)
###-2.Data-###
getwd()
f = dir("./data/")
scelist = list()
for(i in 1:length(f)){
  pda <- Read10X(paste0("./data/",f[[i]]))
  scelist[[i]] <- CreateSeuratObject(counts = pda, 
                                     project = f[[i]])
  print(dim(scelist[[i]]))
}
sce.all = merge(scelist[[1]],scelist[-1])
sce.all = JoinLayers(sce.all)
head(sce.all@meta.data)
#cell numbers
table(sce.all$orig.ident)
sum(table(Idents(sce.all)))
###-3.Quality Control-###
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
sce.all[["percent.hb"]] <- PercentageFeatureSet(sce.all, pattern = "^HB[^(P)]")
head(sce.all@meta.data, 3)
VlnPlot(sce.all, 
        features = c("nFeature_RNA",
                     "nCount_RNA"),
        ncol = 3,pt.size = 0, group.by = "orig.ident")
#Not Show
# "percent.mt",
#"percent.rp",
#"percent.hb"
sce.all = subset(sce.all,percent.mt<25)

###-4.Integrated dimensionless clustering-###
f = "obj.Rdata"
if(!file.exists(f)){
  sce.all = sce.all %>% 
    NormalizeData() %>%  
    FindVariableFeatures() %>%  
    ScaleData(features = rownames(.)) %>%  
    RunPCA(pc.genes = VariableFeatures(.))  %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(dims = 1:15, reduction = "harmony") %>% 
    FindClusters(resolution = 0.5) %>% 
    RunUMAP(dims = 1:15,reduction = "harmony") %>% 
    RunTSNE(dims = 1:15,reduction = "harmony")
  save(sce.all,file = f)
}
load(f)
ElbowPlot(sce.all)
#Cluster
TSNEPlot(sce.all,label = T)
UMAPPlot(sce.all,label = T)
###-5.Anotation-###
ls("package:celldex")

#[1] "BlueprintEncodeData"              "DatabaseImmuneCellExpressionData"
#[3] "defineTextQuery"                  "fetchLatestVersion"              
#[5] "fetchMetadata"                    "fetchReference"                  
#[7] "HumanPrimaryCellAtlasData"        "ImmGenData"                      
#[9] "listReferences"                   "listVersions"                    
#[11] "MonacoImmuneData"                 "MouseRNAseqData"                 
#[13] "NovershternHematopoieticData"     "saveReference"                   
#[15] "searchReferences"                 "surveyReferences" 

f = "ref_HumanPrimaryCellAtlasData.RData"
if(!file.exists(f)){
  ref <- celldex::HumanPrimaryCellAtlasData()
  save(ref,file = f)
}
f = "ref_BlueprintEncode.RData"
if(!file.exists(f)){
  ref <- celldex::HumanPrimaryCellAtlasData()
  save(ref,file = f)
}
f = "ref_DatabaseImmuneCellExpressionData.RData"
if(!file.exists(f)){
  ref <- celldex::DatabaseImmuneCellExpressionData()
  save(ref,file = f)
}
#
f = "ref_ImmGenData.RData"
if(!file.exists(f)){
  ref <- celldex::ImmGenData()
  save(ref,file = f)
}
#
ref <- list(get(load("ref_BlueprintEncode.RData")),
            get(load("ref_ImmGenData.RData")),
            get(load("ref_HumanPrimaryCellAtlasData.RData")),
            get(load("ref_DatabaseImmuneCellExpressionData.RData")))
#
library(BiocParallel)
scRNA = sce.all
test = scRNA@assays$RNA$data
pred.scRNA <- SingleR(test = test, 
                      ref = ref,
                      labels = list(ref[[1]]$label.main,
                                    ref[[2]]$label.main,
                                    ref[[3]]$label.main,
                                    ref[[4]]$label.main), 
                      clusters = scRNA@active.ident)
pred.scRNA$pruned.labels
#查看注释准确性 
#plotScoreHeatmap(pred.scRNA, clusters=pred.scRNA@rownames, fontsize.row = 9,show_colnames = T)
new.cluster.ids <- pred.scRNA$pruned.labels
#new.cluster.ids[is.na(new.cluster.ids)] = "unknown"
names(new.cluster.ids) <- levels(scRNA)
levels(scRNA)
scRNA <- RenameIdents(scRNA,new.cluster.ids)
levels(scRNA)
p2 <- DimPlot(scRNA, reduction = "umap",label = T,pt.size = 0.5)# + NoLegend()
p2

FeaturePlot(scRNA, features = c("Jchain","Cd68", "Themis", "Cd34", "Cd8a", 
                                "Bmp7"))
#Monocyte
VlnPlot(JoinLayers(sce.all), features = c("Ly6a"))
#Memmencychmel progenitor cell 
VlnPlot(JoinLayers(sce.all), features = c( "Cd44"))
#Oteoclast
VlnPlot(JoinLayers(sce.all), features = c( "Acp5",  "Nfatc1"))
#Osteocyte 14
VlnPlot(JoinLayers(sce.all), features = c( "Sost"))
#Dendritic 14
VlnPlot(JoinLayers(sce.all), features = c("Slamf7"))
#Msc 10
VlnPlot(JoinLayers(sce.all), features = c("Cd34"))
#Osteoclast precursor cell
VlnPlot(JoinLayers(sce.all), features = c("Galectin3"))


#
new.cluster.ids <- c("Neutrophils", "Neutrophils", "Monocyte", "Monocyte", "Macrophages",
                     "Eosinohpils", "Eosinohpils", "Eosinohpils", "Monocyte", "Oteoclast",
                     "Mesenchymal stem cell", "Monocyte","Monocyte","Monocyte","Dendritic", 
                     "Monocyte", "Monocyte", "Basophils")
names(new.cluster.ids) <- levels(sce.all)
sce.all <- RenameIdents(sce.all, new.cluster.ids)
DimPlot(JoinLayers(sce.all), reduction = "umap", label = TRUE, pt.size = 0.5)# + NoLegend()
#split
library(tinyarray)
sce.all$seurat_annotations = Idents(sce.all)
table(sce.all$orig.ident)
sce.all$group = ifelse(sce.all$orig.ident %in% c("KO", "WT"), 
                       "treat", "control")
DimPlot(JoinLayers(sce.all), reduction = "umap", group.by = "group")
#
scRNA$seurat_annotations = Idents(scRNA)
table(scRNA$orig.ident)
# 可视化
DimPlot(JoinLayers(sce.all), reduction = "umap", 
        group.by = c("group", "seurat_annotations"))

## 

library(tinyarray)
sce.all$group = ifelse(sce.all$orig.ident %in% c("KO"), "KO","WT")
DimPlot(sce.all, reduction = "umap", group.by = "group")
DimPlot(JoinLayers(sce.all), reduction = "umap", 
        group.by = c("group", "seurat_annotations"))
#
# 每种细胞的数量和比例
cell_counts <- table(Idents(sce.all))
cell.all <- cbind(cell_counts = cell_counts, 
                  cell_Freq = round(prop.table(cell_counts)*100,2))
#各组中每种细胞的数量和比例
cell.num.group <- table(Idents(sce.all), sce.all$group) 
cell.freq.group <- round(prop.table(cell.num.group, margin = 2) *100,2)
cell.all = cbind(cell.all,cell.num.group,cell.freq.group)
cell.all = cell.all[,c(1,3,4,2,5,6)]
colnames(cell.all) = paste(rep(c("all","control","treat"),times = 2),
                           rep(c("count","freq"),each = 3),sep = "_")
cell.all
#DE
if(!require("multtest"))BiocManager::install('multtest')
if(!require("metap"))install.packages('metap')
table(sce.all$seurat_annotations)
sub.markers <- FindConservedMarkers(sce.all, ident.1 = "Monocyte", 
                                    grouping.var = "group", 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25,
                                    verbose = F)
head(sub.markers)
#bubble
markers.to.plot = c("Ly6a", "Slamf7", "Cd34", "Batf3", "Ifi205", "Mycl", "Slpi", "Cx3cr1",
                    "Cd163", "Irf8", "C1qb", "Cd44", "Ctss", "Cybb", "Sost","Ctsk","Mmp9","Acp5","Runx2") #一组感兴趣的基因
#scRNA <- subset(scRNA, seurat_annotations %in% na.omit(scRNA$seurat_annotations))
DotPlot(sce.all, features = markers.to.plot, cols = c("blue", "red"), 
        dot.scale = 8, split.by = "group") +
  RotatedAxis()
FeaturePlot(sce.all, 
            features = c("Cd34", "Cd44"), 
            split.by = "group", 
            max.cutoff = 2, 
            cols = c("grey",  "red"), 
            reduction = "umap")
#
library(patchwork)
plots <- VlnPlot(sce.all, features = c("Cd34", "Cd44"), 
                 split.by = "group", group.by = "seurat_annotations",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
#DEanalysis 
bulk <- AggregateExpression(scRNA, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("seurat_annotations","orig.ident", "group"))

sub <- subset(bulk, seurat_annotations == "CD8+ T-cells")
Idents(sub) <- "group"
de_markers <- FindMarkers(sub, ident.1 = "treat", ident.2 = "control", slot = "counts", test.use = "DESeq2",
                          verbose = F)
de_markers$gene <- rownames(de_markers)
k1 = de_markers$avg_log2FC< -1 & de_markers$p_val <0.01
k2 = de_markers$avg_log2FC> 1 & de_markers$p_val <0.01
de_markers$change <- ifelse(k1,"down",ifelse(k2,"up","not"))
library(ggplot2)
library(ggrepel)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val),color = change)) + 
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = c(1,-1),linetype = 4)+
  geom_hline(yintercept = -log10(0.01),linetype = 4)+
  scale_color_manual(values = c("blue","grey","red"))+
  theme_bw() +
  ylab("-log10(unadjusted p-value)") 

