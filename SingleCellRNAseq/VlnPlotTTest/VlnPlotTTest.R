getwd()
setwd("/media/desk16/tly6103/Downloads/MouseLiver20250128/humanout/")

head(sce.all@meta.data)


sub.healthy.KC  = subset(sce.all,cellType == "KC" & group == "healthy")
sub.cirrhotic.KC  = subset(sce.all,cellType == "KC" & group == "cirrhotic")
sub.cirrhotic.Mo  = subset(sce.all,cellType == "Mo" & group == "cirrhotic")

# add information to identify dataset of origin
sub.healthy.KC$dataset <- 'healthy KC'
sub.cirrhotic.KC$dataset <- 'cirrhosis KC'
sub.cirrhotic.Mo$dataset <- 'cirrhosis MoKC'

# merge
Combined <- merge(
  x = sub.healthy.KC,
  y = list(sub.cirrhotic.KC, sub.cirrhotic.Mo),
  add.cell.ids = c("healthy.KC", "cirrhotic.KC", "sub.cirrhotic.Mo")
)
head(Combined@meta.data)

saveRDS(Combined,"Combined.rds")
Combined <- readRDS("Combined.rds")
head(Combined@meta.data)
#根据需要选择进行比较的分组
my_comparisons <- list( c("healthy KC", "cirrhosis KC"), 
                        c("healthy KC", "cirrhosis MoKC"), 
                        c("cirrhosis KC", "cirrhosis MoKC"))

library(ggpubr)

Combined.PHGDH <- subset(x = Combined, subset = (PHGDH > 0))
markers.to.plot <- c("PHGDH")
pdf('humanPHGDH.pdf',w=6,h=6,family='Times')
VlnPlot(Combined.PHGDH, 
        features = markers.to.plot, 
        #        split.by = "dataset", 
        group.by = c("dataset"),
        pt.size = 1, 
        combine = T)+
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  ylim(0, 2)+
  stat_summary(fun = mean, geom = "point", color = "red", size = 3)
dev.off()


Combined.SHMT2 <- subset(x = Combined, subset = (SHMT2 > 0))
markers.to.plot <- c("SHMT2")
pdf('humanSHMT2.pdf',w=6,h=6,family='Times')
VlnPlot(Combined.SHMT2, 
        features = markers.to.plot, 
        #        split.by = "dataset", 
        group.by = c("dataset"),
        pt.size = 1, 
        combine = T)+
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  ylim(0, 3)+
  stat_summary(fun = mean, geom = "point", color = "red", size = 3)
dev.off()

