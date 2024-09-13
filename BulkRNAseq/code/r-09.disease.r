#候选基因筛选
rm(list=ls()) 
setwd("/data/nas1/yuyangyang/Projects/YQKM-10411-2/")
if(!dir.exists("09.disease")){dir.create("09.disease")}
setwd("09.disease")
getwd()

getwd()
library(ggplot2)
library(forcats)
library("tidyverse")
library("ggrepel")
Data <- read.csv("/data/nas1/yuyangyang/Projects/YQKM-10411-2/geneassociatedisease.csv", 
                 row.names=NULL)
#ALOX5
data <- Data[1:5,]
data %>%
  mutate(Associate.Disease = fct_reorder(Associate.Disease, InferenceScore)) %>%
  ggplot(aes(x= Associate.Disease, y= InferenceScore)) +
  geom_bar(stat="identity", fill="firebrick3")+
  geom_text(aes(label=InferenceScore), vjust=1.6, color="white", size=3.5)+
  labs(title="The Interaction of ALOX5 and  disease", x=" ", y = "Inference Score")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75)) 
ggsave('ALOX5.png',width=8,height=6, dpi = 600)
#ggsave('DNM1L.jpg',width=8,height=6, dpi = 600)
#ggsave('DNM1L.pdf',width=8,height=6, dpi = 600)
dev.off()
#HMOX1
data <- Data[6:10,]
data %>%
  mutate(Associate.Disease = fct_reorder(Associate.Disease, InferenceScore)) %>%
  ggplot(aes(x= Associate.Disease, y= InferenceScore)) +
  geom_bar(stat="identity", fill="firebrick3")+
  geom_text(aes(label=InferenceScore), vjust=1.6, color="white", size=3.5)+
  labs(title="The Interaction of HMOX1 and  disease", x=" ", y = "Inference Score")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75)) 
ggsave('HMOX1.png',width=8,height=6, dpi = 600)
#ggsave('DNM1L.jpg',width=8,height=6, dpi = 600)
#ggsave('DNM1L.pdf',width=8,height=6, dpi = 600)
dev.off()
data <- Data[11:15,]
data %>%
  mutate(Associate.Disease = fct_reorder(Associate.Disease, InferenceScore)) %>%
  ggplot(aes(x= Associate.Disease, y= InferenceScore)) +
  geom_bar(stat="identity", fill="firebrick3")+
  geom_text(aes(label=InferenceScore), vjust=1.6, color="white", size=3.5)+
  labs(title="The Interaction of PLA2G7 and  disease", x=" ", y = "Inference Score")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75)) 
ggsave('PLA2G7.png',width=8,height=6, dpi = 600)
#ggsave('DNM1L.jpg',width=8,height=6, dpi = 600)
#ggsave('DNM1L.pdf',width=8,height=6, dpi = 600)
dev.off()
#MFF
data <- Data[6:10,]##############data
data %>%
  mutate(MaleGenitalDisease = fct_reorder(MaleGenitalDisease, InferenceScore)) %>%
  ggplot(aes(x= MaleGenitalDisease, y= InferenceScore)) +
  geom_bar(stat="identity", fill="firebrick3")+
  geom_text(aes(label=InferenceScore), vjust=1.6, color="white", size=3.5)+
  labs(title="The Interaction of MFF and male genital disease", x=" ", y = "Inference Score")+#title
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75)) 
ggsave('MFF.png',width=8,height=6)
ggsave('MFF.jpg',width=8,height=6, dpi = 600)
ggsave('MFF.pdf',width=8,height=6, dpi = 600)##############name
dev.off()

#NDUFS3
data <- Data[11:15,]##############data
data %>%
  mutate(MaleGenitalDisease = fct_reorder(MaleGenitalDisease, InferenceScore)) %>%
  ggplot(aes(x= MaleGenitalDisease, y= InferenceScore)) +
  geom_bar(stat="identity", fill="firebrick3")+
  geom_text(aes(label=InferenceScore), vjust=1.6, color="white", size=3.5)+
  labs(title="The Interaction of NDUFS3 and male genital disease", x=" ", y = "Inference Score")+#title
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75)) 
ggsave('NDUFS3.png',width=8,height=6)
ggsave('NDUFS3.jpg',width=8,height=6, dpi = 600)
ggsave('NDUFS3.pdf',width=8,height=6, dpi = 600)##############name
dev.off()

#PRKN
data <- Data[16:20,]##############data
data %>%
  mutate(MaleGenitalDisease = fct_reorder(MaleGenitalDisease, InferenceScore)) %>%
  ggplot(aes(x= MaleGenitalDisease, y= InferenceScore)) +
  geom_bar(stat="identity", fill="firebrick3")+
  geom_text(aes(label=InferenceScore), vjust=1.6, color="white", size=3.5)+
  labs(title="The Interaction of PRKN and male genital disease", x=" ", y = "Inference Score")+#title
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75)) 
ggsave('PRKN.png',width=8,height=6)
ggsave('PRKN.jpg',width=8,height=6, dpi = 600)
ggsave('PRKN.pdf',width=8,height=6, dpi = 600)##############name
dev.off()

#RHOT1
data <- Data[21:25,]##############data
data %>%
  mutate(MaleGenitalDisease = fct_reorder(MaleGenitalDisease, InferenceScore)) %>%
  ggplot(aes(x= MaleGenitalDisease, y= InferenceScore)) +
  geom_bar(stat="identity", fill="firebrick3")+
  geom_text(aes(label=InferenceScore), vjust=1.6, color="white", size=3.5)+
  labs(title="The Interaction of RHOT1 and male genital disease", x=" ", y = "Inference Score")+#title
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75)) 
ggsave('RHOT1.png',width=8,height=6)
ggsave('RHOT1.jpg',width=8,height=6, dpi = 600)
ggsave('RHOT1.pdf',width=8,height=6, dpi = 600)##############name
dev.off()

