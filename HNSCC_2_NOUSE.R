library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)


for (i in 1:26) {
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),Read10X(paste("G:/GSE139324/HNSCC_",i,"_PBMC",sep = "")))
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),as.data.frame(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))
  g<-get(paste("HNSCC_",i,"_PBMC.data",sep = ""))
  for (j in 1:length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))) {
    colnames(g)[j]<-paste(colnames(g[j]),"PBMC",i,j,sep = "-")
  }
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),g)
  assign(paste("HNSCC_",i,"_PBMC.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))),rep("PBMC",length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = "")))))))
  a<-get(paste("HNSCC_",i,"_PBMC.metadata",sep = ""))
  b<-as.data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))
  colnames(a)<-c("barcode","group")
  rownames(a)<-b[,1]
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),CreateSeuratObject(counts = g, project = paste("HNSCC_",i,"_PBMC",sep = ""),meta.data = a))
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),NormalizeData(get(paste("HNSCC_",i,"_PBMC",sep = ""))))
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),FindVariableFeatures(get(paste("HNSCC_",i,"_PBMC",sep = "")), selection.method = "vst", nfeatures = 2000))
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),Read10X(paste("G:/GSE139324/HNSCC_",i,"_TIL",sep = "")))
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),as.data.frame(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))
  h<-get(paste("HNSCC_",i,"_TIL.data",sep = ""))
  for (j in 1:length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))) {
    colnames(h)[j]<-paste(colnames(h[j]),"TIL",i,j,sep = "-")
  }
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),h)
  assign(paste("HNSCC_",i,"_TIL.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))),rep("TIL",length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = "")))))))
  c<-get(paste("HNSCC_",i,"_TIL.metadata",sep = ""))
  d<-as.data.frame(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))
  colnames(c)<-c("barcode","group")
  rownames(c)<-d[,1]
  assign(paste("HNSCC_",i,"_TIL",sep = ""),CreateSeuratObject(counts = h, project = paste("HNSCC_",i,"_TIL",sep = ""),meta.data = c))
  assign(paste("HNSCC_",i,"_TIL",sep = ""),NormalizeData(get(paste("HNSCC_",i,"_TIL",sep = ""))))
  assign(paste("HNSCC_",i,"_TIL",sep = ""),FindVariableFeatures(get(paste("HNSCC_",i,"_TIL",sep = "")), selection.method = "vst", nfeatures = 2000))
}

e<-paste("HNSCC_",1,"_PBMC.metadata",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_PBMC.metadata",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_PBMC.metadata<-rbind(HNSCC_1_PBMC.metadata, HNSCC_2_PBMC.metadata, HNSCC_3_PBMC.metadata, HNSCC_4_PBMC.metadata, HNSCC_5_PBMC.metadata, HNSCC_6_PBMC.metadata, HNSCC_7_PBMC.metadata, HNSCC_8_PBMC.metadata, HNSCC_9_PBMC.metadata, HNSCC_10_PBMC.metadata, HNSCC_11_PBMC.metadata, HNSCC_12_PBMC.metadata, HNSCC_13_PBMC.metadata, HNSCC_14_PBMC.metadata, HNSCC_15_PBMC.metadata, HNSCC_16_PBMC.metadata, HNSCC_17_PBMC.metadata, HNSCC_18_PBMC.metadata, HNSCC_19_PBMC.metadata, HNSCC_20_PBMC.metadata, HNSCC_21_PBMC.metadata, HNSCC_22_PBMC.metadata, HNSCC_23_PBMC.metadata, HNSCC_24_PBMC.metadata, HNSCC_25_PBMC.metadata, HNSCC_26_PBMC.metadata)

e<-paste("HNSCC_",1,"_PBMC.data",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_PBMC.data",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_PBMC.data<-cbind(HNSCC_1_PBMC.data, HNSCC_2_PBMC.data, HNSCC_3_PBMC.data, HNSCC_4_PBMC.data, HNSCC_5_PBMC.data, HNSCC_6_PBMC.data, HNSCC_7_PBMC.data, HNSCC_8_PBMC.data, HNSCC_9_PBMC.data, HNSCC_10_PBMC.data, HNSCC_11_PBMC.data, HNSCC_12_PBMC.data, HNSCC_13_PBMC.data, HNSCC_14_PBMC.data, HNSCC_15_PBMC.data, HNSCC_16_PBMC.data, HNSCC_17_PBMC.data, HNSCC_18_PBMC.data, HNSCC_19_PBMC.data, HNSCC_20_PBMC.data, HNSCC_21_PBMC.data, HNSCC_22_PBMC.data, HNSCC_23_PBMC.data, HNSCC_24_PBMC.data, HNSCC_25_PBMC.data, HNSCC_26_PBMC.data)

e<-paste("HNSCC_",1,"_TIL.metadata",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_TIL.metadata",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_TIL.metadata<-rbind(HNSCC_1_TIL.metadata, HNSCC_2_TIL.metadata, HNSCC_3_TIL.metadata, HNSCC_4_TIL.metadata, HNSCC_5_TIL.metadata, HNSCC_6_TIL.metadata, HNSCC_7_TIL.metadata, HNSCC_8_TIL.metadata, HNSCC_9_TIL.metadata, HNSCC_10_TIL.metadata, HNSCC_11_TIL.metadata, HNSCC_12_TIL.metadata, HNSCC_13_TIL.metadata, HNSCC_14_TIL.metadata, HNSCC_15_TIL.metadata, HNSCC_16_TIL.metadata, HNSCC_17_TIL.metadata, HNSCC_18_TIL.metadata, HNSCC_19_TIL.metadata, HNSCC_20_TIL.metadata, HNSCC_21_TIL.metadata, HNSCC_22_TIL.metadata, HNSCC_23_TIL.metadata, HNSCC_24_TIL.metadata, HNSCC_25_TIL.metadata, HNSCC_26_TIL.metadata)

e<-paste("HNSCC_",1,"_TIL.data",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_TIL.data",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_TIL.data<-cbind(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)

HNSCC_PBMC<-CreateSeuratObject(counts = HNSCC_PBMC.data,project = "HNSCC_PBMC",meta.data = HNSCC_PBMC.metadata)
HNSCC_PBMC<-NormalizeData(HNSCC_PBMC)
remove(HNSCC_1_PBMC.data, HNSCC_2_PBMC.data, HNSCC_3_PBMC.data, HNSCC_4_PBMC.data, HNSCC_5_PBMC.data, HNSCC_6_PBMC.data, HNSCC_7_PBMC.data, HNSCC_8_PBMC.data, HNSCC_9_PBMC.data, HNSCC_10_PBMC.data, HNSCC_11_PBMC.data, HNSCC_12_PBMC.data, HNSCC_13_PBMC.data, HNSCC_14_PBMC.data, HNSCC_15_PBMC.data, HNSCC_16_PBMC.data, HNSCC_17_PBMC.data, HNSCC_18_PBMC.data, HNSCC_19_PBMC.data, HNSCC_20_PBMC.data, HNSCC_21_PBMC.data, HNSCC_22_PBMC.data, HNSCC_23_PBMC.data, HNSCC_24_PBMC.data, HNSCC_25_PBMC.data, HNSCC_26_PBMC.data)
remove(HNSCC_1_PBMC.metadata, HNSCC_2_PBMC.metadata, HNSCC_3_PBMC.metadata, HNSCC_4_PBMC.metadata, HNSCC_5_PBMC.metadata, HNSCC_6_PBMC.metadata, HNSCC_7_PBMC.metadata, HNSCC_8_PBMC.metadata, HNSCC_9_PBMC.metadata, HNSCC_10_PBMC.metadata, HNSCC_11_PBMC.metadata, HNSCC_12_PBMC.metadata, HNSCC_13_PBMC.metadata, HNSCC_14_PBMC.metadata, HNSCC_15_PBMC.metadata, HNSCC_16_PBMC.metadata, HNSCC_17_PBMC.metadata, HNSCC_18_PBMC.metadata, HNSCC_19_PBMC.metadata, HNSCC_20_PBMC.metadata, HNSCC_21_PBMC.metadata, HNSCC_22_PBMC.metadata, HNSCC_23_PBMC.metadata, HNSCC_24_PBMC.metadata, HNSCC_25_PBMC.metadata, HNSCC_26_PBMC.metadata)
remove(HNSCC_1_PBMC, HNSCC_2_PBMC, HNSCC_3_PBMC, HNSCC_4_PBMC, HNSCC_5_PBMC, HNSCC_6_PBMC, HNSCC_7_PBMC, HNSCC_8_PBMC, HNSCC_9_PBMC, HNSCC_10_PBMC, HNSCC_11_PBMC, HNSCC_12_PBMC, HNSCC_13_PBMC, HNSCC_14_PBMC, HNSCC_15_PBMC, HNSCC_16_PBMC, HNSCC_17_PBMC, HNSCC_18_PBMC, HNSCC_19_PBMC, HNSCC_20_PBMC, HNSCC_21_PBMC, HNSCC_22_PBMC, HNSCC_23_PBMC, HNSCC_24_PBMC, HNSCC_25_PBMC, HNSCC_26_PBMC)
remove(HNSCC_1_TIL, HNSCC_2_TIL, HNSCC_3_TIL, HNSCC_4_TIL, HNSCC_5_TIL, HNSCC_6_TIL, HNSCC_7_TIL, HNSCC_8_TIL, HNSCC_9_TIL, HNSCC_10_TIL, HNSCC_11_TIL, HNSCC_12_TIL, HNSCC_13_TIL, HNSCC_14_TIL, HNSCC_15_TIL, HNSCC_16_TIL, HNSCC_17_TIL, HNSCC_18_TIL, HNSCC_19_TIL, HNSCC_20_TIL, HNSCC_21_TIL, HNSCC_22_TIL, HNSCC_23_TIL, HNSCC_24_TIL, HNSCC_25_TIL, HNSCC_26_TIL)
remove(a,b,c,d,e,f,g,h,i,j)
HNSCC_TIL<-CreateSeuratObject(counts = HNSCC_TIL.data,project = "HNSCC_TIL",meta.data = HNSCC_TIL.metadata)
HNSCC_TIL<-NormalizeData(HNSCC_TIL)
remove(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)
remove(HNSCC_1_TIL.metadata, HNSCC_2_TIL.metadata, HNSCC_3_TIL.metadata, HNSCC_4_TIL.metadata, HNSCC_5_TIL.metadata, HNSCC_6_TIL.metadata, HNSCC_7_TIL.metadata, HNSCC_8_TIL.metadata, HNSCC_9_TIL.metadata, HNSCC_10_TIL.metadata, HNSCC_11_TIL.metadata, HNSCC_12_TIL.metadata, HNSCC_13_TIL.metadata, HNSCC_14_TIL.metadata, HNSCC_15_TIL.metadata, HNSCC_16_TIL.metadata, HNSCC_17_TIL.metadata, HNSCC_18_TIL.metadata, HNSCC_19_TIL.metadata, HNSCC_20_TIL.metadata, HNSCC_21_TIL.metadata, HNSCC_22_TIL.metadata, HNSCC_23_TIL.metadata, HNSCC_24_TIL.metadata, HNSCC_25_TIL.metadata, HNSCC_26_TIL.metadata)
remove(HNSCC_PBMC.data,HNSCC_PBMC.metadata,HNSCC_TIL.data,HNSCC_TIL.metadata)

immune.anchor <- FindIntegrationAnchors(object.list = list(HNSCC_PBMC,HNSCC_TIL), dims = 1:100,anchor.features = 33694)
immune.combined <- IntegrateData(anchorset = immune.anchor)
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
DefaultAssay(immune.combined) <- "integrated"
VlnPlot(immune.combined,features = c("JMJD1C","NRP1"),pt.size = 0)

