library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)


for (i in 1:26) {
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),Read10X(paste("D:/Jmjd1c_Treg_Tumor/GSE139324/HNSCC_",i,"_PBMC",sep = "")))
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),as.data.frame(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))
  g<-get(paste("HNSCC_",i,"_PBMC.data",sep = ""))
  for (j in 1:length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))) {
    colnames(g)[j]<-paste(colnames(g[j]),"PBMC",i,j,sep = "-")
  }
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),g)
  assign(paste("HNSCC_",i,"_PBMC.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))),rep("PBMC",length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))),rep(i,length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = "")))))))
  a<-get(paste("HNSCC_",i,"_PBMC.metadata",sep = ""))
  b<-as.data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))
  colnames(a)<-c("barcode","group","sample")
  rownames(a)<-b[,1]
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),Read10X(paste("D:/Jmjd1c_Treg_Tumor/GSE139324/HNSCC_",i,"_TIL",sep = "")))
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),as.data.frame(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))
  h<-get(paste("HNSCC_",i,"_TIL.data",sep = ""))
  for (k in 1:length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))) {
    colnames(h)[k]<-paste(colnames(h[k]),"TIL",i,k,sep = "-")
  }
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),h)
  assign(paste("HNSCC_",i,"_TIL.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))),rep("TIL",length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))),rep(i,length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = "")))))))
  c<-get(paste("HNSCC_",i,"_TIL.metadata",sep = ""))
  d<-as.data.frame(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))
  colnames(c)<-c("barcode","group","sample")
  rownames(c)<-d[,1]
}


HNSCC_PBMC.metadata<-rbind(HNSCC_1_PBMC.metadata, HNSCC_2_PBMC.metadata, HNSCC_3_PBMC.metadata, HNSCC_4_PBMC.metadata, HNSCC_5_PBMC.metadata, HNSCC_6_PBMC.metadata, HNSCC_7_PBMC.metadata, HNSCC_8_PBMC.metadata, HNSCC_9_PBMC.metadata, HNSCC_10_PBMC.metadata, HNSCC_11_PBMC.metadata, HNSCC_12_PBMC.metadata, HNSCC_13_PBMC.metadata, HNSCC_14_PBMC.metadata, HNSCC_15_PBMC.metadata, HNSCC_16_PBMC.metadata, HNSCC_17_PBMC.metadata, HNSCC_18_PBMC.metadata, HNSCC_19_PBMC.metadata, HNSCC_20_PBMC.metadata, HNSCC_21_PBMC.metadata, HNSCC_22_PBMC.metadata, HNSCC_23_PBMC.metadata, HNSCC_24_PBMC.metadata, HNSCC_25_PBMC.metadata, HNSCC_26_PBMC.metadata)
colnames(HNSCC_PBMC.metadata)<-c("barcodes","group","sample")
rownames(HNSCC_PBMC.metadata)<-HNSCC_PBMC.metadata[,1]

HNSCC_PBMC.data<-cbind(HNSCC_1_PBMC.data, HNSCC_2_PBMC.data, HNSCC_3_PBMC.data, HNSCC_4_PBMC.data, HNSCC_5_PBMC.data, HNSCC_6_PBMC.data, HNSCC_7_PBMC.data, HNSCC_8_PBMC.data, HNSCC_9_PBMC.data, HNSCC_10_PBMC.data, HNSCC_11_PBMC.data, HNSCC_12_PBMC.data, HNSCC_13_PBMC.data, HNSCC_14_PBMC.data, HNSCC_15_PBMC.data, HNSCC_16_PBMC.data, HNSCC_17_PBMC.data, HNSCC_18_PBMC.data, HNSCC_19_PBMC.data, HNSCC_20_PBMC.data, HNSCC_21_PBMC.data, HNSCC_22_PBMC.data, HNSCC_23_PBMC.data, HNSCC_24_PBMC.data, HNSCC_25_PBMC.data, HNSCC_26_PBMC.data)


HNSCC_TIL.metadata<-rbind(HNSCC_1_TIL.metadata, HNSCC_2_TIL.metadata, HNSCC_3_TIL.metadata, HNSCC_4_TIL.metadata, HNSCC_5_TIL.metadata, HNSCC_6_TIL.metadata, HNSCC_7_TIL.metadata, HNSCC_8_TIL.metadata, HNSCC_9_TIL.metadata, HNSCC_10_TIL.metadata, HNSCC_11_TIL.metadata, HNSCC_12_TIL.metadata, HNSCC_13_TIL.metadata, HNSCC_14_TIL.metadata, HNSCC_15_TIL.metadata, HNSCC_16_TIL.metadata, HNSCC_17_TIL.metadata, HNSCC_18_TIL.metadata, HNSCC_19_TIL.metadata, HNSCC_20_TIL.metadata, HNSCC_21_TIL.metadata, HNSCC_22_TIL.metadata, HNSCC_23_TIL.metadata, HNSCC_24_TIL.metadata, HNSCC_25_TIL.metadata, HNSCC_26_TIL.metadata)
colnames(HNSCC_TIL.metadata)<-c("barcodes","group","sample")
rownames(HNSCC_TIL.metadata)<-HNSCC_TIL.metadata[,1]


HNSCC_TIL.data<-cbind(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)


HNSCC.metadata<-rbind(HNSCC_PBMC.metadata,HNSCC_TIL.metadata)
HNSCC.data<-cbind(HNSCC_PBMC.data,HNSCC_TIL.data)

remove(HNSCC_1_PBMC.data, HNSCC_2_PBMC.data, HNSCC_3_PBMC.data, HNSCC_4_PBMC.data, HNSCC_5_PBMC.data, HNSCC_6_PBMC.data, HNSCC_7_PBMC.data, HNSCC_8_PBMC.data, HNSCC_9_PBMC.data, HNSCC_10_PBMC.data, HNSCC_11_PBMC.data, HNSCC_12_PBMC.data, HNSCC_13_PBMC.data, HNSCC_14_PBMC.data, HNSCC_15_PBMC.data, HNSCC_16_PBMC.data, HNSCC_17_PBMC.data, HNSCC_18_PBMC.data, HNSCC_19_PBMC.data, HNSCC_20_PBMC.data, HNSCC_21_PBMC.data, HNSCC_22_PBMC.data, HNSCC_23_PBMC.data, HNSCC_24_PBMC.data, HNSCC_25_PBMC.data, HNSCC_26_PBMC.data)
remove(HNSCC_1_PBMC.metadata, HNSCC_2_PBMC.metadata, HNSCC_3_PBMC.metadata, HNSCC_4_PBMC.metadata, HNSCC_5_PBMC.metadata, HNSCC_6_PBMC.metadata, HNSCC_7_PBMC.metadata, HNSCC_8_PBMC.metadata, HNSCC_9_PBMC.metadata, HNSCC_10_PBMC.metadata, HNSCC_11_PBMC.metadata, HNSCC_12_PBMC.metadata, HNSCC_13_PBMC.metadata, HNSCC_14_PBMC.metadata, HNSCC_15_PBMC.metadata, HNSCC_16_PBMC.metadata, HNSCC_17_PBMC.metadata, HNSCC_18_PBMC.metadata, HNSCC_19_PBMC.metadata, HNSCC_20_PBMC.metadata, HNSCC_21_PBMC.metadata, HNSCC_22_PBMC.metadata, HNSCC_23_PBMC.metadata, HNSCC_24_PBMC.metadata, HNSCC_25_PBMC.metadata, HNSCC_26_PBMC.metadata)
remove(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)
remove(HNSCC_1_TIL.metadata, HNSCC_2_TIL.metadata, HNSCC_3_TIL.metadata, HNSCC_4_TIL.metadata, HNSCC_5_TIL.metadata, HNSCC_6_TIL.metadata, HNSCC_7_TIL.metadata, HNSCC_8_TIL.metadata, HNSCC_9_TIL.metadata, HNSCC_10_TIL.metadata, HNSCC_11_TIL.metadata, HNSCC_12_TIL.metadata, HNSCC_13_TIL.metadata, HNSCC_14_TIL.metadata, HNSCC_15_TIL.metadata, HNSCC_16_TIL.metadata, HNSCC_17_TIL.metadata, HNSCC_18_TIL.metadata, HNSCC_19_TIL.metadata, HNSCC_20_TIL.metadata, HNSCC_21_TIL.metadata, HNSCC_22_TIL.metadata, HNSCC_23_TIL.metadata, HNSCC_24_TIL.metadata, HNSCC_25_TIL.metadata, HNSCC_26_TIL.metadata)
remove(a,b,c,d,g,h)
remove(HNSCC_PBMC.data,HNSCC_PBMC.metadata,HNSCC_TIL.data,HNSCC_TIL.metadata)

HNSCC <- CreateSeuratObject(counts = HNSCC.data, project = "HNSCC3k",meta.data = HNSCC.metadata,min.cells = 3, min.features = 200)
remove(HNSCC.data,HNSCC.metadata)
HNSCC[["percent.mt"]] <- PercentageFeatureSet(HNSCC, pattern = "^MT-")
VlnPlot(HNSCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by = "group")
HNSCC <- subset(HNSCC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50)
HNSCC <- NormalizeData(HNSCC)
HNSCC <- FindVariableFeatures(HNSCC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HNSCC)
HNSCC <- ScaleData(HNSCC, features = all.genes)
HNSCC <- RunPCA(HNSCC, features = VariableFeatures(object = HNSCC))
ElbowPlot(HNSCC)
HNSCC <- FindNeighbors(HNSCC, dims = 1:20)
HNSCC <- FindClusters(HNSCC, resolution = 0.8)
HNSCC <- RunUMAP(HNSCC, dims = 1:20)
HNSCC <- RunTSNE(HNSCC, dims = 1:20,check_duplicates = FALSE)
DimPlot(HNSCC, reduction = "umap")
DimPlot(HNSCC, reduction = "tsne")
FeaturePlot(HNSCC,features = "FOXP3")
Idents(HNSCC)<-HNSCC@meta.data$sample
DimPlot(HNSCC, reduction = "umap")
DimPlot(HNSCC, reduction = "tsne")
FeaturePlot(HNSCC,features = "FOXP3",reduction = "tsne")
