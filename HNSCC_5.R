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
  assign(paste("HNSCC_",i,"_PBMC.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))),rep("PBMC",length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = "")))))))
  a<-get(paste("HNSCC_",i,"_PBMC.metadata",sep = ""))
  b<-as.data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))
  colnames(a)<-c("barcode","group")
  rownames(a)<-b[,1]
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),CreateSeuratObject(counts = g, project = paste("HNSCC_",i,"_PBMC",sep = ""),meta.data = a))
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),NormalizeData(get(paste("HNSCC_",i,"_PBMC",sep = ""))))
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),FindVariableFeatures(get(paste("HNSCC_",i,"_PBMC",sep = "")), selection.method = "vst", nfeatures = 2000))
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),Read10X(paste("D:/Jmjd1c_Treg_Tumor/GSE139324/HNSCC_",i,"_TIL",sep = "")))
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),as.data.frame(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))
  h<-get(paste("HNSCC_",i,"_TIL.data",sep = ""))
  for (k in 1:length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))) {
    colnames(h)[k]<-paste(colnames(h[k]),"TIL",i,k,sep = "-")
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
colnames(HNSCC_PBMC.metadata)<-c("barcodes","group")
rownames(HNSCC_PBMC.metadata)<-HNSCC_PBMC.metadata[,1]

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
colnames(HNSCC_TIL.metadata)<-c("barcodes","group")
rownames(HNSCC_TIL.metadata)<-HNSCC_TIL.metadata[,1]

e<-paste("HNSCC_",1,"_TIL.data",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_TIL.data",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_TIL.data<-cbind(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)


HNSCC.metadata<-rbind(HNSCC_PBMC.metadata,HNSCC_TIL.metadata)
HNSCC.data<-cbind(HNSCC_PBMC.data,HNSCC_TIL.data)

remove(HNSCC_1_PBMC.data, HNSCC_2_PBMC.data, HNSCC_3_PBMC.data, HNSCC_4_PBMC.data, HNSCC_5_PBMC.data, HNSCC_6_PBMC.data, HNSCC_7_PBMC.data, HNSCC_8_PBMC.data, HNSCC_9_PBMC.data, HNSCC_10_PBMC.data, HNSCC_11_PBMC.data, HNSCC_12_PBMC.data, HNSCC_13_PBMC.data, HNSCC_14_PBMC.data, HNSCC_15_PBMC.data, HNSCC_16_PBMC.data, HNSCC_17_PBMC.data, HNSCC_18_PBMC.data, HNSCC_19_PBMC.data, HNSCC_20_PBMC.data, HNSCC_21_PBMC.data, HNSCC_22_PBMC.data, HNSCC_23_PBMC.data, HNSCC_24_PBMC.data, HNSCC_25_PBMC.data, HNSCC_26_PBMC.data)
remove(HNSCC_1_PBMC.metadata, HNSCC_2_PBMC.metadata, HNSCC_3_PBMC.metadata, HNSCC_4_PBMC.metadata, HNSCC_5_PBMC.metadata, HNSCC_6_PBMC.metadata, HNSCC_7_PBMC.metadata, HNSCC_8_PBMC.metadata, HNSCC_9_PBMC.metadata, HNSCC_10_PBMC.metadata, HNSCC_11_PBMC.metadata, HNSCC_12_PBMC.metadata, HNSCC_13_PBMC.metadata, HNSCC_14_PBMC.metadata, HNSCC_15_PBMC.metadata, HNSCC_16_PBMC.metadata, HNSCC_17_PBMC.metadata, HNSCC_18_PBMC.metadata, HNSCC_19_PBMC.metadata, HNSCC_20_PBMC.metadata, HNSCC_21_PBMC.metadata, HNSCC_22_PBMC.metadata, HNSCC_23_PBMC.metadata, HNSCC_24_PBMC.metadata, HNSCC_25_PBMC.metadata, HNSCC_26_PBMC.metadata)
remove(HNSCC_1_PBMC, HNSCC_2_PBMC, HNSCC_3_PBMC, HNSCC_4_PBMC, HNSCC_5_PBMC, HNSCC_6_PBMC, HNSCC_7_PBMC, HNSCC_8_PBMC, HNSCC_9_PBMC, HNSCC_10_PBMC, HNSCC_11_PBMC, HNSCC_12_PBMC, HNSCC_13_PBMC, HNSCC_14_PBMC, HNSCC_15_PBMC, HNSCC_16_PBMC, HNSCC_17_PBMC, HNSCC_18_PBMC, HNSCC_19_PBMC, HNSCC_20_PBMC, HNSCC_21_PBMC, HNSCC_22_PBMC, HNSCC_23_PBMC, HNSCC_24_PBMC, HNSCC_25_PBMC, HNSCC_26_PBMC)
remove(HNSCC_1_TIL, HNSCC_2_TIL, HNSCC_3_TIL, HNSCC_4_TIL, HNSCC_5_TIL, HNSCC_6_TIL, HNSCC_7_TIL, HNSCC_8_TIL, HNSCC_9_TIL, HNSCC_10_TIL, HNSCC_11_TIL, HNSCC_12_TIL, HNSCC_13_TIL, HNSCC_14_TIL, HNSCC_15_TIL, HNSCC_16_TIL, HNSCC_17_TIL, HNSCC_18_TIL, HNSCC_19_TIL, HNSCC_20_TIL, HNSCC_21_TIL, HNSCC_22_TIL, HNSCC_23_TIL, HNSCC_24_TIL, HNSCC_25_TIL, HNSCC_26_TIL)
remove(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)
remove(HNSCC_1_TIL.metadata, HNSCC_2_TIL.metadata, HNSCC_3_TIL.metadata, HNSCC_4_TIL.metadata, HNSCC_5_TIL.metadata, HNSCC_6_TIL.metadata, HNSCC_7_TIL.metadata, HNSCC_8_TIL.metadata, HNSCC_9_TIL.metadata, HNSCC_10_TIL.metadata, HNSCC_11_TIL.metadata, HNSCC_12_TIL.metadata, HNSCC_13_TIL.metadata, HNSCC_14_TIL.metadata, HNSCC_15_TIL.metadata, HNSCC_16_TIL.metadata, HNSCC_17_TIL.metadata, HNSCC_18_TIL.metadata, HNSCC_19_TIL.metadata, HNSCC_20_TIL.metadata, HNSCC_21_TIL.metadata, HNSCC_22_TIL.metadata, HNSCC_23_TIL.metadata, HNSCC_24_TIL.metadata, HNSCC_25_TIL.metadata, HNSCC_26_TIL.metadata)
remove(a,b,c,d,e,f,g,h,i,j,k)
remove(HNSCC_PBMC.data,HNSCC_PBMC.metadata,HNSCC_TIL.data,HNSCC_TIL.metadata)

for (i in 1:6) {
  assign(paste("HD_",i,"_PBMC.data",sep = ""),Read10X(paste("D:/Jmjd1c_Treg_Tumor/GSE139324/HD_",i,"_PBMC",sep = "")))
  assign(paste("HD_",i,"_PBMC.data",sep = ""),as.data.frame(get(paste("HD_",i,"_PBMC.data",sep = ""))))
  g<-get(paste("HD_",i,"_PBMC.data",sep = ""))
  for (j in 1:length(colnames(get(paste("HD_",i,"_PBMC.data",sep = ""))))) {
    colnames(g)[j]<-paste(colnames(g[j]),"HDPBMC",i,j,sep = "-")
  }
  assign(paste("HD_",i,"_PBMC.data",sep = ""),g)
  assign(paste("HD_",i,"_PBMC.metadata",sep = ""),data.frame(colnames(get(paste("HD_",i,"_PBMC.data",sep = ""))),rep("HD_PBMC",length(colnames(get(paste("HD_",i,"_PBMC.data",sep = "")))))))
  a<-get(paste("HD_",i,"_PBMC.metadata",sep = ""))
  b<-as.data.frame(colnames(get(paste("HD_",i,"_PBMC.data",sep = ""))))
  colnames(a)<-c("barcodes","group")
  rownames(a)<-b[,1]
  assign(paste("HD_",i,"_PBMC.metadata",sep = ""),a)
}

for (i in 1:5) {
  assign(paste("HD_",i,"_Tonsil.data",sep = ""),Read10X(paste("D:/Jmjd1c_Treg_Tumor/GSE139324/HD_",i,"_Tonsil",sep = "")))
  assign(paste("HD_",i,"_Tonsil.data",sep = ""),as.data.frame(get(paste("HD_",i,"_Tonsil.data",sep = ""))))
  h<-get(paste("HD_",i,"_Tonsil.data",sep = ""))
  for (k in 1:length(colnames(get(paste("HD_",i,"_Tonsil.data",sep = ""))))) {
    colnames(h)[k]<-paste(colnames(h[k]),"HDTonsil",i,k,sep = "-")
  }
  assign(paste("HD_",i,"_Tonsil.data",sep = ""),h)
  assign(paste("HD_",i,"_Tonsil.metadata",sep = ""),data.frame(colnames(get(paste("HD_",i,"_Tonsil.data",sep = ""))),rep("HD_Tonsil",length(colnames(get(paste("HD_",i,"_Tonsil.data",sep = "")))))))
  c<-get(paste("HD_",i,"_Tonsil.metadata",sep = ""))
  d<-as.data.frame(colnames(get(paste("HD_",i,"_Tonsil.data",sep = ""))))
  colnames(c)<-c("barcodes","group")
  rownames(c)<-d[,1]
  assign(paste("HD_",i,"_Tonsil.metadata",sep = ""),c)
}

HNSCC.data<-cbind(HNSCC.data,HD_1_PBMC.data,HD_2_PBMC.data,HD_3_PBMC.data,HD_4_PBMC.data,HD_5_PBMC.data,HD_6_PBMC.data,HD_1_Tonsil.data,HD_2_Tonsil.data,HD_3_Tonsil.data,HD_4_Tonsil.data,HD_5_Tonsil.data)
HNSCC.metadata<-rbind(HNSCC.metadata,HD_1_PBMC.metadata,HD_2_PBMC.metadata,HD_3_PBMC.metadata,HD_4_PBMC.metadata,HD_5_PBMC.metadata,HD_6_PBMC.metadata,HD_1_Tonsil.metadata,HD_2_Tonsil.metadata,HD_3_Tonsil.metadata,HD_4_Tonsil.metadata,HD_5_Tonsil.metadata)

HNSCC <- CreateSeuratObject(counts = HNSCC.data, project = "HNSCC3k",meta.data = HNSCC.metadata,min.cells = 3, min.features = 200)
HNSCC[["percent.mt"]] <- PercentageFeatureSet(HNSCC, pattern = "^MT-")
VlnPlot(HNSCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
HNSCC <- subset(HNSCC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50)
HNSCC <- NormalizeData(HNSCC)
HNSCC <- FindVariableFeatures(HNSCC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HNSCC)
remove(HNSCC.data,HD_1_PBMC.data,HD_2_PBMC.data,HD_3_PBMC.data,HD_4_PBMC.data,HD_5_PBMC.data,HD_6_PBMC.data,HD_1_Tonsil.data,HD_2_Tonsil.data,HD_3_Tonsil.data,HD_4_Tonsil.data,HD_5_Tonsil.data,HNSCC.metadata,HD_1_PBMC.metadata,HD_2_PBMC.metadata,HD_3_PBMC.metadata,HD_4_PBMC.metadata,HD_5_PBMC.metadata,HD_6_PBMC.metadata,HD_1_Tonsil.metadata,HD_2_Tonsil.metadata,HD_3_Tonsil.metadata,HD_4_Tonsil.metadata,HD_5_Tonsil.metadata)
remove(a,b,c,d,g,h,i,j,k)
HNSCC <- ScaleData(HNSCC, features = all.genes)
HNSCC <- RunPCA(HNSCC, features = VariableFeatures(object = HNSCC))
ElbowPlot(HNSCC)
HNSCC <- FindNeighbors(HNSCC, dims = 1:20)
HNSCC <- FindClusters(HNSCC, resolution = 0.8)
HNSCC <- RunUMAP(HNSCC, dims = 1:20)
HNSCC <- RunTSNE(HNSCC, dims = 1:20,check_duplicates = FALSE)
DimPlot(HNSCC, reduction = "umap")
DimPlot(HNSCC, reduction = "umap",split.by = "group")
DimPlot(HNSCC, reduction = "tsne")
DimPlot(HNSCC, reduction = "tsne",split.by = "group")

VlnPlot(HNSCC,features = c("FOXP3"),pt.size = 0)
FeaturePlot(HNSCC,features = c("FOXP3"))
Treg<-subset(HNSCC,idents = "4")
Idents(Treg)<-Treg@meta.data$group




library(ggpubr)
Treg_PBMC_TIL<-subset(Treg,idents = c("PBMC","TIL"))
VlnPlot(Treg_PBMC_TIL,features = "JMJD1C")
Treg_PBMC<-subset(Treg_PBMC_TIL,idents = "PBMC")
Treg_TIL<-subset(Treg_PBMC_TIL,idents = "TIL")
K<-as.data.frame(Treg_PBMC@assays$RNA@scale.data)
P<-as.data.frame(t(K))
Jmjd1c_Treg_PBMC<-select(P,starts_with("JMJD1C"))
L<-as.data.frame(Treg_TIL@assays$RNA@scale.data)
N<-as.data.frame(t(L))
Jmjd1c_Treg_TIL<-select(N,starts_with("JMJD1C"))
write.table(Jmjd1c_Treg_PBMC$JMJD1C,file = "D:/Jmjd1c_Treg_Tumor/GSE139324_code/Jmjd1c_Treg_PBMC.txt",sep = "\t")
write.table(Jmjd1c_Treg_TIL$JMJD1C,file = "D:/Jmjd1c_Treg_Tumor/GSE139324_code/Jmjd1c_Treg_TIL.txt",sep = "\t")


FeaturePlot(Treg_PBMC_TIL,features = "JMJD1C",split.by = "group",cols = c("#00008B","#FF69B4"),min.cutoff = "q10",max.cutoff = "q95",pt.size = 0.5,order = TRUE)
FeaturePlot(Treg_PBMC_TIL,features = "JMJD1C",split.by = "group",cols = c("#00008B","#FF69B4"),min.cutoff = "min",max.cutoff = "max",pt.size = 0.1,order = TRUE,reduction = "pca")
FeaturePlot(Treg_PBMC_TIL,features = "JMJD1C",cols = c("#00008B","#FF69B4"),min.cutoff = "q10",max.cutoff = "q95",pt.size = 0.1,order = TRUE,reduction = "umap")

DimPlot(Treg,group.by = "group")
