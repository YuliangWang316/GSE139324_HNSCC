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
  assign(paste("HNSCC_",i,"_PBMC.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))),rep(paste("PBMC",i,sep = ""),length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = "")))))))
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
  assign(paste("HNSCC_",i,"_TIL.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))),rep(paste("TIL",i,sep = ""),length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = "")))))))
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

HNSCC <- CreateSeuratObject(counts = HNSCC.data, project = "HNSCC3k",meta.data = HNSCC.metadata,min.cells = 3, min.features = 200)
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
DimPlot(HNSCC, reduction = "umap",split.by = "group")
DimPlot(HNSCC, reduction = "tsne")
DimPlot(HNSCC, reduction = "tsne",split.by = "group")

remove(HNSCC.data,HNSCC.metadata)
VlnPlot(HNSCC,features = c("FOXP3"),pt.size = 0)
FeaturePlot(HNSCC,features = c("FOXP3"))
Treg<-subset(HNSCC,idents = "9")
VlnPlot(HNSCC,features = "CD8B",pt.size = 0,sort = TRUE)
FeaturePlot(HNSCC,features = c("CD3G","CD8B"),label = TRUE)
CD8T<-subset(HNSCC,idents = c("3","5","7","10"))
#CD8T<-subset(HNSCC,idents = "5")
#remove(HNSCC,all.genes)

Idents(Treg)<-Treg@meta.data$group
Idents(CD8T)<-CD8T@meta.data$group
Idents(HNSCC)<-HNSCC@meta.data$seurat_clusters
new.cluster.ids <- c("0", "1", "2", "CD8_T", "4", "CD8_T",
                     "6", "CD8_T", "8","9","CD8_T","11","12","13","14","15","16","17","18","19","20","21","22","23")
names(new.cluster.ids) <- levels(HNSCC)
HNSCC <- RenameIdents(HNSCC, new.cluster.ids)
HNSCC_markers<-FindAllMarkers(HNSCC,only.pos = TRUE,logfc.threshold = 0.25,min.pct = 0.1)
Idents(HNSCC)<-HNSCC@meta.data$group
CD8T.markers<-filter(HNSCC_markers,p_val_adj<0.00001 & cluster == "CD8_T" & avg_log2FC)
#a<-read.table("D:/Jmjd1c_Treg_Tumor/GSE139324_code/geneset(2).txt",sep = "\t",header = TRUE)
#b<-read.table("D:/Jmjd1c_Treg_Tumor/GSE139324_code/geneset(3).txt",sep = "\t",header = TRUE)
#a_new<-a[2:301,]
#b_new<-b[2:10,]
#CD8T<-AddModuleScore(CD8T,features = a_new,name = "A")
#CD8T<-AddModuleScore(CD8T,features = b_new,name = "C")

#library(dplyr)
#remove(a,b,a_new,b_new)
a_new<-unique(CD8T.markers$gene)
HNSCC<-AddModuleScore(HNSCC,features = a_new,name = "A")
for (i in 1:26) {
  assign(paste("Treg_",i,"_TIL",sep = ""),subset(Treg,idents = paste("TIL",i,sep = "")))
  assign("a",slot(get(paste("Treg_",i,"_TIL",sep = "")),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@scale.data))))
  assign("c",select(b,one_of("JMJD1C")))
  #assign(paste("CD8T_",i,"_TIL",sep = ""),subset(CD8T,idents = paste("TIL",i,sep = "")))
  #assign("d",slot(get(paste("CD8_",i,"_TIL",sep = "")),"meta.data"))
  assign(paste("HNSCC_",i,"_TIL",sep = ""),subset(HNSCC,idents = paste("TIL",i,sep = "")))
  assign("d",slot(get(paste("HNSCC_",i,"_TIL",sep = "")),"meta.data"))
  #assign("d",slot(get(paste("HNSCC_",i,"_TIL",sep = "")),"assays"))
  #assign("e",as.data.frame(t(as.data.frame(d$RNA@scale.data))))
  #assign("f",select(b,one_of("CD8B")))
  A<-select(d,starts_with("A"))
  d$A<-apply(A,MARGIN = 1,median)
  #C<-select(d,starts_with("C"))
  #d$C<-apply(C,MARGIN = 1,median)
  assign("e",as.data.frame(d$A))
  colnames(e)<-"A"
  assign(paste("scaldata_",i,"_TIL",sep = ""),data.frame(sum(c$JMJD1C),sum(e$A)))
  #f<-(length(rownames(d)))/(length(rownames(e)))
  #assign(paste("scaldata_",i,"_TIL",sep = ""),data.frame(sum(c$JMJD1C),f))

}


rm(Treg_1_TIL,Treg_10_TIL,Treg_11_TIL,Treg_12_TIL,Treg_13_TIL,Treg_14_TIL,Treg_15_TIL,Treg_16_TIL,Treg_17_TIL,Treg_18_TIL,Treg_19_TIL,Treg_2_TIL,Treg_21_TIL,Treg_22_TIL,Treg_23_TIL,Treg_24_TIL,Treg_25_TIL,Treg_26_TIL,Treg_3_TIL,Treg_4_TIL,Treg_5_TIL,Treg_6_TIL,Treg_7_TIL,Treg_8_TIL,Treg_9_TIL,Treg_20_TIL,a,b,c,d,e)
scaldata<-rbind(scaldata_1_TIL,scaldata_10_TIL,scaldata_11_TIL,scaldata_12_TIL,scaldata_13_TIL,scaldata_14_TIL,scaldata_15_TIL,scaldata_16_TIL,scaldata_17_TIL,scaldata_18_TIL,scaldata_19_TIL,scaldata_2_TIL,scaldata_21_TIL,scaldata_22_TIL,scaldata_23_TIL,scaldata_24_TIL,scaldata_25_TIL,scaldata_26_TIL,scaldata_3_TIL,scaldata_4_TIL,scaldata_5_TIL,scaldata_6_TIL,scaldata_7_TIL,scaldata_8_TIL,scaldata_9_TIL,scaldata_20_TIL)
rm(scaldata_1_TIL,scaldata_10_TIL,scaldata_11_TIL,scaldata_12_TIL,scaldata_13_TIL,scaldata_14_TIL,scaldata_15_TIL,scaldata_16_TIL,scaldata_17_TIL,scaldata_18_TIL,scaldata_19_TIL,scaldata_2_TIL,scaldata_21_TIL,scaldata_22_TIL,scaldata_23_TIL,scaldata_24_TIL,scaldata_25_TIL,scaldata_26_TIL,scaldata_3_TIL,scaldata_4_TIL,scaldata_5_TIL,scaldata_6_TIL,scaldata_7_TIL,scaldata_8_TIL,scaldata_9_TIL,scaldata_20_TIL,STAT3,PI3K,i)
colnames(scaldata)<-c("JMJD1C","A")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=A)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 

#P<-NULL
#for (i in 1:86) {
#  S<-paste("d$STAT3",i,",",sep = "")
#  P<-paste(P,S,sep = "")
#}



Idents(CD8T)<-CD8T@meta.data$seurat_clusters
VlnPlot(CD8T,features = "CD8B",sort = TRUE)

write.table(scaldata,file = "D:/Jmjd1c_Treg_Tumor/GSE139324_code/1C_Treg_TotaldividedbyCD8.txt",sep = "\t")

