library(Seurat)
library(dplyr)

raw.data <- read.csv('/Users/yajin/Desktop/LuLab/cancer/topgene/TCGA.txt',,sep='\t',row.names = 1)
cancer <- CreateSeuratObject(counts = raw.data, min.cells = 3, min.features = 200)
# filtered samples detecting less than 200 genes
# filtered genes detecting less than 3 samples
cancer
# Visualize QC metrics as a violin plot
VlnPlot(cancer, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot_twoFeatures <- FeatureScatter(cancer, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_twoFeatures

cancer <- NormalizeData(cancer, normalization.method = "LogNormalize", scale.factor = 10000)

# find HVGs
cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cancer), 10)
plot1 <- VariableFeaturePlot(cancer,raster=FALSE,cols=c("#003366","#FF9900"))
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 #+ plot2

all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, features = VariableFeatures(object = cancer))
DimHeatmap(cancer, dims = 1:6, cells = 200, balanced = TRUE)

cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
#JackStrawPlot(cancer, dims = 1:15)
ElbowPlot(cancer)
cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.5)

cancer <- RunUMAP(cancer, dims = 1:10)
DimPlot(cancer, reduction = "umap",label = TRUE,pt.size=2)
