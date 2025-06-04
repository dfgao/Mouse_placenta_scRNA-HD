
# 01.load data It's best to use one resolution----
library(arrow)
library(hdf5r)
library(spacexr)
library(Banksy)
library(SeuratWrappers)
library(viridis)
jet.colors <- colorRampPalette(c("blue", "cyan", "green", "yellow",'red2'))
jet.palette <- jet.colors(50)

localdir <- "/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A3/outs/binned_outputs/"
hd.a3.seu.4res <- Load10X_Spatial(data.dir = '/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A3/outs/', bin.size = c(8,16,20,30))
hd.a4.seu.4res <- Load10X_Spatial(data.dir = '/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A4/outs/', bin.size = c(8,16,20,30))

Assays(hd.a3.seu.4res)
DefaultAssay(hd.a3.seu.4res) <- "Spatial.008um"

VlnPlot(hd.a3.seu.4res, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
SpatialFeaturePlot(hd.a3.seu.4res, features = "nCount_Spatial.008um") + theme(legend.position = "right")
SpatialFeaturePlot(hd.a3.seu.4res, features = "nFeature_Spatial.008um",crop = T,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=300, limit=c(-0, 600), space="Lab")
  # scale_fill_gradientn(colors = jet.colors(5),limit = c(0,500))

Assays(hd.a4.seu.4res)
DefaultAssay(hd.a4.seu.4res) <- "Spatial.008um"

VlnPlot(hd.a4.seu.4res, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
SpatialFeaturePlot(hd.a4.seu.4res, features = "nCount_Spatial.008um") + theme(legend.position = "right")
SpatialFeaturePlot(hd.a4.seu.4res, features = "nFeature_Spatial.008um",crop = T,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=300, limit=c(-0, 600), space="Lab")

# 02.split image ----
test <- Load10X_Spatial(data.dir = '/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A4/outs/', bin.size = 8)
summary(test$nCount_Spatial.008um)
coords <- test@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- test@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x

ggplot( coords,aes(x = y, y = x)) + 
  geom_point(size = .002,color = 'red') + theme_bw()

test$x <- coords$y
test$y <- coords$x
test$REGION <- ifelse(test$x < 14000, 'Age4', 'KO2')

SpatialDimPlot(test, group.by = "REGION", label = TRUE)
test_left <- subset(test, subset = REGION == "Age4")
SpatialDimPlot(test_left, group.by = "REGION", label = TRUE)

rm(test,coords,test_left)

## A3
DefaultAssay(hd.a3.seu.4res) <- "Spatial.008um"
coords <- hd.a3.seu.4res@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- hd.a3.seu.4res@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x
ggplot( coords,aes(x = y, y = x)) + 
  geom_point(size = .002,color = 'red') + theme_bw()

hd.a3.seu.4res$coord_spatial_0.08_x <- coords$y
hd.a3.seu.4res$coord_spatial_0.08_y <- coords$x
hd.a3.seu.4res$region <- ifelse(hd.a3.seu.4res$coord_spatial_0.08_x < 13000,'Young1', 'Age3')
SpatialDimPlot(hd.a3.seu.4res, group.by = "region", label = T)
Age3.hd <- subset(hd.a3.seu.4res, subset = region == "Age3")
Young1.hd <- subset(hd.a3.seu.4res, subset = region == "Young1")

p1 <- SpatialFeaturePlot(Age3.hd, features = "nCount_Spatial.008um",crop = F,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=300, limit=c(-0, 600), space="Lab",
                       name="Age3.nUMIs.8um")
p2 <- SpatialFeaturePlot(Young1.hd, features = "nCount_Spatial.008um",crop = F,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=300, limit=c(-0, 600), space="Lab",
                       name="Young1.nUMIs.8um")

## A4
DefaultAssay(hd.a4.seu.4res) <- "Spatial.008um"
coords <- hd.a4.seu.4res@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- hd.a4.seu.4res@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x
ggplot( coords,aes(x = y, y = x)) + 
  geom_point(size = .001,color = 'red') + theme_bw()

hd.a4.seu.4res$coord_spatial_0.08_x <- coords$y
hd.a4.seu.4res$coord_spatial_0.08_y <- coords$x
hd.a4.seu.4res$region <- ifelse(hd.a4.seu.4res$coord_spatial_0.08_x < 14000,'Age4', 'KO2')
SpatialDimPlot(hd.a4.seu.4res, group.by = "region", label = T)
Age4.hd <- subset(hd.a4.seu.4res, subset = region == "Age4")
KO2.hd <- subset(hd.a4.seu.4res, subset = region == "KO2")

p3 <- SpatialFeaturePlot(KO2.hd, features = "nCount_Spatial.008um",crop = F,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=300, limit=c(-0, 600), space="Lab",
                       name="KO2.nUMIs.8um")
p4 <- SpatialFeaturePlot(Age4.hd, features = "nCount_Spatial.008um",crop = F,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=300, limit=c(-0, 600), space="Lab",
                       name="Age4.nUMIs.8um")

VlnPlot(Age4.hd, features = "nCount_Spatial.008um", pt.size = 0) +  NoLegend() +   ggtitle('Age4.nUMIs.8um') + ylim(0,1000) 
VlnPlot(KO2.hd, features = "nCount_Spatial.008um", pt.size = 0) + NoLegend() + ggtitle('KO2.nUMIs.8um')+ ylim(0,1000)
VlnPlot(Young1.hd, features = "nCount_Spatial.008um", pt.size = 0) + NoLegend() + ggtitle('Young1.nUMIs.8um')+ ylim(0,1000)
VlnPlot(Age3.hd, features = "nCount_Spatial.008um", pt.size = 0) + NoLegend() + ggtitle('Age3.nUMIs.8um')+ ylim(0,1000)
summary(na.omit(Age4.hd$nCount_Spatial.008um))
summary(na.omit(KO2.hd$nCount_Spatial.008um))
summary(na.omit(Young1.hd$nCount_Spatial.008um))
summary(na.omit(Age3.hd$nCount_Spatial.008um))

A3.meta <- hd.a3.seu.4res@meta.data
A4.meta <- hd.a4.seu.4res@meta.data

# 2um bin A3
hd.a3.2um <- Load10X_Spatial(data.dir = '/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A3/outs/', bin.size = 2)
coords <- hd.a3.2um@images[["slice1.002um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- hd.a3.2um@images[["slice1.002um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x
ggplot( coords,aes(x = y, y = x)) +
  geom_point(size = .002,color = 'red') + theme_bw()

hd.a3.2um$coord_spatial_0.02_x <- coords$y
hd.a3.2um$coord_spatial_0.02_y <- coords$x
hd.a3.2um$region <- ifelse(hd.a3.2um$coord_spatial_0.02_x < 13000,'Young1', 'Age3')
SpatialDimPlot(hd.a3.2um, group.by = "region", label = T)
Age3.2um.hd <- subset(hd.a3.2um, subset = region == "Age3")
Young1.2um.hd <- subset(hd.a3.2um, subset = region == "Young1")

p1 <- SpatialFeaturePlot(Age3.2um.hd, features = "nCount_Spatial.002um",crop = F,image.alpha = '.2') +
  theme(legend.position = "right") +
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=30, limit=c(0, 60), space="Lab",
                       name="Age3.nUMIs.2um")
p2 <- SpatialFeaturePlot(Young1.2um.hd, features = "nCount_Spatial.002um",crop = F,image.alpha = '.2') +
  theme(legend.position = "right") +
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=30, limit=c(-0, 60), space="Lab",
                       name="Young1.nUMIs.2um")

# 2um bin A4
hd.a4.2um <- Load10X_Spatial(data.dir = '/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/A4/outs/', bin.size = 2)
coords <- hd.a4.2um@images[["slice1.002um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- hd.a4.2um@images[["slice1.002um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x
ggplot( coords,aes(x = y, y = x)) + 
  geom_point(size = .001,color = 'red') + theme_bw()

hd.a4.2um$coord_spatial_0.02_x <- coords$y
hd.a4.2um$coord_spatial_0.02_y <- coords$x
hd.a4.2um$region <- ifelse(hd.a4.2um$coord_spatial_0.02_x < 14000,'Age4', 'KO2')
SpatialDimPlot(hd.a4.2um, group.by = "region", label = T)
Age4.2um.hd <- subset(hd.a4.2um, subset = region == "Age4")
KO2.2um.hd <- subset(hd.a4.2um, subset = region == "KO2")

p3 <- SpatialFeaturePlot(KO2.2um.hd, features = "nCount_Spatial.002um",crop = F,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=30, limit=c(-0, 60), space="Lab",
                       name="KO2.nUMIs.2um")
p4 <- SpatialFeaturePlot(Age4.2um.hd, features = "nCount_Spatial.002um",crop = F,image.alpha = '.2') + 
  theme(legend.position = "right") + 
  scale_fill_gradient2(low="blue", high="red", mid="yellow",
                       midpoint=30, limit=c(-0, 60), space="Lab",
                       name="Age4.nUMIs.2um")

VlnPlot(Age4.2um.hd, features = "nCount_Spatial.002um", pt.size = 0) +  NoLegend() +   ggtitle('Age4.nUMIs.2um') + ylim(0,60) 
VlnPlot(KO2.2um.hd, features = "nCount_Spatial.002um", pt.size = 0) + NoLegend() + ggtitle('KO2.nUMIs.2um')+ ylim(0,60)
VlnPlot(Young1.2um.hd, features = "nCount_Spatial.002um", pt.size = 0) + NoLegend() + ggtitle('Young1.nUMIs.2um')+ ylim(0,60)
VlnPlot(Age3.2um.hd, features = "nCount_Spatial.002um", pt.size = 0) + NoLegend() + ggtitle('Age3.nUMIs.2um')+ ylim(0,60)
summary(na.omit(Age4.2um.hd$nCount_Spatial.002um))
summary(na.omit(KO2.2um.hd$nCount_Spatial.002um))
summary(na.omit(Young1.2um.hd$nCount_Spatial.002um))
summary(na.omit(Age3.2um.hd$nCount_Spatial.002um))

# 03.define tissue domains by key markers ----
plan(multisession, workers=1)
DefaultAssay(Age3.hd) <- 'Spatial.008um'
Age3.hd <- NormalizeData(Age3.hd)

DefaultAssay(Young1.hd) <- 'Spatial.008um'
Young1.hd <- NormalizeData(Young1.hd)

DefaultAssay(KO2.hd) <- 'Spatial.008um'
KO2.hd <- NormalizeData(KO2.hd)

## LZ
SeuratPipe::spatial_feature_plot(Age3.hd, features = c('Lepr','Stra6'),col_pal = 'BrBG',alpha = 0.9,) + ggtitle("LZ markers (8um)")
SpatialFeaturePlot(Age3.hd, features = c('Lepr','Stra6'),image.alpha = 0.1)
SpatialFeaturePlot(Young1.hd, features = c('Lepr','Stra6'),image.alpha = 0.1)
SpatialFeaturePlot(KO2.hd, features = c('Lepr','Stra6'),image.alpha = 0.1)

## JZ
SpatialFeaturePlot(Age3.hd, features = c('Prl8a9','Prl3a1'),image.alpha = 0.1)
SpatialFeaturePlot(Young1.hd, features = c('Prl8a9','Prl3a1'),image.alpha = 0.1)
SpatialFeaturePlot(KO2.hd, features = c('Prl8a9','Prl3a1'),image.alpha = 0.1)

## DB
SpatialFeaturePlot(Age3.hd, features = c('Inhba','Cryab'),image.alpha = 0.1)
SpatialFeaturePlot(Young1.hd, features = c('Inhba','Cryab'),image.alpha = 0.1)
SpatialFeaturePlot(KO2.hd, features = c('Inhba','Cryab'),image.alpha = 0.1,)

# 04. define tissue domain by baskey ----

## Age3
plan(multisession, workers=30)
Age3.hd <- FindVariableFeatures(Age3.hd)
Age3.hd <- RunBanksy(Age3.hd,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Spatial.008um", slot = "data", features = "variable",
                    k_geom = 30)
DefaultAssay(Age3.hd) <- "BANKSY"
Age3.hd <- RunPCA(Age3.hd, assay = "BANKSY", 
                  reduction.name = "pca.banksy", 
                  features = rownames(Age3.hd), npcs = 30) %>% 
  FindNeighbors( reduction = "pca.banksy", dims = 1:30) %>% 
  FindClusters( cluster.name = "banksy_cluster", resolution = 0.5)
Idents(Age3.hd) <- "banksy_cluster"
SpatialDimPlot(Age3.hd, group.by = "banksy_cluster",
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(use.cols,npg.cols))

## Young1 
Young1.hd <- FindVariableFeatures(Young1.hd)
Young1.hd <- RunBanksy(Young1.hd,
                     lambda = 0.8, verbose = TRUE,
                     assay = "Spatial.008um", slot = "data", features = "variable",
                     k_geom = 30)
DefaultAssay(Young1.hd) <- "BANKSY"
Young1.hd <- RunPCA(Young1.hd, assay = "BANKSY", 
                  reduction.name = "pca.banksy", 
                  features = rownames(Young1.hd), npcs = 30) %>% 
  FindNeighbors( reduction = "pca.banksy", dims = 1:30) %>% 
  FindClusters( cluster.name = "banksy_cluster", resolution = 0.5)
Idents(Young1.hd) <- "banksy_cluster"
SpatialDimPlot(Young1.hd, group.by = "banksy_cluster",
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(use.cols,npg.cols))

## KO2
KO2.hd <- FindVariableFeatures(KO2.hd)
KO2.hd <- RunBanksy(KO2.hd,
                     lambda = 0.8, verbose = TRUE,
                     assay = "Spatial.008um", slot = "data", features = "variable",
                     k_geom = 30)
DefaultAssay(KO2.hd) <- "BANKSY"
KO2.hd <- RunPCA(KO2.hd, assay = "BANKSY", 
                  reduction.name = "pca.banksy", 
                  features = rownames(KO2.hd), npcs = 30) %>% 
  FindNeighbors( reduction = "pca.banksy", dims = 1:30) %>% 
  FindClusters( cluster.name = "banksy_cluster", resolution = 0.5)
Idents(KO2.hd) <- "banksy_cluster"
SpatialDimPlot(KO2.hd, group.by = "banksy_cluster",
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(use.cols,npg.cols))

save.image(compress = F)

# 05.define tissue domain by SNN ----

## Age3
DefaultAssay(Age3.hd) <- "Spatial.008um"
Age3.hd <- ScaleData(Age3.hd)
Age3.hd <- SketchData(
  object = Age3.hd,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(Age3.hd) <- "sketch"
Age3.hd <- FindVariableFeatures(Age3.hd) %>% 
  ScaleData() %>% 
  RunPCA( assay = "sketch", reduction.name = "pca.sketch") %>% 
  FindNeighbors( assay = "sketch", reduction = "pca.sketch", dims = 1:30) %>% 
  FindClusters( cluster.name = "seurat_cluster.sketched", resolution = 1) %>% 
  RunUMAP( reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:30)
Age3.hd <- ProjectData(
  object = Age3.hd,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:30,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(Age3.hd) <- "sketch"
Idents(Age3.hd) <- "seurat_cluster.sketched"
p1 <- DimPlot(Age3.hd, reduction = "umap.sketch", label = F) +
  ggtitle("Sketched clustering (50,000 cells)") + 
  theme(legend.position = "bottom")

DefaultAssay(Age3.hd) <- "Spatial.008um"
Idents(Age3.hd) <- "seurat_cluster.projected"
p2 <- DimPlot(Age3.hd, reduction = "full.umap.sketch", label = F) + 
  ggtitle("Projected clustering (full dataset)") + 
  theme(legend.position = "bottom")

p1 | p2

SpatialDimPlot(Age3.hd, 
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(use.cols,npg.cols))

## Young1
DefaultAssay(Young1.hd) <- "Spatial.008um"
Young1.hd <- ScaleData(Young1.hd)
Young1.hd <- SketchData(
  object = Young1.hd,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(Young1.hd) <- "sketch"
Young1.hd <- FindVariableFeatures(Young1.hd) %>% 
  ScaleData() %>% 
  RunPCA( assay = "sketch", reduction.name = "pca.sketch") %>% 
  FindNeighbors( assay = "sketch", reduction = "pca.sketch", dims = 1:30) %>% 
  FindClusters( cluster.name = "seurat_cluster.sketched", resolution = 1) %>% 
  RunUMAP( reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:30)
Young1.hd <- ProjectData(
  object = Young1.hd,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:30,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(Young1.hd) <- "sketch"
Idents(Young1.hd) <- "seurat_cluster.sketched"
p1 <- DimPlot(Young1.hd, reduction = "umap.sketch", label = F) +
  ggtitle("Sketched clustering (50,000 cells)") + 
  theme(legend.position = "bottom")

DefaultAssay(Young1.hd) <- "Spatial.008um"
Idents(Young1.hd) <- "seurat_cluster.projected"
p2 <- DimPlot(Young1.hd, reduction = "full.umap.sketch", label = F) + 
  ggtitle("Projected clustering (full dataset)") + 
  theme(legend.position = "bottom")

p1 | p2

SpatialDimPlot(Young1.hd, 
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(use.cols,npg.cols))

## KO2
DefaultAssay(KO2.hd) <- "Spatial.008um"
KO2.hd <- ScaleData(KO2.hd)
KO2.hd <- SketchData(
  object = KO2.hd,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(KO2.hd) <- "sketch"
KO2.hd <- FindVariableFeatures(KO2.hd) %>% 
  ScaleData() %>% 
  RunPCA( assay = "sketch", reduction.name = "pca.sketch") %>% 
  FindNeighbors( assay = "sketch", reduction = "pca.sketch", dims = 1:30) %>% 
  FindClusters( cluster.name = "seurat_cluster.sketched", resolution = 0.2) %>% 
  RunUMAP( reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:30)
KO2.hd <- ProjectData(
  object = KO2.hd,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:30,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(KO2.hd) <- "sketch"
Idents(KO2.hd) <- "seurat_cluster.sketched"
p1 <- DimPlot(KO2.hd, reduction = "umap.sketch", label = F) +
  ggtitle("Sketched clustering (50,000 cells)") + 
  theme(legend.position = "bottom")

DefaultAssay(KO2.hd) <- "Spatial.008um"
Idents(KO2.hd) <- "seurat_cluster.projected"
p2 <- DimPlot(KO2.hd, reduction = "full.umap.sketch", label = F) + 
  ggtitle("Projected clustering (full dataset)") + 
  theme(legend.position = "bottom")

p1 | p2

SpatialDimPlot(KO2.hd, 
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) 

# 06.split JZ&DB -----
## handle select by loupe browser
Age3.jz_db.barcodes <- read.csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/age3_jz_db.csv',header = T)
Young1.jz_db.barcodes <- read.csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/young1_jz_db.csv',header = T)
ko2.jz_db.barcodes <- read.csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/05.10x.hd/02.align/ko2_jz_db.csv',header = T)
### Age3
DefaultAssay(Age3.hd) <- "Spatial.008um"
Age3.jz_db.hd <- subset(Age3.hd, cells = Age3.jz_db.barcodes$Barcode)
SpatialDimPlot(Age3.jz_db.hd, group.by = 'region',
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(use.cols,npg.cols))

### Young1
DefaultAssay(Young1.hd) <- "Spatial.008um"
Young1.jz_db.hd <- subset(Young1.hd, cells = Young1.jz_db.barcodes$Barcode)
SpatialDimPlot(Young1.jz_db.hd, group.by = 'region',
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(npg.cols))

### KO2
DefaultAssay(KO2.hd) <- "Spatial.008um"
KO2.jz_db.hd <- subset(KO2.hd, cells = ko2.jz_db.barcodes$Barcode)
SpatialDimPlot(KO2.jz_db.hd, group.by = 'region',
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(npg.cols[2]))

# 07.deconvolution ----

## CARD
library(CARD)

### Age3
coords <- Age3.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- Age3.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x
ggplot( coords,aes(x = y, y = x)) + 
  geom_point(size = .0001,color = 'red') + theme_bw()

Age3_CARD_obj = createCARDObject(
  sc_count = Clean_sct.inte.rm.lpt@assays$RNA$counts,
  sc_meta = Clean_sct.inte.rm.lpt@meta.data,
  spatial_count = Age3.jz_db.hd@assays$Spatial.008um$counts,
  spatial_location = coords,
  ct.varname = "cell.type.percise.new",
  ct.select = unique(Clean_sct.inte.rm.lpt$cell.type.percise.new),
  sample.varname = 'sample',
  minCountGene = 50,
  minCountSpot = 5) 

Age3_CARD_obj = CARD_deconvolution(CARD_object = Age3_CARD_obj)
Age3_CARD_obj@Proportion_CARD

ct.visualize = c( unique(Clean_sct.inte.rm.lpt$cell.type.percise.new))
coords <- Age3.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- Age3.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
coords <- coords[rownames(Age3_CARD_obj@spatial_location),]

CARD.visualize.prop(
  proportion = Age3_CARD_obj@Proportion_CARD,
  # spatial_location = Age3_CARD_obj@spatial_location, 
  spatial_location = coords, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4,                                 ### number of columns in the figure panel
  pointSize = 0.0001) 

age3.card.socre <- Age3_CARD_obj@Proportion_CARD
threshold <- 0.5
greater_than_thre <- age3.card.socre > threshold
age3.card.ct <- colSums(greater_than_thre)


### Young1
coords <- Young1.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- Young1.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x
ggplot( coords,aes(x = y, y = x)) + 
  geom_point(size = .0001,color = 'red') + theme_bw()

Young1_CARD_obj = createCARDObject(
  sc_count = Clean_sct.inte.rm.lpt@assays$RNA$counts,
  sc_meta = Clean_sct.inte.rm.lpt@meta.data,
  spatial_count = Young1.jz_db.hd@assays$Spatial.008um$counts,
  spatial_location = coords,
  ct.varname = "cell.type.percise.new",
  ct.select = unique(Clean_sct.inte.rm.lpt$cell.type.percise.new),
  sample.varname = 'sample',
  minCountGene = 50,
  minCountSpot = 5) 

Young1_CARD_obj = CARD_deconvolution(CARD_object = Young1_CARD_obj)
# Young1_CARD_obj@Proportion_CARD

ct.visualize = c( unique(Clean_sct.inte.rm.lpt$cell.type.percise.new))
# coords <- Young1.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
# rownames(coords) <- Young1.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
# coords <- coords[rownames(Young1_CARD_obj@spatial_location),]

CARD.visualize.prop(
  proportion = Young1_CARD_obj@Proportion_CARD,        
  spatial_location = Young1_CARD_obj@spatial_location,
  # spatial_location = coords, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4,                                 ### number of columns in the figure panel
  pointSize = 0.0001)

young1.card.socre <- Young1_CARD_obj@Proportion_CARD
threshold <- 0.5
greater_than_thre <- young1.card.socre > threshold
young1.card.ct <- colSums(greater_than_thre)
young1.card.ct/74892

### KO2
coords <- KO2.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
rownames(coords) <- KO2.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
coords$x <- -coords$x
ggplot( coords,aes(x = y, y = x)) + 
  geom_point(size = .0001,color = 'red') + theme_bw()

KO2_CARD_obj = createCARDObject(
  sc_count = Clean_sct.inte.rm.lpt@assays$RNA$counts,
  sc_meta = Clean_sct.inte.rm.lpt@meta.data,
  spatial_count = KO2.jz_db.hd@assays$Spatial.008um$counts,
  spatial_location = coords,
  ct.varname = "cell.type.percise.new",
  ct.select = unique(Clean_sct.inte.rm.lpt$cell.type.percise.new),
  sample.varname = 'sample',
  minCountGene = 50,
  minCountSpot = 5) 

KO2_CARD_obj = CARD_deconvolution(CARD_object = KO2_CARD_obj)
# KO2_CARD_obj@Proportion_CARD

ct.visualize = c( unique(Clean_sct.inte.rm.lpt$cell.type.percise.new))
# coords <- KO2.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@coords %>% as.data.frame()
# rownames(coords) <- KO2.jz_db.hd@images[["slice1.008um"]]@boundaries[["centroids"]]@cells
# coords <- coords[rownames(KO2_CARD_obj@spatial_location),]

CARD.visualize.prop(
  proportion = KO2_CARD_obj@Proportion_CARD,        
  spatial_location = KO2_CARD_obj@spatial_location,
  # spatial_location = coords, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4,                                 ### number of columns in the figure panel
  pointSize = 0.0001) 

ko2.card.socre <- KO2_CARD_obj@Proportion_CARD
threshold <- 0.5
greater_than_thre <- ko2.card.socre > threshold
ko2.card.ct <- colSums(greater_than_thre)
ko2.card.ct/83200

## RCTD
# library(doSNOW)
# library(foreach)

ref <- Clean_sct.inte.rm.lpt
Idents(ref)
counts <- ref@assays[["RNA"]]@layers[["counts"]]
counts <- ref[["RNA"]]
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 30000)

### Age3
DefaultAssay(Age3.jz_db.hd) <- "Spatial.008um"
Age3.jz_db.hd <- FindVariableFeatures(Age3.jz_db.hd)
Age3.jz_db.hd <- SketchData(
  object = Age3.jz_db.hd,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(Age3.jz_db.hd) <- "sketch"
Age3.jz_db.hd <- ScaleData(Age3.jz_db.hd) %>% 
  RunPCA( assay = "sketch", reduction.name = "pca.Age3.jz_db.hd.sketch", verbose = T) %>% 
  FindNeighbors( reduction = "pca.Age3.jz_db.hd.sketch", dims = 1:50) %>% 
  RunUMAP( reduction = "pca.Age3.jz_db.hd.sketch", reduction.name = "umap.Age3.jz_db.hd.sketch", return.model = T, dims = 1:30, verbose = T)
DimPlot(Age3.jz_db.hd)
FeaturePlot(Age3.jz_db.hd,features = 'Ptprc',reduction= "umap.Age3.jz_db.hd.sketch")

counts_hd <- Age3.jz_db.hd[["sketch"]]$counts
Age3.jz_db.hd_cells_hd <- colnames(Age3.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Age3.jz_db.hd)[Age3.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 30)
Age3_RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
Age3.jz_db.hd.back <- Age3.jz_db.hd
Age3.jz_db.hd <- Age3.jz_db.hd.back
### add meta 
Age3.jz_db.hd <- AddMetaData(Age3.jz_db.hd, metadata = Age3_RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
Age3.jz_db.hd$first_type <- as.character(Age3.jz_db.hd$first_type)
Age3.jz_db.hd$first_type[is.na(Age3.jz_db.hd$first_type)] <- "Unknown"
Age3.jz_db.hd <- ProjectData(
  object = Age3.jz_db.hd,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Age3.jz_db.hd.sketch",
  umap.model = "umap.Age3.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Age3.jz_db.hd) <- "Spatial.008um"
Age3.jz_db.hd[[]][,"full_first_type"]

Age3.jz_db.hd$full_first_type <- factor(Age3.jz_db.hd$full_first_type, levels = c('GCs','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells','Unknown'))
SpatialDimPlot(Age3.jz_db.hd,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()
# Idents(Age3.jz_db.hd) <- "full_first_type"
# cells <- CellsByIdentities(Age3.jz_db.hd)
# excitatory_names <- c('GCs','Syns','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells')
SpatialDimPlot(Age3.jz_db.hd, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(Age3.jz_db.hd, cells.highlight = cells[c('Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 1,image.alpha = 0.3)

### Young1
DefaultAssay(Young1.jz_db.hd) <- "Spatial.008um"
Young1.jz_db.hd <- FindVariableFeatures(Young1.jz_db.hd)
Young1.jz_db.hd <- SketchData(
  object = Young1.jz_db.hd,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(Young1.jz_db.hd) <- "sketch"
Young1.jz_db.hd <- ScaleData(Young1.jz_db.hd) %>% 
  RunPCA( assay = "sketch", reduction.name = "pca.Young1.jz_db.hd.sketch", verbose = T) %>% 
  FindNeighbors( reduction = "pca.Young1.jz_db.hd.sketch", dims = 1:50) %>% 
  RunUMAP( reduction = "pca.Young1.jz_db.hd.sketch", reduction.name = "umap.Young1.jz_db.hd.sketch", return.model = T, dims = 1:30, verbose = T)
DimPlot(Young1.jz_db.hd,reduction= "umap.Young1.jz_db.hd.sketch")
FeaturePlot(Young1.jz_db.hd,features = 'Ptprc',reduction= "umap.Young1.jz_db.hd.sketch")

counts_hd <- Young1.jz_db.hd[["sketch"]]$counts
Young1.jz_db.hd_cells_hd <- colnames(Young1.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Young1.jz_db.hd)[Young1.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 30)
Young1_RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
Young1.jz_db.hd.back <- Young1.jz_db.hd


### add meta 
Young1.jz_db.hd <- AddMetaData(Young1.jz_db.hd, metadata = Young1_RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
Young1.jz_db.hd$first_type <- as.character(Young1.jz_db.hd$first_type)
Young1.jz_db.hd$first_type[is.na(Young1.jz_db.hd$first_type)] <- "Unknown"
Young1.jz_db.hd <- ProjectData(
  object = Young1.jz_db.hd,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Young1.jz_db.hd.sketch",
  umap.model = "umap.Young1.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Young1.jz_db.hd) <- "Spatial.008um"
Young1.jz_db.hd[[]][,"full_first_type"]
table(Young1.jz_db.hd$full_first_type)

SpatialDimPlot(Young1.jz_db.hd, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(Young1.jz_db.hd, cells.highlight = cells[c('Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 1,image.alpha = 0.3)


### KO2
DefaultAssay(KO2.jz_db.hd) <- "Spatial.008um"
KO2.jz_db.hd <- FindVariableFeatures(KO2.jz_db.hd)
KO2.jz_db.hd <- SketchData(
  object = KO2.jz_db.hd,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(KO2.jz_db.hd) <- "sketch"
KO2.jz_db.hd <- ScaleData(KO2.jz_db.hd) %>% 
  RunPCA( assay = "sketch", reduction.name = "pca.KO2.jz_db.hd.sketch", verbose = T) %>% 
  FindNeighbors( reduction = "pca.KO2.jz_db.hd.sketch", dims = 1:50) %>% 
  RunUMAP( reduction = "pca.KO2.jz_db.hd.sketch", reduction.name = "umap.KO2.jz_db.hd.sketch", return.model = T, dims = 1:30, verbose = T)
DimPlot(KO2.jz_db.hd,reduction= "umap.KO2.jz_db.hd.sketch")
FeaturePlot(KO2.jz_db.hd,features = 'Ptprc',reduction= "umap.KO2.jz_db.hd.sketch")

counts_hd <- KO2.jz_db.hd[["sketch"]]$counts
KO2.jz_db.hd_cells_hd <- colnames(KO2.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd)[KO2.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 30)
KO2_RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
KO2.jz_db.hd.back <- KO2.jz_db.hd

### add meta 
KO2.jz_db.hd <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd$first_type <- as.character(KO2.jz_db.hd$first_type)
KO2.jz_db.hd$first_type[is.na(KO2.jz_db.hd$first_type)] <- "Unknown"
KO2.jz_db.hd <- ProjectData(
  object = KO2.jz_db.hd,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(KO2.jz_db.hd) <- "Spatial.008um"
KO2.jz_db.hd[[]][,"full_first_type"]
table(KO2.jz_db.hd$full_first_type)
Idents(KO2.jz_db.hd) <- 'full_first_type'

SpatialDimPlot(KO2.jz_db.hd, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(KO2.jz_db.hd, cells.highlight = cells[c('Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 1,image.alpha = 0.3)

# 07.1 RCTD NOT USE SKETCH ----
Age3.jz_db.hd.ns <- subset(Age3.hd, cells = Age3.jz_db.barcodes$Barcode)
SpatialDimPlot(Age3.jz_db.hd.ns, group.by = 'region',
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(use.cols,npg.cols))

### Young1
DefaultAssay(Young1.hd) <- "Spatial.008um"
Young1.jz_db.hd.ns <- subset(Young1.hd, cells = Young1.jz_db.barcodes$Barcode)
SpatialDimPlot(Young1.jz_db.hd.ns, group.by = 'region',
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(npg.cols))

### KO2
DefaultAssay(KO2.hd) <- "Spatial.008um"
KO2.jz_db.hd.ns <- subset(KO2.hd, cells = ko2.jz_db.barcodes$Barcode)
SpatialDimPlot(KO2.jz_db.hd.ns, group.by = 'region',
               images = 'slice1.008um',image.alpha = 0.1,
               label = T, repel = T, label.size = 4,crop = F) + 
  scale_fill_manual(values = c(npg.cols[2]))


## RCTD
library(spacexr)

ref <- Clean_sct.inte.rm.lpt
Idents(ref)
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

### Age3
DefaultAssay(Age3.jz_db.hd.ns) <- "Spatial.008um"
Age3.jz_db.hd.ns <- FindVariableFeatures(Age3.jz_db.hd.ns)

Age3.jz_db.hd.ns <- ScaleData(Age3.jz_db.hd.ns) %>% 
  RunPCA() %>% 
  FindNeighbors( reduction = "pca", dims = 1:50) %>% 
  RunUMAP( reduction = "pca", reduction.name = "umap", return.model = T, dims = 1:30, verbose = T)
DimPlot(Age3.jz_db.hd.ns)
FeaturePlot(Age3.jz_db.hd.ns,features = 'S100a9',reduction= "umap",order = T, cols = c('gray90','red'))

counts_hd <- Age3.jz_db.hd.ns[["Spatial.008um"]]$counts
Age3.jz_db.hd.ns_cells_hd <- colnames(Age3.jz_db.hd.ns[["Spatial.008um"]])
coords <- GetTissueCoordinates(Age3.jz_db.hd.ns)[Age3.jz_db.hd.ns_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 80)
Age3.ns_RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

Age3.jz_db.hd.ns.back <- Age3.jz_db.hd.ns
Age3.jz_db.hd.ns <- Age3.jz_db.hd.ns.back
### add meta
Age3.jz_db.hd.ns <- AddMetaData(Age3.jz_db.hd.ns, metadata = Age3.ns_RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
Age3.jz_db.hd.ns$first_type <- as.character(Age3.jz_db.hd.ns$first_type)
# Age3.jz_db.hd.ns$first_type[is.na(Age3.jz_db.hd.ns$first_type)] <- "Unknown"
DefaultAssay(Age3.jz_db.hd.ns) <- "Spatial.008um"
table(Age3.jz_db.hd.ns$first_type)/ncol(Age3.jz_db.hd.ns)
table(Age3.jz_db.hd.ns$second_type)

Idents(Age3.jz_db.hd.ns) <- 'second_type'
cells <- CellsByIdentities(Age3.jz_db.hd.ns)

Age3.jz_db.hd.ns$first_type <- factor(Age3.jz_db.hd.ns$first_type, levels = c('GCs','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells'))
SpatialDimPlot(Age3.jz_db.hd.ns,group.by = 'first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()
SpatialDimPlot(Age3.jz_db.hd.ns, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(Age3.jz_db.hd.ns, cells.highlight = cells[c('Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)

### Young1
DefaultAssay(Young1.jz_db.hd.ns) <- "Spatial.008um"
Young1.jz_db.hd.ns <- FindVariableFeatures(Young1.jz_db.hd.ns)

Young1.jz_db.hd.ns <- ScaleData(Young1.jz_db.hd.ns) %>% 
  RunPCA() %>% 
  FindNeighbors( reduction = "pca", dims = 1:50) %>% 
  RunUMAP( reduction = "pca", reduction.name = "umap", return.model = T, dims = 1:30, verbose = T)
DimPlot(Young1.jz_db.hd.ns)
FeaturePlot(Young1.jz_db.hd.ns,features = 'S100a9',reduction= "umap",order = T, cols = c('gray90','red'))

counts_hd <- Young1.jz_db.hd.ns[["Spatial.008um"]]$counts
Young1.jz_db.hd.ns_cells_hd <- colnames(Young1.jz_db.hd.ns[["Spatial.008um"]])
coords <- GetTissueCoordinates(Young1.jz_db.hd.ns)[Young1.jz_db.hd.ns_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 80)
Young1.ns_RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
Young1.jz_db.hd.ns.back <- Young1.jz_db.hd.ns

### add meta
Young1.jz_db.hd.ns <- AddMetaData(Young1.jz_db.hd.ns, metadata = Young1.ns_RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
Young1.jz_db.hd.ns$first_type <- as.character(Young1.jz_db.hd.ns$first_type)
# Young1.jz_db.hd.ns$first_type[is.na(Young1.jz_db.hd.ns$first_type)] <- "Unknown"
DefaultAssay(Young1.jz_db.hd.ns) <- "Spatial.008um"
table(Young1.jz_db.hd.ns$first_type )/ncol(Young1.jz_db.hd.ns)
table(Young1.jz_db.hd.ns$second_type )


Idents(Young1.jz_db.hd.ns) <- 'second_type'
cells <- CellsByIdentities(Young1.jz_db.hd.ns)

# Young1.jz_db.hd.ns$first_type <- factor(Young1.jz_db.hd.ns$first_type, levels = c('GCs','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells'))
SpatialDimPlot(Young1.jz_db.hd.ns,group.by = 'first_type',label = F,repel = T,image.alpha = .1,crop = F) + scale_fill_igv()
SpatialDimPlot(Young1.jz_db.hd.ns, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(Young1.jz_db.hd.ns, cells.highlight = cells[c('Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)

SpatialFeaturePlot(Young1.jz_db.hd.ns, features = c('Ms4a1','Ptprc'),image.alpha = 0.3,)
SpatialFeaturePlot(Young1.jz_db.hd.ns, features = c('Prl8a9','Stra6'),image.alpha = 0.1)
### KO2
DefaultAssay(KO2.jz_db.hd.ns) <- "Spatial.008um"
KO2.jz_db.hd.ns <- FindVariableFeatures(KO2.jz_db.hd.ns)

KO2.jz_db.hd.ns <- ScaleData(KO2.jz_db.hd.ns) %>% 
  RunPCA() %>% 
  FindNeighbors( reduction = "pca", dims = 1:50) %>% 
  RunUMAP( reduction = "pca", reduction.name = "umap", return.model = T, dims = 1:30, verbose = T)
DimPlot(KO2.jz_db.hd.ns)
FeaturePlot(KO2.jz_db.hd.ns,features = 'S100a9',reduction= "umap",order = T, cols = c('gray90','red'))

counts_hd <- KO2.jz_db.hd.ns[["Spatial.008um"]]$counts
KO2.jz_db.hd.ns_cells_hd <- colnames(KO2.jz_db.hd.ns[["Spatial.008um"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd.ns)[KO2.jz_db.hd.ns_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 80)
KO2.ns_RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
KO2.jz_db.hd.ns.back <- KO2.jz_db.hd.ns

### add meta
KO2.jz_db.hd.ns <- AddMetaData(KO2.jz_db.hd.ns, metadata = KO2.ns_RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd.ns$first_type <- as.character(KO2.jz_db.hd.ns$first_type)
KO2.jz_db.hd.ns$first_type[is.na(KO2.jz_db.hd.ns$first_type)] <- "Unknown"

KO2.jz_db.hd.ns$second_type <- as.character(KO2.jz_db.hd.ns$second_type)
KO2.jz_db.hd.ns$second_type[is.na(KO2.jz_db.hd.ns$second_type)] <- "Unknown"

DefaultAssay(KO2.jz_db.hd.ns) <- "Spatial.008um"
table(KO2.jz_db.hd.ns$first_type )/ncol(KO2.jz_db.hd.ns)
table(KO2.jz_db.hd.ns$second_type )/ncol(KO2.jz_db.hd.ns)
Idents(KO2.jz_db.hd.ns) <- 'first_type'
cells <- CellsByIdentities(KO2.jz_db.hd.ns)

# KO2.jz_db.hd.ns$first_type <- factor(KO2.jz_db.hd.ns$first_type, levels = c('GCs','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells'))
SpatialDimPlot(KO2.jz_db.hd.ns,group.by = 'second_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()
SpatialDimPlot(KO2.jz_db.hd.ns,group.by = 'first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()

SpatialDimPlot(KO2.jz_db.hd.ns, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(KO2.jz_db.hd.ns, cells.highlight = cells[c('Macrophages')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)

SpatialFeaturePlot(KO2.jz_db.hd.ns, features = c('Ms4a1','Ptprc'),image.alpha = 0.3)
SpatialFeaturePlot(KO2.jz_db.hd.ns, features = c('Lepr','Stra6'),image.alpha = 0.1)

DimPlot(age.seu)

# 01.age devon on age seu-----
age3.db.8.barcode <- read.csv('../05.10x.hd/02.align/age3_8_db.csv')

## RCTD
Idents(age.seu) <- age.seu$cell.type.percise.new
ref <- age.seu
Idents(ref)
counts <- GetAssayData(age.seu,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

### Age3:FindVariableFeatures & SketchData & ScaleData had run in 07.split.hd.R

DefaultAssay(Age3.jz_db.hd) <- "sketch"
counts_hd <- Age3.jz_db.hd[["sketch"]]$counts
Age3.jz_db.hd_cells_hd <- colnames(Age3.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Age3.jz_db.hd)[Age3.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Age3_RCTD_age.seu <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
Age3.jz_db.hd.age.seu <- AddMetaData(Age3.jz_db.hd, metadata = Age3_RCTD_age.seu@results$results_df)

# project RCTD labels from sketched cells to all cells
Age3.jz_db.hd.age.seu$first_type <- as.character(Age3.jz_db.hd.age.seu$first_type)
Age3.jz_db.hd.age.seu$first_type[is.na(Age3.jz_db.hd.age.seu$first_type)] <- "Unknown"
Age3.jz_db.hd.age.seu <- ProjectData(
  object = Age3.jz_db.hd.age.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Age3.jz_db.hd.sketch",
  umap.model = "umap.Age3.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Age3.jz_db.hd.age.seu) <- "Spatial.008um"
Age3.jz_db.hd.age.seu[[]][,"full_first_type"]

Age3.jz_db.hd.age.seu$full_first_type <- factor(Age3.jz_db.hd.age.seu$full_first_type, levels = c('GCs','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells','Unknown'))
SpatialDimPlot(Age3.jz_db.hd.age.seu,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()
# Idents(Age3.jz_db.hd) <- "full_first_type"
# cells <- CellsByIdentities(Age3.jz_db.hd)
# excitatory_names <- c('GCs','Syns','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells')
SpatialDimPlot(Age3.jz_db.hd.age.seu, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
table(Age3.jz_db.hd.age.seu$full_first_type)/ncol(Age3.jz_db.hd.age.seu)

#### age use >5% | 1% cell type 
## RCTD
age.seu.5pct <- subset(age.seu, subset = (cell.type.percise.new == 'GCs' | 
                                            cell.type.percise.new == 'DSCs' |
                                            cell.type.percise.new == 'MSCs'|
                                          cell.type.percise.new == 'Macrophages' |
                                            cell.type.percise.new == 'Neutrophils' |
                                            cell.type.percise.new == 'DCs'))
Idents(age.seu.5pct) <- age.seu.5pct$cell.type.percise.new
ref <- age.seu.5pct
Idents(ref)
counts <- GetAssayData(age.seu.5pct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(Age3.jz_db.hd) <- "sketch"
counts_hd <- Age3.jz_db.hd[["sketch"]]$counts
Age3.jz_db.hd_cells_hd <- colnames(Age3.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Age3.jz_db.hd)[Age3.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Age3_RCTD_age.seu.5pct <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
Age3.jz_db.hd.age.seu.5pct <- AddMetaData(Age3.jz_db.hd, metadata = Age3_RCTD_age.seu.5pct@results$results_df)

# project RCTD labels from sketched cells to all cells
Age3.jz_db.hd.age.seu.5pct$first_type <- as.character(Age3.jz_db.hd.age.seu.5pct$first_type)
Age3.jz_db.hd.age.seu.5pct$first_type[is.na(Age3.jz_db.hd.age.seu.5pct$first_type)] <- "Unknown"
Age3.jz_db.hd.age.seu.5pct <- ProjectData(
  object = Age3.jz_db.hd.age.seu.5pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Age3.jz_db.hd.sketch",
  umap.model = "umap.Age3.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Age3.jz_db.hd.age.seu.5pct) <- "Spatial.008um"
Age3.jz_db.hd.age.seu.5pct[[]][,"full_first_type"]

SpatialDimPlot(Age3.jz_db.hd.age.seu.5pct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()
table(Age3.jz_db.hd.age.seu.5pct$full_first_type)/ncol(Age3.jz_db.hd.age.seu.5pct)

## 2pct
age.seu.2pct <- subset(age.seu, subset = (cell.type.percise.new == 'EryCs' |
                                            cell.type.percise.new == 'GMPs' |
                                            cell.type.percise.new == 'Mast cells'),invert = T)
Idents(age.seu.2pct) <- age.seu.2pct$cell.type.percise.new
ref <- age.seu.2pct
Idents(ref)
counts <- GetAssayData(age.seu.2pct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(Age3.jz_db.hd) <- "sketch"
counts_hd <- Age3.jz_db.hd[["sketch"]]$counts
Age3.jz_db.hd_cells_hd <- colnames(Age3.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Age3.jz_db.hd)[Age3.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Age3_RCTD_age.seu.2pct <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
Age3.jz_db.hd.age.seu.2pct <- AddMetaData(Age3.jz_db.hd, metadata = Age3_RCTD_age.seu.2pct@results$results_df)

# project RCTD labels from sketched cells to all cells
Age3.jz_db.hd.age.seu.2pct$first_type <- as.character(Age3.jz_db.hd.age.seu.2pct$first_type)
Age3.jz_db.hd.age.seu.2pct$first_type[is.na(Age3.jz_db.hd.age.seu.2pct$first_type)] <- "Unknown"
Age3.jz_db.hd.age.seu.2pct <- ProjectData(
  object = Age3.jz_db.hd.age.seu.2pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Age3.jz_db.hd.sketch",
  umap.model = "umap.Age3.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Age3.jz_db.hd.age.seu.2pct) <- "Spatial.008um"
Age3.jz_db.hd.age.seu.2pct[[]][,"full_first_type"]
Age3.jz_db.hd.age.seu.2pct$full_first_type <- factor(Age3.jz_db.hd.age.seu.2pct$full_first_type, levels = c('DSCs','GCs','ECs','Macrophages','MSCs','S-TGCs','SpTs','T cells','Unknown'))
SpatialDimPlot(Age3.jz_db.hd.age.seu.2pct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_nejm()
table(Age3.jz_db.hd.age.seu.2pct$full_first_type)/ncol(Age3.jz_db.hd.age.seu.2pct)

## 5 cell types
age.seu.5ct <- subset(age.seu, subset = (cell.type.percise.new == 'GCs' |
                                           cell.type.percise.new == 'DSCs' |
                                           cell.type.percise.new == 'S-TGCs'|
                                           cell.type.percise.new == 'SpTs' |
                                           cell.type.percise.new == 'ECs' ),invert = F)
Idents(age.seu.5ct) <- age.seu.5ct$cell.type.percise.new
ref <- age.seu.5ct
Idents(ref)
counts <- GetAssayData(age.seu.5ct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(Age3.jz_db.hd) <- "sketch"
counts_hd <- Age3.jz_db.hd[["sketch"]]$counts
Age3.jz_db.hd_cells_hd <- colnames(Age3.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Age3.jz_db.hd)[Age3.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Age3_RCTD_age.seu.5ct <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
Age3.jz_db.hd.age.seu.5ct <- AddMetaData(Age3.jz_db.hd, metadata = Age3_RCTD_age.seu.5ct@results$results_df)

# project RCTD labels from sketched cells to all cells
Age3.jz_db.hd.age.seu.5ct$first_type <- as.character(Age3.jz_db.hd.age.seu.5ct$first_type)
Age3.jz_db.hd.age.seu.5ct$first_type[is.na(Age3.jz_db.hd.age.seu.5ct$first_type)] <- "Unknown"
Age3.jz_db.hd.age.seu.5ct <- ProjectData(
  object = Age3.jz_db.hd.age.seu.5ct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Age3.jz_db.hd.sketch",
  umap.model = "umap.Age3.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Age3.jz_db.hd.age.seu.5ct) <- "Spatial.008um"
Age3.jz_db.hd.age.seu.5ct[[]][,"full_first_type"]
Age3.jz_db.hd.age.seu.5ct$full_first_type <- factor(Age3.jz_db.hd.age.seu.5ct$full_first_type, levels = c('DSCs','GCs','ECs','S-TGCs','SpTs','Unknown'))
SpatialDimPlot(Age3.jz_db.hd.age.seu.5ct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1,crop = F) + scale_fill_nejm()
table(Age3.jz_db.hd.age.seu.5ct$full_first_type)/ncol(Age3.jz_db.hd.age.seu.5ct)

Idents(Age3.jz_db.hd.age.seu.5ct) <- "full_first_type"
# cells <- CellsByIdentities(Age3.jz_db.hd.age.seu.5ct)
excitatory_names <- c('DSCs','GCs','ECs','S-TGCs','SpTs','Unknown')
SpatialDimPlot(Age3.jz_db.hd.age.seu.5ct, cells.highlight = cells[c('DSCs','GCs','ECs','S-TGCs','SpTs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3,crop = F)

# 02.young devon on young seu-----
save.image(compress = F)
Idents(young.seu) <- young.seu$cell.type.percise.new
ref <- young.seu
Idents(ref)
counts <- GetAssayData(young.seu,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(Young1.jz_db.hd) <- "sketch"

counts_hd <- Young1.jz_db.hd[["sketch"]]$counts
Young1.jz_db.hd_cells_hd <- colnames(Young1.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Young1.jz_db.hd)[Young1.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Young1_RCTD_young.seu <- run.RCTD(RCTD, doublet_mode = "doublet")
# Young1.jz_db.hd.back <- Young1.jz_db.hd
save.image(compress = F)

### add meta 
Young1.jz_db.hd.young.seu <- AddMetaData(Young1.jz_db.hd, metadata = Young1_RCTD_young.seu@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
Young1.jz_db.hd.young.seu$first_type <- as.character(Young1.jz_db.hd.young.seu$first_type)
Young1.jz_db.hd.young.seu$first_type[is.na(Young1.jz_db.hd.young.seu$first_type)] <- "Unknown"
Young1.jz_db.hd.young.seu <- ProjectData(
  object = Young1.jz_db.hd.young.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Young1.jz_db.hd.sketch",
  umap.model = "umap.Young1.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)

Young1.jz_db.hd.young.seu$second_type <- as.character(Young1.jz_db.hd.young.seu$second_type)
Young1.jz_db.hd.young.seu$second_type[is.na(Young1.jz_db.hd.young.seu$second_type)] <- "Unknown"
Young1.jz_db.hd.young.seu <- ProjectData(
  object = Young1.jz_db.hd.young.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Young1.jz_db.hd.sketch",
  umap.model = "umap.Young1.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_second_type = "second_type")
)

DefaultAssay(Young1.jz_db.hd.young.seu) <- "Spatial.008um"
Young1.jz_db.hd.young.seu[[]][,"full_first_type"]
table(Young1.jz_db.hd.young.seu$full_first_type)/ncol(Young1.jz_db.hd.young.seu)
table(Young1.jz_db.hd.young.seu$full_second_type)/ncol(Young1.jz_db.hd.young.seu)

SpatialDimPlot(Young1.jz_db.hd.young.seu, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(Young1.jz_db.hd.young.seu, cells.highlight = cells[c('Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 1,image.alpha = 0.3)

#### young use >5% cell type 
## RCTD
young.seu.5pct <- subset(young.seu, subset = (cell.type.percise.new == 'GCs' | 
                                            cell.type.percise.new == 'DSCs' |
                                            cell.type.percise.new == 'MSCs'|
                                            cell.type.percise.new == 'Macrophyoungs' |
                                            cell.type.percise.new == 'Neutrophils' |
                                            cell.type.percise.new == 'DCs'))
Idents(young.seu.5pct) <- young.seu.5pct$cell.type.percise.new
ref <- young.seu.5pct
Idents(ref)
counts <- GetAssayData(young.seu.5pct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(Young1.jz_db.hd) <- "sketch"

counts_hd <- Young1.jz_db.hd[["sketch"]]$counts
Young1.jz_db.hd_cells_hd <- colnames(Young1.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Young1.jz_db.hd)[Young1.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Young1_RCTD_young.seu.5pct <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
Young1.jz_db.hd.young.seu.5pct <- AddMetaData(Young1.jz_db.hd, metadata = Young1_RCTD_young.seu.5pct@results$results_df)

# project RCTD labels from sketched cells to all cells
Young1.jz_db.hd.young.seu.5pct$first_type <- as.character(Young1.jz_db.hd.young.seu.5pct$first_type)
Young1.jz_db.hd.young.seu.5pct$first_type[is.na(Young1.jz_db.hd.young.seu.5pct$first_type)] <- "Unknown"
Young1.jz_db.hd.young.seu.5pct <- ProjectData(
  object = Young1.jz_db.hd.young.seu.5pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Young1.jz_db.hd.sketch",
  umap.model = "umap.Young1.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Young1.jz_db.hd.young.seu.5pct) <- "Spatial.008um"
Young1.jz_db.hd.young.seu.5pct[[]][,"full_first_type"]
Young1.jz_db.hd.young.seu.5pct$full_first_type <- factor(Young1.jz_db.hd.young.seu.5pct$full_first_type, levels = c('DSCs','GCs','Macrophages','MSCs','DCs','Neutrophils','Unknown'))
SpatialDimPlot(Young1.jz_db.hd.young.seu.5pct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()
table(Young1.jz_db.hd.young.seu.5pct$full_first_type)/ncol(Young1.jz_db.hd.young.seu.5pct)

## 2pct
young.seu.2pct <- subset(young.seu, subset = (cell.type.percise.new == 'EryCs' |
                                            cell.type.percise.new == 'GMPs' |
                                            cell.type.percise.new == 'Mast cells'),invert = T)
Idents(young.seu.2pct) <- young.seu.2pct$cell.type.percise.new
ref <- young.seu.2pct
Idents(ref)
counts <- GetAssayData(young.seu.2pct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(Young1.jz_db.hd) <- "sketch"
counts_hd <- Young1.jz_db.hd[["sketch"]]$counts
Young1.jz_db.hd_cells_hd <- colnames(Young1.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Young1.jz_db.hd)[Young1.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Young1_RCTD_young.seu.2pct <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
Young1.jz_db.hd.young.seu.2pct <- AddMetaData(Young1.jz_db.hd, metadata = Young1_RCTD_young.seu.2pct@results$results_df)

# project RCTD labels from sketched cells to all cells
Young1.jz_db.hd.young.seu.2pct$first_type <- as.character(Young1.jz_db.hd.young.seu.2pct$first_type)
Young1.jz_db.hd.young.seu.2pct$first_type[is.na(Young1.jz_db.hd.young.seu.2pct$first_type)] <- "Unknown"
Young1.jz_db.hd.young.seu.2pct <- ProjectData(
  object = Young1.jz_db.hd.young.seu.2pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Young1.jz_db.hd.sketch",
  umap.model = "umap.Young1.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Young1.jz_db.hd.young.seu.2pct) <- "Spatial.008um"
Young1.jz_db.hd.young.seu.2pct[[]][,"full_first_type"]
Young1.jz_db.hd.young.seu.2pct$full_first_type <- factor(Young1.jz_db.hd.young.seu.2pct$full_first_type, levels = c('DSCs','GCs','ECs','Macrophages','S-TGCs','SpTs','Unknown'))
SpatialDimPlot(Young1.jz_db.hd.young.seu.2pct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_nejm()
table(Young1.jz_db.hd.young.seu.2pct$full_first_type)/ncol(Young1.jz_db.hd.young.seu.2pct)

## 5 cell types
young.seu.5ct <- subset(young.seu, subset = (cell.type.percise.new == 'GCs' |
                                           cell.type.percise.new == 'DSCs' |
                                           cell.type.percise.new == 'S-TGCs'|
                                           cell.type.percise.new == 'SpTs' |
                                           cell.type.percise.new == 'ECs' ),invert = F)
Idents(young.seu.5ct) <- young.seu.5ct$cell.type.percise.new
ref <- young.seu.5ct
Idents(ref)
counts <- GetAssayData(young.seu.5ct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(Young1.jz_db.hd) <- "sketch"
counts_hd <- Young1.jz_db.hd[["sketch"]]$counts
Young1.jz_db.hd_cells_hd <- colnames(Young1.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(Young1.jz_db.hd)[Young1.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
Young1_RCTD_young.seu.5ct <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
Young1.jz_db.hd.young.seu.5ct <- AddMetaData(Young1.jz_db.hd, metadata = Young1_RCTD_young.seu.5ct@results$results_df)

# project RCTD labels from sketched cells to all cells
Young1.jz_db.hd.young.seu.5ct$first_type <- as.character(Young1.jz_db.hd.young.seu.5ct$first_type)
Young1.jz_db.hd.young.seu.5ct$first_type[is.na(Young1.jz_db.hd.young.seu.5ct$first_type)] <- "Unknown"
Young1.jz_db.hd.young.seu.5ct <- ProjectData(
  object = Young1.jz_db.hd.young.seu.5ct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.Young1.jz_db.hd.sketch",
  umap.model = "umap.Young1.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(Young1.jz_db.hd.young.seu.5ct) <- "Spatial.008um"
Young1.jz_db.hd.young.seu.5ct[[]][,"full_first_type"]
Young1.jz_db.hd.young.seu.5ct$full_first_type <- factor(Young1.jz_db.hd.young.seu.5ct$full_first_type, levels = c('DSCs','GCs','ECs','S-TGCs','SpTs','Unknown'))
SpatialDimPlot(Young1.jz_db.hd.young.seu.5ct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()
table(Young1.jz_db.hd.young.seu.5ct$full_first_type)/ncol(Young1.jz_db.hd.young.seu.5ct)

# 03.ko devon on age or young seu-----
DefaultAssay(KO2.jz_db.hd) <- "sketch"

counts_hd <- KO2.jz_db.hd[["sketch"]]$counts
KO2.jz_db.hd_cells_hd <- colnames(KO2.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd)[KO2.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
KO2_RCTD_young.seu <- run.RCTD(RCTD, doublet_mode = "doublet")
# KO2.jz_db.hd.back <- KO2.jz_db.hd
save.image(compress = F)

### add meta 
KO2.jz_db.hd.young.seu <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD_young.seu@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd.young.seu$first_type <- as.character(KO2.jz_db.hd.young.seu$first_type)
KO2.jz_db.hd.young.seu$first_type[is.na(KO2.jz_db.hd.young.seu$first_type)] <- "Unknown"
KO2.jz_db.hd.young.seu <- ProjectData(
  object = KO2.jz_db.hd.young.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(KO2.jz_db.hd.young.seu) <- "Spatial.008um"
KO2.jz_db.hd.young.seu[[]][,"full_first_type"]
table(KO2.jz_db.hd.young.seu$full_first_type)/ncol(KO2.jz_db.hd.young.seu)

KO2.jz_db.hd.young.seu$second_type <- as.character(KO2.jz_db.hd.young.seu$second_type)
KO2.jz_db.hd.young.seu$second_type[is.na(KO2.jz_db.hd.young.seu$second_type)] <- "Unknown"
KO2.jz_db.hd.young.seu <- ProjectData(
  object = KO2.jz_db.hd.young.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_second_type = "second_type")
)
table(KO2.jz_db.hd.young.seu$full_second_type)/ncol(KO2.jz_db.hd.young.seu)

SpatialDimPlot(KO2.jz_db.hd, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.3)
SpatialDimPlot(KO2.jz_db.hd, cells.highlight = cells[c('Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 1,image.alpha = 0.3)


# age
ref <- age.seu
Idents(ref)
counts <- GetAssayData(age.seu,assay = 'RNA', layer = 'counts')
# counts <- ref[["RNA"]]
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(KO2.jz_db.hd) <- "sketch"

RCTD <- create.RCTD(query, reference, max_cores = 90)
KO2_RCTD_age.seu <- run.RCTD(RCTD, doublet_mode = "doublet")

### add meta 
KO2.jz_db.hd.age.seu <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD_age.seu@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd.age.seu$first_type <- as.character(KO2.jz_db.hd.age.seu$first_type)
KO2.jz_db.hd.age.seu$first_type[is.na(KO2.jz_db.hd.age.seu$first_type)] <- "Unknown"
KO2.jz_db.hd.age.seu <- ProjectData(
  object = KO2.jz_db.hd.age.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(KO2.jz_db.hd.age.seu) <- "Spatial.008um"
KO2.jz_db.hd.age.seu[[]][,"full_first_type"]
table(KO2.jz_db.hd.age.seu$full_first_type)/ncol(KO2.jz_db.hd.age.seu)

KO2.jz_db.hd.age.seu$second_type <- as.character(KO2.jz_db.hd.age.seu$second_type)
KO2.jz_db.hd.age.seu$second_type[is.na(KO2.jz_db.hd.age.seu$second_type)] <- "Unknown"
KO2.jz_db.hd.age.seu <- ProjectData(
  object = KO2.jz_db.hd.age.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_second_type = "second_type")
)
table(KO2.jz_db.hd.age.seu$full_second_type)/ncol(KO2.jz_db.hd.age.seu)

# young 5pct
DefaultAssay(KO2.jz_db.hd) <- "sketch"

counts_hd <- KO2.jz_db.hd[["sketch"]]$counts
KO2.jz_db.hd_cells_hd <- colnames(KO2.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd)[KO2.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
KO2_RCTD_young.seu.5pct <- run.RCTD(RCTD, doublet_mode = "doublet")

KO2.jz_db.hd.young.seu.5pct <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD_young.seu.5pct@results$results_df)

## young 2pct

DefaultAssay(KO2.jz_db.hd) <- "sketch"
counts_hd <- KO2.jz_db.hd[["sketch"]]$counts
KO2.jz_db.hd_cells_hd <- colnames(KO2.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd)[KO2.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
KO2_RCTD_young.seu.2pct <- run.RCTD(RCTD, doublet_mode = "doublet")

KO2.jz_db.hd.young.seu.2pct <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD_young.seu.2pct@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd.young.seu.2pct$first_type <- as.character(KO2.jz_db.hd.young.seu.2pct$first_type)
KO2.jz_db.hd.young.seu.2pct$first_type[is.na(KO2.jz_db.hd.young.seu.2pct$first_type)] <- "Unknown"
KO2.jz_db.hd.young.seu.2pct <- ProjectData(
  object = KO2.jz_db.hd.young.seu.2pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(KO2.jz_db.hd.young.seu.2pct) <- "Spatial.008um"
KO2.jz_db.hd.young.seu.2pct[[]][,"full_first_type"]
table(KO2.jz_db.hd.young.seu.2pct$full_first_type)/ncol(KO2.jz_db.hd.young.seu.2pct)

KO2.jz_db.hd.young.seu.2pct$second_type <- as.character(KO2.jz_db.hd.young.seu.2pct$second_type)
KO2.jz_db.hd.young.seu.2pct$second_type[is.na(KO2.jz_db.hd.young.seu.2pct$second_type)] <- "Unknown"
KO2.jz_db.hd.young.seu.2pct <- ProjectData(
  object = KO2.jz_db.hd.young.seu.2pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_second_type = "second_type")
)
table(KO2.jz_db.hd.young.seu.2pct$full_second_type)/ncol(KO2.jz_db.hd.age.seu)

KO2.jz_db.hd.young.seu.2pct$full_first_type <- factor(KO2.jz_db.hd.young.seu.2pct$full_first_type, levels = c('DSCs','GCs','MSCs','DCs','Neutrophils','Unknown'))
SpatialDimPlot(KO2.jz_db.hd.young.seu.2pct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_nejm()

# age 5|2 pct

Idents(age.seu.2pct) <- age.seu.2pct$cell.type.percise.new
ref <- age.seu.2pct
Idents(ref)
counts <- GetAssayData(age.seu.2pct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(KO2.jz_db.hd) <- "sketch"

counts_hd <- KO2.jz_db.hd[["sketch"]]$counts
KO2.jz_db.hd_cells_hd <- colnames(KO2.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd)[KO2.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
KO2_RCTD_age.seu.2pct <- run.RCTD(RCTD, doublet_mode = "doublet")
save.image(compress = F)

KO2.jz_db.hd.age.seu.2pct <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD_age.seu.2pct@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd.age.seu.2pct$first_type <- as.character(KO2.jz_db.hd.age.seu.2pct$first_type)
KO2.jz_db.hd.age.seu.2pct$first_type[is.na(KO2.jz_db.hd.age.seu.2pct$first_type)] <- "Unknown"
KO2.jz_db.hd.age.seu.2pct <- ProjectData(
  object = KO2.jz_db.hd.age.seu.2pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(KO2.jz_db.hd.age.seu.2pct) <- "Spatial.008um"
KO2.jz_db.hd.age.seu.2pct[[]][,"full_first_type"]
table(KO2.jz_db.hd.age.seu.2pct$full_first_type)/ncol(KO2.jz_db.hd.age.seu.2pct)

KO2.jz_db.hd.age.seu.2pct$full_first_type <- factor(KO2.jz_db.hd.age.seu.2pct$full_first_type, levels = c('DSCs','GCs','MSCs','DCs','Neutrophils','Unknown'))
SpatialDimPlot(KO2.jz_db.hd.age.seu.2pct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()

SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'Prl8a9')

Idents(age.seu.2pct) <- age.seu.2pct$cell.type.percise.new
ref <- age.seu.2pct
Idents(ref)
counts <- GetAssayData(age.seu.2pct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(KO2.jz_db.hd) <- "sketch"

counts_hd <- KO2.jz_db.hd[["sketch"]]$counts
KO2.jz_db.hd_cells_hd <- colnames(KO2.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd)[KO2.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
KO2_RCTD_age.seu.2pct <- run.RCTD(RCTD, doublet_mode = "doublet")
save.image(compress = F)

KO2.jz_db.hd.age.seu.2pct <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD_age.seu.2pct@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd.age.seu.2pct$first_type <- as.character(KO2.jz_db.hd.age.seu.2pct$first_type)
KO2.jz_db.hd.age.seu.2pct$first_type[is.na(KO2.jz_db.hd.age.seu.2pct$first_type)] <- "Unknown"
KO2.jz_db.hd.age.seu.2pct <- ProjectData(
  object = KO2.jz_db.hd.age.seu.2pct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(KO2.jz_db.hd.age.seu.2pct) <- "Spatial.008um"
KO2.jz_db.hd.age.seu.2pct[[]][,"full_first_type"]
table(KO2.jz_db.hd.age.seu.2pct$full_first_type)/ncol(KO2.jz_db.hd.age.seu.2pct)

KO2.jz_db.hd.age.seu.2pct$full_first_type <- factor(KO2.jz_db.hd.age.seu.2pct$full_first_type, levels = c('DSCs','GCs','MSCs','DCs','Neutrophils','Unknown'))
SpatialDimPlot(KO2.jz_db.hd.age.seu.2pct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1) + scale_fill_igv()

SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'Prl8a9')

## 5ct
Clean_sct.inte.rm.lpt.5ct <- subset(Clean_sct.inte.rm.lpt, subset = (cell.type.percise.new == 'GCs' |
                                                           cell.type.percise.new == 'DSCs' |
                                                           cell.type.percise.new == 'S-TGCs'|
                                                           cell.type.percise.new == 'SpTs' |
                                                           cell.type.percise.new == 'ECs' ),invert = F)

Idents(Clean_sct.inte.rm.lpt.5ct) <- Clean_sct.inte.rm.lpt.5ct$cell.type.percise.new
ref <- Clean_sct.inte.rm.lpt.5ct
Idents(ref)
counts <- GetAssayData(Clean_sct.inte.rm.lpt.5ct,assay = 'RNA', layer = 'counts')
cluster <- as.factor(ref$cell.type.percise.new)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI,n_max_cells = 5000)

DefaultAssay(KO2.jz_db.hd) <- "sketch"

counts_hd <- KO2.jz_db.hd[["sketch"]]$counts
KO2.jz_db.hd_cells_hd <- colnames(KO2.jz_db.hd[["sketch"]])
coords <- GetTissueCoordinates(KO2.jz_db.hd)[KO2.jz_db.hd_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 90)
KO2_RCTD_all.seu.5ct <- run.RCTD(RCTD, doublet_mode = "doublet")
save.image(compress = F)

KO2.jz_db.hd.all.seu.5ct <- AddMetaData(KO2.jz_db.hd, metadata = KO2_RCTD_all.seu.5ct@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
KO2.jz_db.hd.all.seu.5ct$first_type <- as.character(KO2.jz_db.hd.all.seu.5ct$first_type)
KO2.jz_db.hd.all.seu.5ct$first_type[is.na(KO2.jz_db.hd.all.seu.5ct$first_type)] <- "Unknown"
KO2.jz_db.hd.all.seu.5ct <- ProjectData(
  object = KO2.jz_db.hd.all.seu.5ct,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.KO2.jz_db.hd.sketch",
  umap.model = "umap.KO2.jz_db.hd.sketch",
  dims = 1:30,
  refdata = list(full_first_type = "first_type")
)
DefaultAssay(KO2.jz_db.hd.all.seu.5ct) <- "Spatial.008um"
KO2.jz_db.hd.all.seu.5ct[[]][,"full_first_type"]
table(KO2.jz_db.hd.all.seu.5ct$full_first_type)/ncol(KO2.jz_db.hd.all.seu.5ct)

KO2.jz_db.hd.all.seu.5ct$full_first_type <- factor(KO2.jz_db.hd.all.seu.5ct$full_first_type, levels = c('DSCs','GCs','ECs','S-TGCs','SpTs','Unknown'))
SpatialDimPlot(KO2.jz_db.hd.all.seu.5ct,group.by = 'full_first_type',label = F,repel = T,image.alpha = .1,crop = F) + scale_fill_igv()


# 04.addmodulescore ----
table(rm.lpt.ct.deg$cluster)

## macrop
macrop.gene <- list(c(rm.lpt.ct.deg[rm.lpt.ct.deg$cluster == 'Macrophages',]$gene))
macrop.gene <- list(c(rm.lpt.ct.deg.top10[rm.lpt.ct.deg.top10$cluster == 'Macrophages',]$gene))
summary(Age3.jz_db.hd.age.seu$macrop.score1)
summary(Young1.jz_db.hd.young.seu$macrop.score1)
summary(KO2.jz_db.hd.age.seu$macrop.score1)

VlnPlot(KO2.jz_db.hd.age.seu, features = 'macrop.score1',pt.size = 0, group.by = )

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = macrop.gene, name = 'macrop.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'macrop.score1',crop = F,image.alpha = 0.2, 
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$macrop.score1, probs = 0.999))
Age3.jz_db.hd.age.seu$macrop.like <- ifelse(
  Age3.jz_db.hd.age.seu$macrop.score1 > quantile(Age3.jz_db.hd.age.seu$macrop.score1, probs = 0.9),'macrop.like','other')

Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = macrop.gene, name = 'macrop.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'macrop.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$macrop.score1, probs = 0.999))

DefaultAssay(KO2.jz_db.hd.age.seu) <- 'Spatial.008um'
KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = macrop.gene, name = 'macrop.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'macrop.score1',crop = F,image.alpha = 0.2)

KO2.jz_db.hd.young.seu <- AddModuleScore(KO2.jz_db.hd.young.seu, features = macrop.gene, name = 'macrop.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'macrop.score1',crop = T,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.young.seu$macrop.score1, probs = 0.99))

## macrop 5ct
Age3.jz_db.hd.age.seu.5ct <- AddModuleScore(Age3.jz_db.hd.age.seu.5ct, features = macrop.gene, name = 'macrop.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu.5ct,features = 'macrop.score1',crop = F,image.alpha = 0.2, 
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu.5ct$macrop.score1, probs = 0.999))
Age3.jz_db.hd.age.seu.5ct$macrop.like <- ifelse(
  Age3.jz_db.hd.age.seu.5ct$macrop.score1 > quantile(Age3.jz_db.hd.age.seu.5ct$macrop.score1, probs = 0.9),'macrop.like','other')

Young1.jz_db.hd.young.seu.5ct <- AddModuleScore(Young1.jz_db.hd.young.seu.5ct, features = macrop.gene, name = 'macrop.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu.5ct,features = 'macrop.score1',crop = F,image.alpha = 0.2, 
                   max.cutoff = quantile(Young1.jz_db.hd.young.seu.5ct$macrop.score1, probs = 0.999))
Young1.jz_db.hd.young.seu.5ct$macrop.like <- ifelse(
  Young1.jz_db.hd.young.seu.5ct$macrop.score1 > quantile(Young1.jz_db.hd.young.seu.5ct$macrop.score1, probs = 0.9),'macrop.like','other')

KO2.jz_db.hd.all.seu.5ct <- AddModuleScore(KO2.jz_db.hd.all.seu.5ct, features = macrop.gene, name = 'macrop.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.all.seu.5ct,features = 'macrop.score1',crop = F,image.alpha = 0.2, 
                   max.cutoff = quantile(KO2.jz_db.hd.all.seu.5ct$macrop.score1, probs = 0.999))
KO2.jz_db.hd.all.seu.5ct$macrop.like <- ifelse(
  KO2.jz_db.hd.all.seu.5ct$macrop.score1 > quantile(KO2.jz_db.hd.all.seu.5ct$macrop.score1, probs = 0.9),'macrop.like','other')

# filter na
### age
Age3.jz_db.hd.age.seu.5ct.filter.na.meta <- Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct@meta.data$full_first_type != 'Unknown',]
Age3.jz_db.hd.age.seu.5ct.filter.na.meta$Barcode <- rownames(Age3.jz_db.hd.age.seu.5ct.filter.na.meta)
Age3.jz_db.hd.age.seu.5ct.filter.na.meta <- Age3.jz_db.hd.age.seu.5ct.filter.na.meta[,c(2,3,10:14,25,27:29)]

Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop <- Age3.jz_db.hd.age.seu.5ct@meta.data
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$guess.ct <- NA
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$guess.ct <- as.character(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$full_first_type)
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$guess.ct[which(str_detect(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$macrop.like, "^macrop.like"))] <- "Macrophages.like"
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop[Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$guess.ct != 'Unknown',]
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$Barcode <- rownames(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop)
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop[,c(2,3,10:14,25,27:30)]

Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon <- Age3.jz_db.hd.age.seu.5ct@meta.data
Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon <- Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon[Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon$macrop.like != 'other' & Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon$full_first_type != 'Unknown',]
Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon$Barcode <- rownames(Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon)
Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon <- Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon[,c(2,3,10:14,25,27:29)]

write.csv(Age3.jz_db.hd.age.seu.5ct.filter.na.meta, file = '../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.5ct.filter.na.meta.csv',quote = F,col.names = T)
write.csv(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop, file = '../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.csv',quote = F,col.names = T)
write.csv(Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon, file = '../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon.csv',quote = F,col.names = T)
table(Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon$full_first_type)/nrow(Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon)
table(Age3.jz_db.hd.age.seu.5ct.filter.na.only.macrop.decon$full_first_type)/table(Age3.jz_db.hd.age.seu.5ct.filter.na.meta$full_first_type)

ggplot(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop, aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=guess.ct)) + 
  geom_point(size=.1) 

### young
Young1.jz_db.hd.young.seu.5ct.filter.na.meta <- Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct@meta.data$full_first_type != 'Unknown',]
Young1.jz_db.hd.young.seu.5ct.filter.na.meta$Barcode <- rownames(Young1.jz_db.hd.young.seu.5ct.filter.na.meta)
Young1.jz_db.hd.young.seu.5ct.filter.na.meta <- Young1.jz_db.hd.young.seu.5ct.filter.na.meta[,c(2,3,10:14,30,32:34)]

Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop <- Young1.jz_db.hd.young.seu.5ct@meta.data
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$guess.ct <- NA
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$guess.ct <- as.character(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$full_first_type)
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$guess.ct[which(str_detect(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$macrop.like, "^macrop.like"))] <- "Macrophages.like"
Young1.jz_db.hd.young.seu.5ct$guess.ct <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$guess.ct
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop[Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$guess.ct != 'Unknown',]
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$Barcode <- rownames(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop)
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop[,c(2,3,10:14,30,32:35)]

Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon <- Young1.jz_db.hd.young.seu.5ct@meta.data
Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon <- Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon[Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon$macrop.like != 'other' & Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon$full_first_type != 'Unknown',]
Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon$Barcode <- rownames(Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon)
Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon <- Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon[,c(2,3,10:14,30,32:34)]

write.csv(Young1.jz_db.hd.young.seu.5ct.filter.na.meta, file = '../05.10x.hd/03.py.analysis/Young1.jz_db.hd.young.seu.5ct.filter.na.meta.csv',quote = F,col.names = T)
write.csv(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop, file = '../05.10x.hd/03.py.analysis/Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.csv',quote = F,col.names = T)
write.csv(Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon, file = '../05.10x.hd/03.py.analysis/Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon.csv',quote = F,col.names = T)
table(Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon$full_first_type)/nrow(Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon)
table(Young1.jz_db.hd.young.seu.5ct.filter.na.only.macrop.decon$full_first_type)/table(Young1.jz_db.hd.young.seu.5ct.filter.na.meta$full_first_type)

### ko
KO2.jz_db.hd.all.seu.5ct.filter.na.meta <- KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct@meta.data$full_first_type != 'Unknown',]
KO2.jz_db.hd.all.seu.5ct.filter.na.meta$Barcode <- rownames(KO2.jz_db.hd.all.seu.5ct.filter.na.meta)
KO2.jz_db.hd.all.seu.5ct.filter.na.meta <- KO2.jz_db.hd.all.seu.5ct.filter.na.meta[,c(2,3,10:14,29,31:33)]

KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop <- KO2.jz_db.hd.all.seu.5ct@meta.data
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct <- NA
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct <- as.character(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$full_first_type)
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct[which(str_detect(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$macrop.like, "^macrop.like"))] <- "Macrophages.like"
table(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct)

KO2.jz_db.hd.all.seu.5ct$guess.ct <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop[KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct != 'Unknown',]
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$Barcode <- rownames(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop)
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop[,c(2,3,10:14,29,31:34)]

KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon <- KO2.jz_db.hd.all.seu.5ct@meta.data
KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon <- KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon[KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon$macrop.like != 'other' & KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon$full_first_type != 'Unknown',]
KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon$Barcode <- rownames(KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon)
KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon <- KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon[,c(2,3,10:14,29,31:33)]

write.csv(KO2.jz_db.hd.all.seu.5ct.filter.na.meta, file = '../05.10x.hd/03.py.analysis/KO2.jz_db.hd.all.seu.5ct.filter.na.meta.csv',quote = F,col.names = T)
write.csv(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop, file = '../05.10x.hd/03.py.analysis/KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.csv',quote = F,col.names = T)
write.csv(KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon, file = '../05.10x.hd/03.py.analysis/KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon.csv',quote = F,col.names = T)
table(KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon$full_first_type)/nrow(KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon)
table(KO2.jz_db.hd.all.seu.5ct.filter.na.only.macrop.decon$full_first_type)/table(KO2.jz_db.hd.all.seu.5ct.filter.na.meta$full_first_type)

Idents(Age3.jz_db.hd.age.seu.5ct) <- "full_first_type"
cells <- CellsByIdentities(Age3.jz_db.hd.age.seu.5ct)
SpatialDimPlot(Age3.jz_db.hd.age.seu.5ct, cells.highlight = cells[c('ECs','SpTs','S-TGCs','GCs','DSCs','Unknown')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 4,image.alpha = 0.2,crop = F)

Age3.jz_db.hd.age.seu.5ct$guess.ct <- NA
Age3.jz_db.hd.age.seu.5ct$guess.ct <- as.character(Age3.jz_db.hd.age.seu.5ct$full_first_type)
Age3.jz_db.hd.age.seu.5ct$guess.ct[which(str_detect(Age3.jz_db.hd.age.seu.5ct$macrop.like, "^macrop.like"))] <- "Macrophages.like"
Idents(Age3.jz_db.hd.age.seu.5ct) <- "guess.ct"
cells <- CellsByIdentities(Age3.jz_db.hd.age.seu.5ct)
SpatialDimPlot(Age3.jz_db.hd.age.seu.5ct, cells.highlight = cells[c('ECs')], cols.highlight = c("red3", "grey80"), facet.highlight = T, combine = T, ncol = 1,image.alpha = 0.2,crop = F,label.size = 1)

# get neibor age -----

# strict ED: √((Δx)^2+(Δy)^2)≤2
## AGE 
age3.tp <- arrow::read_parquet(
  "../05.10x.hd/02.align/A3/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet"
)
rownames(age3.tp) <- age3.tp$barcode
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop %>% 
  left_join(y = age3.tp[,c("barcode","array_row","array_col")], by = c('Barcode' = 'barcode'))

coords <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]

mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

### rough
library(RANN)
res <- nn2(data  = as.matrix(dsc[,c("array_row","array_col")]),
           query = as.matrix(mac[,c("array_row","array_col")]),
           searchtype="radius", radius=2)

mac$neighbors_DSCs <- lapply(seq_len(nrow(mac)), function(i) {
  idxs <- res$nn.idx[i, ]         # 第 i 行：这个 Mac 拓展半径内所有候选 DSC 索引
  idxs <- idxs[idxs > 0]          # 去掉填充用的 0
  dsc$Barcode[idxs]               # 映射回 DSC 的 barcode
})

mac <- mac %>%
  mutate(
    neighbor_DSCs_str = sapply(neighbors_DSCs, function(x) {
      if (length(x)==0) return(NA_character_)
      paste(x, collapse = ";")
    })
  )
unique(unlist(mac$neighbors_DSCs))

### percise
radius <- 2
eps    <- 1e-6
res <- nn2(data  = as.matrix(dsc[,c("array_row","array_col")]),
           query = as.matrix(mac[,c("array_row","array_col")]),
           searchtype="radius", radius=radius + eps)

mac$neighbors_euc <- lapply(seq_len(nrow(mac)), function(i) {
  idxs <- res$nn.idx[i, ]
  idxs <- idxs[idxs > 0]
  # 进一步手动计算欧氏距离并严格筛一次
  keep <- which(
    sqrt((dsc$array_row[idxs] - mac$array_row[i])^2 +
           (dsc$array_col[idxs] - mac$array_col[i])^2)
    <= radius
  )
  dsc$Barcode[idxs[keep]]
})

Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor) <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_euc)),]$cell_type <- "Neibor.DSCs"

ggplot(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()
ggplot(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=-array_col, y=-array_row, color=cell_type)) + 
  geom_point(size=.02) + scale_color_npg() + theme_bw()

table(Age3.jz_db.hd.age.seu.5ct$guess.ct)
Age3.jz_db.hd.age.seu.5ct$cell_type_euc <- Age3.jz_db.hd.age.seu.5ct$guess.ct
Age3.jz_db.hd.age.seu.5ct@meta.data[unique(unlist(mac$neighbors_euc)),]$cell_type_euc <- "Near_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_euc)

## YOUNG

Young1.tp <- arrow::read_parquet(
  "../05.10x.hd/02.align/A3/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet"
)
rownames(Young1.tp) <- Young1.tp$barcode
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop %>% 
  left_join(y = Young1.tp[,c("barcode","array_row","array_col")], by = c('Barcode' = 'barcode'))

coords <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]

mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius <- 2
eps    <- 1e-6
res <- nn2(data  = as.matrix(dsc[,c("array_row","array_col")]),
           query = as.matrix(mac[,c("array_row","array_col")]),
           searchtype="radius", radius=radius + eps)

mac$neighbors_euc <- lapply(seq_len(nrow(mac)), function(i) {
  idxs <- res$nn.idx[i, ]
  idxs <- idxs[idxs > 0]
  # 进一步手动计算欧氏距离并严格筛一次
  keep <- which(
    sqrt((dsc$array_row[idxs] - mac$array_row[i])^2 +
           (dsc$array_col[idxs] - mac$array_col[i])^2)
    <= radius
  )
  dsc$Barcode[idxs[keep]]
})

Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor) <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_euc)),]$cell_type <- "Neibor.DSCs"

ggplot(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

table(Young1.jz_db.hd.young.seu.5ct$guess.ct)
table(Young1.jz_db.hd.young.seu.5ct$full_first_type)
Young1.jz_db.hd.young.seu.5ct$cell_type_euc <- Young1.jz_db.hd.young.seu.5ct$guess.ct
Young1.jz_db.hd.young.seu.5ct@meta.data[unique(unlist(mac$neighbors_euc)),]$cell_type_euc <- "Near_mac_DSCs"
table(Young1.jz_db.hd.young.seu.5ct$cell_type_euc)

## KO

KO2.tp <- arrow::read_parquet(
  "../05.10x.hd/02.align/A4/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet"
)
rownames(KO2.tp) <- KO2.tp$barcode
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop %>% 
  left_join(y = KO2.tp[,c("barcode","array_row","array_col")], by = c('Barcode' = 'barcode'))

coords <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
table(coords$guess.ct)
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius <- 2
eps    <- 1e-6
res <- nn2(data  = as.matrix(dsc[,c("array_row","array_col")]),
           query = as.matrix(mac[,c("array_row","array_col")]),
           searchtype="radius", radius=radius + eps)

mac$neighbors_euc <- lapply(seq_len(nrow(mac)), function(i) {
  idxs <- res$nn.idx[i, ]
  idxs <- idxs[idxs > 0]
  # 进一步手动计算欧氏距离并严格筛一次
  keep <- which(
    sqrt((dsc$array_row[idxs] - mac$array_row[i])^2 +
           (dsc$array_col[idxs] - mac$array_col[i])^2)
    <= radius
  )
  dsc$Barcode[idxs[keep]]
})

KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor) <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_euc)),]$cell_type <- "Neibor.DSCs"
table(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type)
ggplot(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()
ggplot(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=-array_col, y=-array_row, color=cell_type)) + 
  geom_point(size=.02) + scale_color_npg() + theme_bw()

# KO2.jz_db.hd.all.seu.5ct$guess.ct <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type
# table(KO2.jz_db.hd.all.seu.5ct$full_first_type)
KO2.jz_db.hd.all.seu.5ct$cell_type_euc <- KO2.jz_db.hd.all.seu.5ct$guess.ct
KO2.jz_db.hd.all.seu.5ct@meta.data[unique(unlist(mac$neighbors_euc)),]$cell_type_euc <- "Near_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_euc)
table(KO2.jz_db.hd.all.seu.5ct$guess.ct)

# Chebyshev: |Δx|≤2 & |Δy|≤2

## age
coords <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_bin <- 1
mac$neighbors_sq <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) <= radius_bin,
           abs(array_row - yi) <= radius_bin) %>%
    pull(Barcode)
})

Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor) <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_sq)),]$cell_type <- "Neibor.DSCs"

ggplot(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

table(Age3.jz_db.hd.age.seu.5ct$guess.ct)
Age3.jz_db.hd.age.seu.5ct$cell_type_sq <- Age3.jz_db.hd.age.seu.5ct$guess.ct
Age3.jz_db.hd.age.seu.5ct@meta.data[unique(unlist(mac$neighbors_sq)),]$cell_type_sq <- "Near_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_sq)

## young

coords <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_bin <- 1
mac$neighbors_sq <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) <= radius_bin,
           abs(array_row - yi) <= radius_bin) %>%
    pull(Barcode)
})

Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor) <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_sq)),]$cell_type <- "Neibor.DSCs"

ggplot(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

table(Young1.jz_db.hd.young.seu.5ct$guess.ct)
table(Young1.jz_db.hd.young.seu.5ct$full_first_type)
Young1.jz_db.hd.young.seu.5ct$cell_type_sq <- Young1.jz_db.hd.young.seu.5ct$guess.ct
Young1.jz_db.hd.young.seu.5ct@meta.data[unique(unlist(mac$neighbors_sq)),]$cell_type_sq <- "Near_mac_DSCs"
table(Young1.jz_db.hd.young.seu.5ct$cell_type_sq)

## KO
coords <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_bin <- 1
mac$neighbors_sq <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) <= radius_bin,
           abs(array_row - yi) <= radius_bin) %>%
    pull(Barcode)
})

KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor) <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_sq)),]$cell_type <- "Neibor.DSCs"
table(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type)
ggplot(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

KO2.jz_db.hd.all.seu.5ct$cell_type_sq <- KO2.jz_db.hd.all.seu.5ct$guess.ct
KO2.jz_db.hd.all.seu.5ct@meta.data[unique(unlist(mac$neighbors_sq)),]$cell_type_sq <- "Near_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_sq)
table(KO2.jz_db.hd.all.seu.5ct$guess.ct)
# Manhattan: |Δx|+|Δy|≤2

## age
coords <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_m  <- 2
mac$neighbors_mh <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) + abs(array_row - yi) <= radius_m) %>%
    pull(Barcode)
})

Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor) <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_mh)),]$cell_type <- "Neibor.DSCs"

ggplot(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

table(Age3.jz_db.hd.age.seu.5ct$guess.ct)
Age3.jz_db.hd.age.seu.5ct$cell_type_mh <- Age3.jz_db.hd.age.seu.5ct$guess.ct
Age3.jz_db.hd.age.seu.5ct@meta.data[unique(unlist(mac$neighbors_mh)),]$cell_type_mh <- "Near_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_mh)

## young

coords <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_m  <- 2
mac$neighbors_mh <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) + abs(array_row - yi) <= radius_m) %>%
    pull(Barcode)
})

Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor) <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_mh)),]$cell_type <- "Neibor.DSCs"

ggplot(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

table(Young1.jz_db.hd.young.seu.5ct$guess.ct)
table(Young1.jz_db.hd.young.seu.5ct$full_first_type)
Young1.jz_db.hd.young.seu.5ct$cell_type_mh <- Young1.jz_db.hd.young.seu.5ct$guess.ct
Young1.jz_db.hd.young.seu.5ct@meta.data[unique(unlist(mac$neighbors_mh)),]$cell_type_mh <- "Near_mac_DSCs"
table(Young1.jz_db.hd.young.seu.5ct$cell_type_mh)

## KO
coords <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_m  <- 2
mac$neighbors_mh <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) + abs(array_row - yi) <= radius_m) %>%
    pull(Barcode)
})

KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor) <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_mh)),]$cell_type <- "Neibor.DSCs"
table(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type)
ggplot(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

KO2.jz_db.hd.all.seu.5ct$cell_type_mh <- KO2.jz_db.hd.all.seu.5ct$guess.ct
KO2.jz_db.hd.all.seu.5ct@meta.data[unique(unlist(mac$neighbors_mh)),]$cell_type_mh <- "Near_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_mh)
table(KO2.jz_db.hd.all.seu.5ct$guess.ct)

# DSC degs in hd ----
Idents(dsc.8um) <- dsc.8um$sample
dsc.8um.age_young_degs <- FindMarkers(dsc.8um, ident.1 = 'Age',ident.2 = 'Young',logfc.threshold = 0.2) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = 'gene')
dsc.8um.age_young_degs$type <- factor(ifelse(dsc.8um.age_young_degs$avg_log2FC > 0,'Age','Young'))
table(dsc.8um.age_young_degs$type )

dsc.8um.age_ko_degs <- FindMarkers(dsc.8um, ident.1 = 'Age',ident.2 = 'KO',logfc.threshold = 0.2) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = 'gene')
dsc.8um.age_ko_degs$type <- factor(ifelse(dsc.8um.age_ko_degs$avg_log2FC > 0,'Age','KO'))
table(dsc.8um.age_ko_degs$type )
dsc.hd.8um.deg <- list(dsc.8um.age_young_degs = dsc.8um.age_young_degs,
                       dsc.8um.age_ko_degs = dsc.8um.age_ko_degs)
rio::export(dsc.hd.8um.deg, file = 'dsc.hd.8um.deg.xlsx')

## DSC RCTD enriched addmodulescoe? ----
dsc.gene <- list(rm.lpt.ct.deg[rm.lpt.ct.deg$cluster == 'DSCs',]$gene)
rm.lpt.ct.deg.top11 <- rm.lpt.ct.deg %>% dplyr::group_by(cluster) %>% top_n(n = 11,wt = avg_log2FC)
dsc.gene <- list(rm.lpt.ct.deg.top11[rm.lpt.ct.deg.top11$cluster == 'DSCs',]$gene)

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = dsc.gene, name = 'dsc.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'dsc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$dsc.score1, probs = 0.999))
Age3.jz_db.hd.age.seu$dsc.like <- ifelse(
  Age3.jz_db.hd.age.seu$dsc.score1 > quantile(Age3.jz_db.hd.age.seu$dsc.score1, probs = 0.8),'dsc.like','other')

ams <- rownames(Age3.jz_db.hd.age.seu@meta.data[Age3.jz_db.hd.age.seu$dsc.like == 'dsc.like',])
rctd <- rownames( Age3.jz_db.hd.age.seu@meta.data[Age3.jz_db.hd.age.seu$full_first_type == 'DSCs',])
all.char <- union(ams,rctd)

N <- length(all.char)
n <- length(intersect(ams,rctd))
k <- length(ams)
K <- length(rctd)


phyper(n-1, K, N-k ,n)
phyper(k-1, n, N-n ,m, lower.tail = T)

overlap_table <- matrix(c(k, n - k, m - k, N - n - m + k), nrow = 2, 
                        dimnames = list(c("In Group2", "Not in Group2"), c("In Group1", "Not in Group1")))
fisher_result <- fisher.test(overlap_table)
cat("p-value from Fisher's test:", fisher_result$p.value, "\n")


bin.mat <- data.frame(
  ams = as.numeric(all.char %in% ams),
  rctd = as.numeric(all.char %in% rctd)
)
bin.mat$match <- ifelse(bin.mat$ams == bin.mat$rctd, 'yes', 'no')
table(bin.mat$match)
rownames(bin.mat) <- all.char
cor(bin.mat$ams, bin.mat$rctd,method = 'spearman')
cor(ams,rctd)

jacd <- sum(bin.mat$ams & bin.mat$rctd)/sum(bin.mat$ams | bin.mat$rctd)



Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = dsc.gene, name = 'dsc.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'dsc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$dsc.score1, probs = 0.999))

Young1.jz_db.hd.young.seu$dsc.like <- ifelse(
  Young1.jz_db.hd.young.seu$dsc.score1 > quantile(Young1.jz_db.hd.young.seu$dsc.score1, probs = 0.85),'dsc.like','other')
table(Young1.jz_db.hd.young.seu$dsc.like)

ams <- rownames(Young1.jz_db.hd.young.seu@meta.data[Young1.jz_db.hd.young.seu$dsc.like == 'dsc.like',])
rctd <- rownames( Young1.jz_db.hd.young.seu@meta.data[Young1.jz_db.hd.young.seu$full_first_type == 'DSCs',])
all.char <- union(ams,rctd)

N <- length(all.char)
k <- length(intersect(ams,rctd))
n <- length(ams)
m <- length(rctd)
phyper(k-1, m, N-m ,n, lower.tail = T)
phyper(k-1, n, N-n ,m, lower.tail = T)

overlap_table <- matrix(c(k, n - k, m - k, N - n - m + k), nrow = 2, 
                        dimnames = list(c("In Group2", "Not in Group2"), c("In Group1", "Not in Group1")))
fisher_result <- fisher.test(overlap_table)
cat("p-value from Fisher's test:", fisher_result$p.value, "\n")

bin.mat <- data.frame(
  ams = as.numeric(all.char %in% ams),
  rctd = as.numeric(all.char %in% rctd)
)
bin.mat$match <- ifelse(bin.mat$ams == bin.mat$rctd, 'yes', 'no')
table(bin.mat$match)
rownames(bin.mat) <- all.char
cor(bin.mat$ams, bin.mat$rctd,method = 'spearman')
# cor(ams,rctd)
jacd <- sum(bin.mat$ams & bin.mat$rctd)/sum(bin.mat$ams | bin.mat$rctd)
jacd


DefaultAssay(KO2.jz_db.hd.age.seu) <- 'Spatial.008um'
KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = dsc.gene, name = 'dsc.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'dsc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.age.seu$dsc.score1, probs = 0.999))

KO2.jz_db.hd.age.seu$dsc.like <- ifelse(
  KO2.jz_db.hd.age.seu$dsc.score1 > quantile(KO2.jz_db.hd.age.seu$dsc.score1, probs = 0.98),'dsc.like','other')
table(KO2.jz_db.hd.age.seu$dsc.like)

ams <- rownames(KO2.jz_db.hd.age.seu@meta.data[KO2.jz_db.hd.age.seu$dsc.like == 'dsc.like',])
rctd <- rownames( KO2.jz_db.hd.age.seu@meta.data[KO2.jz_db.hd.age.seu$full_first_type == 'DSCs',])
all.char <- union(ams,rctd)

N <- length(all.char)
k <- length(intersect(ams,rctd))
n <- length(ams)
m <- length(rctd)
phyper(k-1, m, N-m ,n, lower.tail = T)
phyper(k-1, n, N-n ,m, lower.tail = T)

overlap_table <- matrix(c(k, n - k, m - k, N - n - m + k), nrow = 2, 
                        dimnames = list(c("In Group2", "Not in Group2"), c("In Group1", "Not in Group1")))
fisher_result <- fisher.test(overlap_table)
cat("p-value from Fisher's test:", fisher_result$p.value, "\n")

bin.mat <- data.frame(
  ams = as.numeric(all.char %in% ams),
  rctd = as.numeric(all.char %in% rctd)
)
bin.mat$match <- ifelse(bin.mat$ams == bin.mat$rctd, 'yes', 'no')
table(bin.mat$match)
rownames(bin.mat) <- all.char
cor(bin.mat$ams, bin.mat$rctd,method = 'spearman')
# cor(ams,rctd)
jacd <- sum(bin.mat$ams & bin.mat$rctd)/sum(bin.mat$ams | bin.mat$rctd)
jacd


KO2.jz_db.hd.young.seu <- AddModuleScore(KO2.jz_db.hd.young.seu, features = dsc.gene, name = 'dsc.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'dsc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.young.seu$dsc.score1, probs = 0.99))

## GCs
gc.gene <- list(rm.lpt.ct.deg[rm.lpt.ct.deg$cluster == 'GCs',]$gene)
gc.gene <- list(rm.lpt.ct.deg.top10[rm.lpt.ct.deg.top10$cluster == 'GCs',]$gene)

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = gc.gene, name = 'gc.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'gc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$gc.score1, probs = 0.999))
Age3.jz_db.hd.age.seu$gc.like <- ifelse(
  Age3.jz_db.hd.age.seu$gc.score1 > quantile(Age3.jz_db.hd.age.seu$gc.score1, probs = 0.95),'gc.like','other')
table(Age3.jz_db.hd.age.seu$gc.like )
table(Age3.jz_db.hd.age.seu$full_first_type)

ams <- rownames(Age3.jz_db.hd.age.seu@meta.data[Age3.jz_db.hd.age.seu$gc.like == 'gc.like',])
rctd <- rownames( Age3.jz_db.hd.age.seu@meta.data[Age3.jz_db.hd.age.seu$full_first_type == 'GCs',])
all.char <- union(ams,rctd)

N <- length(all.char)
k <- length(intersect(ams,rctd))
n <- length(ams)
m <- length(rctd)
phyper(k-1, m, N-m ,n, lower.tail = T)
phyper(k-1, n, N-n ,m, lower.tail = T)

overlap_table <- matrix(c(k, n - k, m - k, N - n - m + k), nrow = 2, 
                        dimnames = list(c("In Group2", "Not in Group2"), c("In Group1", "Not in Group1")))
fisher_result <- fisher.test(overlap_table)
cat("p-value from Fisher's test:", fisher_result$p.value, "\n")

bin.mat <- data.frame(
  ams = as.numeric(all.char %in% ams),
  rctd = as.numeric(all.char %in% rctd)
)
bin.mat$match <- ifelse(bin.mat$ams == bin.mat$rctd, 'yes', 'no')
table(bin.mat$match)
rownames(bin.mat) <- all.char
cor(bin.mat$ams, bin.mat$rctd,method = 'spearman')
cor(ams,rctd)

jacd <- sum(bin.mat$ams & bin.mat$rctd)/sum(bin.mat$ams | bin.mat$rctd)
jacd


Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = gc.gene, name = 'gc.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'gc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$gc.score1, probs = 0.999))

Young1.jz_db.hd.young.seu$gc.like <- ifelse(
  Young1.jz_db.hd.young.seu$gc.score1 > quantile(Young1.jz_db.hd.young.seu$gc.score1, probs = 0.95),'gc.like','other')
table(Young1.jz_db.hd.young.seu$gc.like )
table(Young1.jz_db.hd.young.seu$full_first_type)

ams <- rownames(Young1.jz_db.hd.young.seu@meta.data[Young1.jz_db.hd.young.seu$gc.like == 'gc.like',])
rctd <- rownames( Young1.jz_db.hd.young.seu@meta.data[Young1.jz_db.hd.young.seu$full_first_type == 'GCs',])
all.char <- union(ams,rctd)

N <- length(all.char)
k <- length(intersect(ams,rctd))
n <- length(ams)
m <- length(rctd)
phyper(k-1, m, N-m ,n, lower.tail = T)
phyper(k-1, n, N-n ,m, lower.tail = T)

overlap_table <- matrix(c(k, n - k, m - k, N - n - m + k), nrow = 2, 
                        dimnames = list(c("In Group2", "Not in Group2"), c("In Group1", "Not in Group1")))
fisher_result <- fisher.test(overlap_table)
cat("p-value from Fisher's test:", fisher_result$p.value, "\n")

bin.mat <- data.frame(
  ams = as.numeric(all.char %in% ams),
  rctd = as.numeric(all.char %in% rctd)
)
bin.mat$match <- ifelse(bin.mat$ams == bin.mat$rctd, 'yes', 'no')
table(bin.mat$match)
rownames(bin.mat) <- all.char
cor(bin.mat$ams, bin.mat$rctd,method = 'spearman')
cor(ams,rctd)

jacd <- sum(bin.mat$ams & bin.mat$rctd)/sum(bin.mat$ams | bin.mat$rctd)
jacd


DefaultAssay(KO2.jz_db.hd.age.seu) <- 'Spatial.008um'
KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = gc.gene, name = 'gc.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'gc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.age.seu$gc.score1, probs = 0.999))

KO2.jz_db.hd.young.seu <- AddModuleScore(KO2.jz_db.hd.young.seu, features = gc.gene, name = 'gc.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'gc.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.young.seu$gc.score1, probs = 0.99))

## SpTs
SpT.gene <- list(rm.lpt.ct.deg[rm.lpt.ct.deg$cluster == 'SpTs',]$gene)
SpT.gene <- list(rm.lpt.ct.deg.top10[rm.lpt.ct.deg.top10$cluster == 'SpTs',]$gene)

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$SpT.score1, probs = 0.999))
Age3.jz_db.hd.age.seu$SpT.like <- ifelse(
  Age3.jz_db.hd.age.seu$SpT.score1 > quantile(Age3.jz_db.hd.age.seu$SpT.score1, probs = 0.6),'SpT.like','other')

table(Age3.jz_db.hd.age.seu$SpT.like )
table(Age3.jz_db.hd.age.seu$full_first_type)

ams <- rownames(Age3.jz_db.hd.age.seu@meta.data[Age3.jz_db.hd.age.seu$gc.like == 'gc.like',])
rctd <- rownames( Age3.jz_db.hd.age.seu@meta.data[Age3.jz_db.hd.age.seu$full_first_type == 'GCs',])
all.char <- union(ams,rctd)

N <- length(all.char)
k <- length(intersect(ams,rctd))
n <- length(ams)
m <- length(rctd)
phyper(k-1, m, N-m ,n, lower.tail = T)
phyper(k-1, n, N-n ,m, lower.tail = T)

overlap_table <- matrix(c(k, n - k, m - k, N - n - m + k), nrow = 2, 
                        dimnames = list(c("In Group2", "Not in Group2"), c("In Group1", "Not in Group1")))
fisher_result <- fisher.test(overlap_table)
cat("p-value from Fisher's test:", fisher_result$p.value, "\n")

bin.mat <- data.frame(
  ams = as.numeric(all.char %in% ams),
  rctd = as.numeric(all.char %in% rctd)
)
bin.mat$match <- ifelse(bin.mat$ams == bin.mat$rctd, 'yes', 'no')
table(bin.mat$match)
rownames(bin.mat) <- all.char
cor(bin.mat$ams, bin.mat$rctd,method = 'spearman')
cor(ams,rctd)
jacd <- sum(bin.mat$ams & bin.mat$rctd)/sum(bin.mat$ams | bin.mat$rctd)
jacd


Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$SpT.score1, probs = 0.999))

DefaultAssay(KO2.jz_db.hd.age.seu) <- 'Spatial.008um'
KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.age.seu$SpT.score1, probs = 0.999))

KO2.jz_db.hd.young.seu <- AddModuleScore(KO2.jz_db.hd.young.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.young.seu$SpT.score1, probs = 0.99))
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'Prl8a9',crop = F,image.alpha = 0.2)


## SpTs
SpT.gene <- list(rm.lpt.ct.deg[rm.lpt.ct.deg$cluster == 'SpTs',]$gene)
SpT.gene <- list(rm.lpt.ct.deg.top10[rm.lpt.ct.deg.top10$cluster == 'SpTs',]$gene)

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$SpT.score1, probs = 0.999))
Age3.jz_db.hd.age.seu$SpT.like <- ifelse(
  Age3.jz_db.hd.age.seu$SpT.score1 > quantile(Age3.jz_db.hd.age.seu$SpT.score1, probs = 0.9),'SpT.like','other')

Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$SpT.score1, probs = 0.999))

DefaultAssay(KO2.jz_db.hd.age.seu) <- 'Spatial.008um'
KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.age.seu$SpT.score1, probs = 0.999))

KO2.jz_db.hd.young.seu <- AddModuleScore(KO2.jz_db.hd.young.seu, features = SpT.gene, name = 'SpT.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'SpT.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.young.seu$SpT.score1, probs = 0.99))
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'Prl8a9',crop = F,image.alpha = 0.2)

## ECs
EC.gene <- list(rm.lpt.ct.deg[rm.lpt.ct.deg$cluster == 'ECs',]$gene)
EC.gene <- list(rm.lpt.ct.deg.top10[rm.lpt.ct.deg.top10$cluster == 'ECs',]$gene)

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = EC.gene, name = 'EC.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'EC.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$EC.score1, probs = 0.999))
Age3.jz_db.hd.age.seu$EC.like <- ifelse(
  Age3.jz_db.hd.age.seu$EC.score1 > quantile(Age3.jz_db.hd.age.seu$EC.score1, probs = 0.9),'EC.like','other')

Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = EC.gene, name = 'EC.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'EC.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$EC.score1, probs = 0.999))

DefaultAssay(KO2.jz_db.hd.age.seu) <- 'Spatial.008um'
KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = EC.gene, name = 'EC.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'EC.score1',crop = F,image.alpha = 0.2)

KO2.jz_db.hd.young.seu <- AddModuleScore(KO2.jz_db.hd.young.seu, features = EC.gene, name = 'EC.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'EC.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.young.seu$EC.score1, probs = 0.99))
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'Prl8a9',crop = F,image.alpha = 0.2)

## S.TGCs
S.TGC.gene <- list(rm.lpt.ct.deg[rm.lpt.ct.deg$cluster == 'S-TGCs',]$gene)
rm.lpt.ct.deg.top12 <- rm.lpt.ct.deg %>% dplyr::group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC)
S.TGC.gene <- list(rm.lpt.ct.deg.top12[rm.lpt.ct.deg.top12$cluster == 'S-TGCs',]$gene)

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = S.TGC.gene, name = 'S.TGC.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'S.TGC.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$S.TGC.score1, probs = 0.999))
Age3.jz_db.hd.age.seu$S.TGC.like <- ifelse(
  Age3.jz_db.hd.age.seu$S.TGC.score1 > quantile(Age3.jz_db.hd.age.seu$S.TGC.score1, probs = 0.9),'S.TGC.like','other')

Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = S.TGC.gene, name = 'S.TGC.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'S.TGC.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(Age3.jz_db.hd.age.seu$S.TGC.score1, probs = 0.999))

DefaultAssay(KO2.jz_db.hd.age.seu) <- 'Spatial.008um'
KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = S.TGC.gene, name = 'S.TGC.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'S.TGC.score1',crop = F,image.alpha = 0.2)

KO2.jz_db.hd.young.seu <- AddModuleScore(KO2.jz_db.hd.young.seu, features = S.TGC.gene, name = 'S.TGC.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'S.TGC.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = quantile(KO2.jz_db.hd.young.seu$S.TGC.score1, probs = 0.99))
SpatialFeaturePlot(KO2.jz_db.hd.young.seu,features = 'Prl8a9',crop = F,image.alpha = 0.2)

SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'Cd68',crop = F,image.alpha = 0.2)

## SASP 

Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = sasp.geneuse, name = 'SASP.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'SASP.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.3)

Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = sasp.geneuse, name = 'SASP.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'SASP.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.3)

KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = sasp.geneuse, name = 'SASP.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'SASP.score1',crop = F,image.alpha = 0.2,max.cutoff = 0.3)

## infm
Age3.jz_db.hd.age.seu <- AddModuleScore(Age3.jz_db.hd.age.seu, features = infm.geneuse, name = 'INFM.score',slot = 'data')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu,features = 'INFM.score1',crop = F,image.alpha = 0.2,max.cutoff = 0.1)

Young1.jz_db.hd.young.seu <- AddModuleScore(Young1.jz_db.hd.young.seu, features = infm.geneuse, name = 'INFM.score',slot = 'data')
SpatialFeaturePlot(Young1.jz_db.hd.young.seu,features = 'INFM.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.1)

KO2.jz_db.hd.age.seu <- AddModuleScore(KO2.jz_db.hd.age.seu, features = infm.geneuse, name = 'INFM.score',slot = 'data')
SpatialFeaturePlot(KO2.jz_db.hd.age.seu,features = 'INFM.score1',crop = F,image.alpha = 0.2,max.cutoff = 0.1)

summary(Age3.jz_db.hd.age.seu$INFM.score1)
summary(Young1.jz_db.hd.young.seu$INFM.score1)
summary(KO2.jz_db.hd.age.seu$INFM.score1)

# 05.to squidpy -----
library(SeuratDisk)
SaveH5Seurat(Age3.jz_db.hd.age.seu, filename = '../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.h5Seurat',overwrite = T)
Convert('../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.h5Seurat',dest = 'h5ad',overwrite = T)

Age3.jz_db.hd.age.seu.tmp <- DietSeurat(Age3.jz_db.hd.age.seu,assays = 'Spatial.008um')
SpatialFeaturePlot(Age3.jz_db.hd.age.seu.tmp,features = 'Cd68',crop = F,image.alpha = 0.2)
SaveH5Seurat(Age3.jz_db.hd.age.seu.tmp, filename = '../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.h5Seurat',overwrite = T)
Convert('../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.h5Seurat',dest = 'h5ad',overwrite = T)

sceasy::convertFormat(Age3.jz_db.hd.age.seu,from = 'seurat',to='anndata',outFile = '../05.10x.hd/03.py.analysis/Age3.jz_db.hd.age.seu.h5ad',assay = 'Spatial.008um')

# 06.dis DSC ----

# Chebyshev: |Δx|≤2 & |Δy|≤2

## age
coords <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_bin <- 1
mac$neighbors_sq <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) <= radius_bin,
           abs(array_row - yi) <= radius_bin) %>%
    pull(Barcode)
})

Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor) <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_sq)),]$cell_type <- "Neibor.DSCs"

ggplot(Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

table(Age3.jz_db.hd.age.seu.5ct$guess.ct)
Age3.jz_db.hd.age.seu.5ct$cell_type_sq <- Age3.jz_db.hd.age.seu.5ct$guess.ct
Age3.jz_db.hd.age.seu.5ct@meta.data[unique(unlist(mac$neighbors_sq)),]$cell_type_sq <- "Near_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_sq)

## young

coords <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_bin <- 1
mac$neighbors_sq <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) <= radius_bin,
           abs(array_row - yi) <= radius_bin) %>%
    pull(Barcode)
})

Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor) <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_sq)),]$cell_type <- "Neibor.DSCs"

ggplot(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

table(Young1.jz_db.hd.young.seu.5ct$guess.ct)
table(Young1.jz_db.hd.young.seu.5ct$full_first_type)
Young1.jz_db.hd.young.seu.5ct$cell_type_sq <- Young1.jz_db.hd.young.seu.5ct$guess.ct
Young1.jz_db.hd.young.seu.5ct@meta.data[unique(unlist(mac$neighbors_sq)),]$cell_type_sq <- "Near_mac_DSCs"
table(Young1.jz_db.hd.young.seu.5ct$cell_type_sq)

## KO
coords <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

radius_bin <- 1
mac$neighbors_sq <- lapply(seq_len(nrow(mac)), function(i) {
  xi <- mac$array_col[i]; yi <- mac$array_row[i]
  dsc %>%
    filter(abs(array_col - xi) <= radius_bin,
           abs(array_row - yi) <= radius_bin) %>%
    pull(Barcode)
})

KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type = KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$guess.ct
rownames(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor) <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$Barcode
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor[unique(unlist(mac$neighbors_sq)),]$cell_type <- "Neibor.DSCs"
table(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor$cell_type)
ggplot(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor, 
       aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type)) + 
  geom_point(size=.02) + scale_color_igv() + theme_bw()

KO2.jz_db.hd.all.seu.5ct$cell_type_sq <- KO2.jz_db.hd.all.seu.5ct$guess.ct
KO2.jz_db.hd.all.seu.5ct@meta.data[unique(unlist(mac$neighbors_sq)),]$cell_type_sq <- "Near_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_sq)
table(KO2.jz_db.hd.all.seu.5ct$guess.ct)

## distance---KO
library(fields)
coords <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop.add.neibor[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

ds_mat  <- as.matrix(dsc[, c("array_col", "array_row")])
mac_mat <- as.matrix(mac[, c("array_col", "array_row")])
dist_mat <- rdist(ds_mat, mac_mat)
rownames(dist_mat) <- dsc$Barcode
colnames(dist_mat) <- mac$Barcode

ko.dsc_dist_mac <- data.frame(dsc_dist_mac = rowSums(dist_mat))
ko.dsc_dist_mac <- ko.dsc_dist_mac[order(ko.dsc_dist_mac$dsc_dist_mac),,drop = F ]

counts <- rowSums(dist_mat <= 5)
dsc$mac_within5bin <- counts
table(dsc$mac_within5bin)

## method1:sum
KO2.jz_db.hd.all.seu.5ct$cell_type_dis_euc <- KO2.jz_db.hd.all.seu.5ct$guess.ct
KO2.jz_db.hd.all.seu.5ct@meta.data[rownames(ko.dsc_dist_mac[1:table(KO2.jz_db.hd.all.seu.5ct$cell_type_sq)[5],,drop = F]),]$cell_type_dis_euc <- "Near_mac_DSCs"
KO2.jz_db.hd.all.seu.5ct@meta.data[rownames(ko.dsc_dist_mac[(nrow(ko.dsc_dist_mac)-table(KO2.jz_db.hd.all.seu.5ct$cell_type_sq)[5]):nrow(ko.dsc_dist_mac) ,,drop = F]),]$cell_type_dis_euc <- "Dis_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_euc)

KO2.jz_db.hd.all.seu.5ct$cell_type_dis_euc <- factor(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_euc,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct$cell_type_dis_euc %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_euc")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_euc)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method2:5bin
KO2.jz_db.hd.all.seu.5ct$cell_type_dis_5bin <- KO2.jz_db.hd.all.seu.5ct$guess.ct
dsc.order <- dsc[order(dsc$mac_within5bin),]
KO2.jz_db.hd.all.seu.5ct@meta.data[dsc.order[(nrow(dsc.order)-333):nrow(dsc.order),]$Barcode,]$cell_type_dis_5bin <- "Near_mac_DSCs"
KO2.jz_db.hd.all.seu.5ct@meta.data[dsc.order[1:1317,]$Barcode,]$cell_type_dis_5bin <- "Dis_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_5bin)

KO2.jz_db.hd.all.seu.5ct$cell_type_dis_5bin <- factor(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_5bin,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct$cell_type_dis_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_5bin)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method3: method1-dis + sq-near
KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_euc <- KO2.jz_db.hd.all.seu.5ct$cell_type_sq

KO2.jz_db.hd.all.seu.5ct@meta.data[setdiff(rownames(ko.dsc_dist_mac[(nrow(ko.dsc_dist_mac)-table(KO2.jz_db.hd.all.seu.5ct$cell_type_sq)[5]):nrow(ko.dsc_dist_mac) ,,drop = F]),rownames(KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_euc =='Near_mac_DSCs',]) ),]$cell_type_dis_sq_euc <- "Dis_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_euc)

KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_euc <- factor(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_euc,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_euc %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_euc")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_euc)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method4: method2-dis + sq-near
# KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_5bin <- KO2.jz_db.hd.all.seu.5ct$cell_type_sq
# 
# KO2.jz_db.hd.all.seu.5ct@meta.data[setdiff(dsc.order[1:1317,]$Barcode,rownames(KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_5bin =='Near_mac_DSCs',]) ),]$cell_type_dis_sq_5bin <- "Dis_mac_DSCs"
table(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_5bin)

# KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_5bin <- factor(KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_5bin,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
# KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>% 
#   ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) + 
#   geom_point(size=.02) + scale_color_tron() + theme_bw() +
#   guides(color = guide_legend(override.aes = list(size = 5)))
KO2.jz_db.hd.all.seu.5ct@meta.data[, c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) + 
  geom_point(size=.02) +
  scale_color_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF','#FF410DFF',rep('#D0DFE6FF',5))) +
  # theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(x = '',y = '') + 
  theme(
    panel.background = element_blank(),  
    plot.background = element_blank(),   
    axis.line = element_line(color = "black"),       
    panel.border = element_rect(color = "black", fill = NA)  
  )

KO2.jz_db.hd.all.seu.5ct@meta.data[KO2.jz_db.hd.all.seu.5ct$cell_type_dis_sq_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) + 
  geom_point(size=.02) +
  scale_color_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF','#FF410DFF')) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  labs(x = '',y = '') + 
  theme(
    panel.background = element_blank(),  
    plot.background = element_blank(),   
    axis.line = element_line(color = "black"),       
    panel.border = element_rect(color = "black", fill = NA)  
  )

#### young
## distance
library(fields)
coords <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

ds_mat  <- as.matrix(dsc[, c("array_col", "array_row")])
mac_mat <- as.matrix(mac[, c("array_col", "array_row")])
dist_mat <- rdist(ds_mat, mac_mat)
rownames(dist_mat) <- dsc$Barcode
colnames(dist_mat) <- mac$Barcode

young.dsc_dist_mac <- data.frame(dsc_dist_mac = rowSums(dist_mat))
young.dsc_dist_mac <- young.dsc_dist_mac[order(young.dsc_dist_mac$dsc_dist_mac),,drop = F ]

counts <- rowSums(dist_mat <= 5)
dsc$mac_within5bin <- counts
table(dsc$mac_within5bin)

## method1:sum
Young1.jz_db.hd.young.seu.5ct$cell_type_dis_euc <- Young1.jz_db.hd.young.seu.5ct$guess.ct
Young1.jz_db.hd.young.seu.5ct@meta.data[rownames(young.dsc_dist_mac[1:table(Young1.jz_db.hd.young.seu.5ct$cell_type_sq)[5],,drop = F]),]$cell_type_dis_euc <- "Near_mac_DSCs"
Young1.jz_db.hd.young.seu.5ct@meta.data[rownames(young.dsc_dist_mac[(nrow(young.dsc_dist_mac)-table(Young1.jz_db.hd.young.seu.5ct$cell_type_sq)[5]):nrow(young.dsc_dist_mac) ,,drop = F]),]$cell_type_dis_euc <- "Dis_mac_DSCs"
table(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_euc)

Young1.jz_db.hd.young.seu.5ct$cell_type_dis_euc <- factor(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_euc,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct$cell_type_dis_euc %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_euc")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_euc)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method2:5bin
Young1.jz_db.hd.young.seu.5ct$cell_type_dis_5bin <- Young1.jz_db.hd.young.seu.5ct$guess.ct
dsc.order <- dsc[order(dsc$mac_within5bin),]
Young1.jz_db.hd.young.seu.5ct@meta.data[dsc.order[(nrow(dsc.order)-2605):nrow(dsc.order),]$Barcode,]$cell_type_dis_5bin <- "Near_mac_DSCs"
Young1.jz_db.hd.young.seu.5ct@meta.data[dsc.order[1:2606,]$Barcode,]$cell_type_dis_5bin <- "Dis_mac_DSCs"
table(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_5bin)

Young1.jz_db.hd.young.seu.5ct$cell_type_dis_5bin <- factor(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_5bin,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct$cell_type_dis_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_5bin)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method3: method1-dis + sq-near
Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_euc <- Young1.jz_db.hd.young.seu.5ct$cell_type_sq

Young1.jz_db.hd.young.seu.5ct@meta.data[setdiff( dsc[dsc$mac_within5bin == '0',]$Barcode,rownames(Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_euc =='Near_mac_DSCs',]) ),]$cell_type_dis_sq_euc <- "Dis_mac_DSCs"
table(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_euc)

Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_euc <- factor(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_euc,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_euc %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_euc")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_euc)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method4: method2-dis + sq-near
# Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_5bin <- Young1.jz_db.hd.young.seu.5ct$cell_type_sq
# 
# Young1.jz_db.hd.young.seu.5ct@meta.data[setdiff(dsc[dsc$mac_within5bin == '0',]$Barcode,rownames(Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_5bin =='Near_mac_DSCs',]) ),]$cell_type_dis_sq_5bin <- "Dis_mac_DSCs"
# table(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_5bin)
# 
# Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_5bin <- factor(Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_5bin,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
# Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>% 
#   ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) + 
#   geom_point(size=.02) + scale_color_tron() + theme_bw() +
#   guides(color = guide_legend(override.aes = list(size = 5)))

Young1.jz_db.hd.young.seu.5ct@meta.data[, c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) + 
  geom_point(size=.02) +
  scale_color_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF','#FF410DFF',rep('#D0DFE6FF',5))) +
  # theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(x = '',y = '') + 
  theme(
    panel.background = element_blank(),  
    plot.background = element_blank(),   
    axis.line = element_line(color = "black"),       
    panel.border = element_rect(color = "black", fill = NA)  
  )

Young1.jz_db.hd.young.seu.5ct@meta.data[Young1.jz_db.hd.young.seu.5ct$cell_type_dis_sq_euc %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_euc")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_euc)) + 
  geom_point(size=.02) +
  scale_color_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF','#FF410DFF')) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  labs(x = '',y = '') + 
  theme(
    panel.background = element_blank(),  
    plot.background = element_blank(),   
    axis.line = element_line(color = "black"),       
    panel.border = element_rect(color = "black", fill = NA)  
  )

#### age
## distance
library(fields)
coords <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop[,c('Barcode',"array_row","array_col","guess.ct")]
mac <- coords %>% filter(guess.ct=="Macrophages.like")
dsc <- coords %>% filter(guess.ct=="DSCs")

ds_mat  <- as.matrix(dsc[, c("array_col", "array_row")])
mac_mat <- as.matrix(mac[, c("array_col", "array_row")])
dist_mat <- rdist(ds_mat, mac_mat)
rownames(dist_mat) <- dsc$Barcode
colnames(dist_mat) <- mac$Barcode

age.dsc_dist_mac <- data.frame(dsc_dist_mac = rowSums(dist_mat))
age.dsc_dist_mac <- age.dsc_dist_mac[order(age.dsc_dist_mac$dsc_dist_mac),,drop = F ]

counts <- rowSums(dist_mat <= 5)
dsc$mac_within5bin <- counts
table(dsc$mac_within5bin)

## method1:sum
Age3.jz_db.hd.age.seu.5ct$cell_type_dis_euc <- Age3.jz_db.hd.age.seu.5ct$guess.ct
Age3.jz_db.hd.age.seu.5ct@meta.data[rownames(age.dsc_dist_mac[1:table(Age3.jz_db.hd.age.seu.5ct$cell_type_sq)[5],,drop = F]),]$cell_type_dis_euc <- "Near_mac_DSCs"
Age3.jz_db.hd.age.seu.5ct@meta.data[rownames(age.dsc_dist_mac[(nrow(age.dsc_dist_mac)-table(Age3.jz_db.hd.age.seu.5ct$cell_type_sq)[5]):nrow(age.dsc_dist_mac) ,,drop = F]),]$cell_type_dis_euc <- "Dis_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_euc)

Age3.jz_db.hd.age.seu.5ct$cell_type_dis_euc <- factor(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_euc,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct$cell_type_dis_euc %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_euc")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_euc)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method2:5bin
Age3.jz_db.hd.age.seu.5ct$cell_type_dis_5bin <- Age3.jz_db.hd.age.seu.5ct$guess.ct
dsc.order <- dsc[order(dsc$mac_within5bin),]
Age3.jz_db.hd.age.seu.5ct@meta.data[dsc.order[(nrow(dsc.order)-4963):nrow(dsc.order),]$Barcode,]$cell_type_dis_5bin <- "Near_mac_DSCs"
Age3.jz_db.hd.age.seu.5ct@meta.data[dsc.order[1:4964,]$Barcode,]$cell_type_dis_5bin <- "Dis_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_5bin)

Age3.jz_db.hd.age.seu.5ct$cell_type_dis_5bin <- factor(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_5bin,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct$cell_type_dis_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_5bin)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method3: method1-dis + sq-near
Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_euc <- Age3.jz_db.hd.age.seu.5ct$cell_type_sq

Age3.jz_db.hd.age.seu.5ct@meta.data[setdiff( dsc[dsc$mac_within5bin == '0',]$Barcode,rownames(Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_euc =='Near_mac_DSCs',]) ),]$cell_type_dis_sq_euc <- "Dis_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_euc)

Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_euc <- factor(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_euc,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_euc %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_euc")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_euc)) + 
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

## method4: method2-dis + sq-near
Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_5bin <- Age3.jz_db.hd.age.seu.5ct$cell_type_sq

Age3.jz_db.hd.age.seu.5ct@meta.data[setdiff(dsc[dsc$mac_within5bin == '0',]$Barcode,rownames(Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_5bin =='Near_mac_DSCs',]) ),]$cell_type_dis_sq_5bin <- "Dis_mac_DSCs"
table(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_5bin)

Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_5bin <- factor(Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_5bin,levels = c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like','ECs','GCs','S-TGCs','SpTs','Unknown'))
Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>%
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) +
  geom_point(size=.02) + scale_color_tron() + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 5)))

Age3.jz_db.hd.age.seu.5ct@meta.data[, c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) + 
  geom_point(size=.02) +
  scale_color_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF','#FF410DFF',rep('#D0DFE6FF',5))) +
  # theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(
    panel.background = element_blank(),  
    plot.background = element_blank(),   
    axis.line = element_line(color = "black"),       
    panel.border = element_rect(color = "black", fill = NA)  
  )

Age3.jz_db.hd.age.seu.5ct@meta.data[Age3.jz_db.hd.age.seu.5ct$cell_type_dis_sq_5bin %in% c('Dis_mac_DSCs','DSCs','Near_mac_DSCs','Macrophages.like'), c("coord_spatial_0.08_x", "coord_spatial_0.08_y", "cell_type_dis_sq_5bin")] %>% 
  ggplot(aes(x=coord_spatial_0.08_x, y=coord_spatial_0.08_y, color=cell_type_dis_sq_5bin)) + 
  geom_point(size=.02) +
  scale_color_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF','#FF410DFF')) +
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  labs(x = '',y = '') + 
  theme(
    panel.background = element_blank(),  
    plot.background = element_blank(),   
    axis.line = element_line(color = "black"),       
    panel.border = element_rect(color = "black", fill = NA)  
  )


# 01.split region -----
Age3.db.8.barcode <- read.csv('../05.10x.hd/02.align/age3_8_db.csv')
Young1.db.8.barcode <- read.csv('../05.10x.hd/02.align/young1_8_db.csv')
KO2.db.8.barcode <- read.csv('../05.10x.hd/02.align/ko2_8_db.csv')

## age3
dim(Age3.hd)
Age3.hd.8um <- Seurat::DietSeurat(Age3.hd,assays = 'Spatial.008um')
Age3.hd.8um$sub_region <- 'LZ'
Age3.hd.8um$Barcode <- rownames(Age3.hd.8um@meta.data)
Age3.hd.8um@meta.data[Age3.jz_db.barcodes$Barcode,]$sub_region <- 'JZ'
Age3.hd.8um@meta.data[Age3.db.8.barcode$Barcode,]$sub_region <- 'DB'
table(Age3.hd.8um$sub_region)
p1 <- SpatialDimPlot(Age3.hd.8um,group.by = 'sub_region',label = F,repel = T,image.alpha = .2,crop = F,alpha = 0.5) + scale_fill_jama()
Age3.hd.8um$sub_region <- factor(Age3.hd.8um$sub_region, levels = c('DB','JZ','LZ'))
View(Age3.hd.8um@meta.data)

## young1 
dim(Young1.hd)
Young1.hd.8um <- Seurat::DietSeurat(Young1.hd,assays = 'Spatial.008um')
Young1.hd.8um$sub_region <- 'LZ'
Young1.hd.8um$Barcode <- rownames(Young1.hd.8um@meta.data)
Young1.hd.8um@meta.data[Young1.jz_db.barcodes$Barcode,]$sub_region <- 'JZ'
Young1.hd.8um@meta.data[Young1.db.8.barcode$Barcode,]$sub_region <- 'DB'
table(Young1.hd.8um$sub_region)
p2 <- SpatialDimPlot(Young1.hd.8um,group.by = 'sub_region',label = F,repel = T,image.alpha = .2,crop = F,alpha = 0.5) + scale_fill_jama()
Young1.hd.8um$sub_region <- factor(Young1.hd.8um$sub_region, levels = c('DB','JZ','LZ'))

## KO2 
dim(KO2.hd)
KO2.hd.8um <- Seurat::DietSeurat(KO2.hd,assays = 'Spatial.008um')
KO2.hd.8um$sub_region <- 'LZ'
KO2.hd.8um$Barcode <- rownames(KO2.hd.8um@meta.data)
KO2.hd.8um@meta.data[ko2.jz_db.barcodes$Barcode,]$sub_region <- 'JZ'
KO2.hd.8um@meta.data[KO2.db.8.barcode$Barcode,]$sub_region <- 'DB'
table(KO2.hd.8um$sub_region)
p3 <- SpatialDimPlot(KO2.hd.8um,group.by = 'sub_region',label = F,repel = T,image.alpha = .2,crop = F,alpha = 0.5) + scale_fill_jama()
KO2.hd.8um$sub_region <- factor(KO2.hd.8um$sub_region, levels = c('DB','JZ','LZ'))

p1 + p2 + p3

# 02.calculate SASP in sub-region ----

## Age
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = sasp.geneuse, name = 'SASP.score',slot = 'data')
p1 <- SpatialFeaturePlot(Age3.hd.8um,features = 'SASP.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.3)
colnames(Age3.hd.8um@meta.data)
plot.data <- Age3.hd.8um@meta.data[,c(15,17)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("LZ", "JZ"),c('JZ','DB'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## Young1
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = sasp.geneuse, name = 'SASP.score',slot = 'data')
p2 <- SpatialFeaturePlot(Young1.hd.8um,features = 'SASP.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.3)
colnames(Young1.hd.8um@meta.data)
plot.data <- Young1.hd.8um@meta.data[,c(19,21)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = sasp.geneuse, name = 'SASP.score',slot = 'data')
p3 <- SpatialFeaturePlot(KO2.hd.8um,features = 'SASP.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.3)
p1 + p2 + p3
colnames(KO2.hd.8um@meta.data)
plot.data <- KO2.hd.8um@meta.data[,c(19,21)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

# 03.calculate inflammation in sub-region ----

## Age
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = infm.geneuse, name = 'INFM.score',slot = 'data')
p1 <- SpatialFeaturePlot(Age3.hd.8um,features = 'INFM.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.1)
colnames(Age3.hd.8um@meta.data)
plot.data <- Age3.hd.8um@meta.data[,c(15,19)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="INFM.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('INFM score')  +
  theme_bw() +
  ylim(-0.1,0.2) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("LZ", "JZ"),c('JZ','DB'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## Young1
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = infm.geneuse, name = 'INFM.score',slot = 'data')
p2 <- SpatialFeaturePlot(Young1.hd.8um,features = 'INFM.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.1)
colnames(Young1.hd.8um@meta.data)
plot.data <- Young1.hd.8um@meta.data[,c(19,23)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="INFM.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('INFM score')  +
  theme_bw() +
  ylim(-0.1,0.25) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = infm.geneuse, name = 'INFM.score',slot = 'data')
p3 <- SpatialFeaturePlot(KO2.hd.8um,features = 'INFM.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.1)
p1 + p2 + p3
colnames(KO2.hd.8um@meta.data)
plot.data <- KO2.hd.8um@meta.data[,c(19,23)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="INFM.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('INFM score')  +
  theme_bw() +
  ylim(-0.1,0.25) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

# 04.calculate of NAD in sub-region -----

## Age
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = nad.geneuse, name = 'NAD.score',slot = 'data')
p1 <- SpatialFeaturePlot(Age3.hd.8um,features = 'NAD.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.4)
colnames(Age3.hd.8um@meta.data)
plot.data <- Age3.hd.8um@meta.data[,c(15,20)]
plot.data %>% dplyr::group_by(sub_region) %>% dplyr::summarise(quantile = median(NAD.score1))
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="NAD.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NAD score')  +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("LZ", "JZ"),c('JZ','DB'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## Young
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = nad.geneuse, name = 'NAD.score',slot = 'data')
p2 <- SpatialFeaturePlot(Young1.hd.8um,features = 'NAD.score1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.4)
colnames(Young1.hd.8um@meta.data)
plot.data <- Young1.hd.8um@meta.data[,c(19,24)]
plot.data %>% dplyr::group_by(sub_region) %>% dplyr::summarise(quantile = median(NAD.score1))
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="NAD.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NAD score')  +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("LZ", "JZ"),c('JZ','DB'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = nad.geneuse, name = 'NAD.score',slot = 'data')
p3 <- SpatialFeaturePlot(KO2.hd.8um,features = 'NAD.score1',crop = F,image.alpha = 0.2,
                   # max.cutoff = 0.1
                   )
p1 + p2 + p3
colnames(KO2.hd.8um@meta.data)
plot.data <- KO2.hd.8um@meta.data[,c(19,24)]
plot.data %>% dplyr::group_by(sub_region) %>% dplyr::summarise(quantile = median(NAD.score1))
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="NAD.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NAD score')  +
  theme_bw() +
  ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

# 05.sub region score -----
## SASP
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,17)], Young1.hd.8um@meta.data[,c(19,21)], KO2.hd.8um@meta.data[,c(19,21)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-SASP
DB.SASP <- plot.data[plot.data$sub_region == 'DB',]
DB.SASP %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(SASP.score1)) # FC:Age-Young 1.344045

p1 <- ggpubr::ggviolin(DB.SASP, x="sample", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                       draw_quantiles = T,
                        alpha = 0.8) +
  ylab('SASP score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

DB.SASP$sample <- factor(DB.SASP$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(DB.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',add = 'boxplot',
                 palette = npg.cols[c(2,1,3)],width = 0.8,
                 add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('') + theme_bw()
p2 <- ggpubr::ggviolin(DB.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',
                       palette = npg.cols[c(2,1,3)],width = 0.8,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('') + theme_bw()
p1 + p2

### JZ
JZ.SASP <- plot.data[plot.data$sub_region == 'JZ',]
JZ.SASP %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(SASP.score1)) # FC:Age-Young 1.055921
p1 <- ggpubr::ggboxplot(JZ.SASP, x="sample", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ggtitle('JZ') +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")
ggpubr::ggviolin(JZ.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',
                 palette = npg.cols[1:3],width = 0.8,
                 add.params = list(fill = 'gray90')) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('')

JZ.SASP$sample <- factor(JZ.SASP$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(JZ.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',add = 'boxplot',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('') + theme_bw()
p2 <- ggpubr::ggviolin(JZ.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('') + theme_bw()
p1 + p2

### LZ
LZ.SASP <- plot.data[plot.data$sub_region == 'LZ',]
LZ.SASP %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(SASP.score1)) # FC:Age-Young 1.156863
p1 <- ggpubr::ggboxplot(LZ.SASP, x="sample", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ggtitle('LZ') +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

ggpubr::ggviolin(LZ.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',
                 palette = npg.cols[1:3],width = 0.8,
                 add.params = list(fill = 'gray90')) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('')

LZ.SASP$sample <- factor(LZ.SASP$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(LZ.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',add = 'boxplot',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('') + theme_bw()
p2 <- ggpubr::ggviolin(LZ.SASP, x = 'sample', y = 'SASP.score1',fill = 'sample',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('SASP score') + xlab('') + theme_bw()
p1 + p2

## INFM
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,19)], Young1.hd.8um@meta.data[,c(19,23)], KO2.hd.8um@meta.data[,c(19,23)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-INFM
DB.INFM <- plot.data[plot.data$sub_region == 'DB',]
DB.INFM %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(INFM.score1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.INFM, x="sample", y="INFM.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('INFM score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

ggpubr::ggviolin(DB.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',
                 palette = npg.cols[1:3],width = 0.8,
                 add.params = list(fill = 'gray90')) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()

DB.INFM$sample <- factor(DB.INFM$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(DB.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',add = 'boxplot',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()
p2 <- ggpubr::ggviolin(DB.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()
p1 + p2

### JZ-INFM
JZ.INFM <- plot.data[plot.data$sub_region == 'JZ',]
JZ.INFM %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(INFM.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.INFM, x="sample", y="INFM.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('INFM score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

ggpubr::ggviolin(JZ.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',outlier.size = .1,
                 palette = npg.cols[1:3],width = 0.8,
                 add.params = list(fill = 'gray90')) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()

JZ.INFM$sample <- factor(JZ.INFM$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(JZ.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',add = 'boxplot',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()
p2 <- ggpubr::ggviolin(JZ.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()
p1 + p2

### LZ-INFM
LZ.INFM <- plot.data[plot.data$sub_region == 'LZ',]
LZ.INFM %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(INFM.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.INFM, x="sample", y="INFM.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('INFM score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

ggpubr::ggviolin(LZ.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',
                 palette = npg.cols[1:3],width = 0.8,
                 add.params = list(fill = 'gray90')) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()

LZ.INFM$sample <- factor(LZ.INFM$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(LZ.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',add = 'boxplot',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()
p2 <- ggpubr::ggviolin(LZ.INFM, x = 'sample', y = 'INFM.score1',fill = 'sample',
                       palette = npg.cols[c(2,1,3)],width = 0.7,
                       add.params = list(fill = 'gray90',alpha=0.6)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "wilcox.test") + 
  ylab('Inflammation reponse score') + xlab('') + theme_bw()
p1 + p2

## NAD
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,20)], Young1.hd.8um@meta.data[,c(19,24)], KO2.hd.8um@meta.data[,c(19,24)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-NAD
DB.NAD <- plot.data[plot.data$sub_region == 'DB',]
DB.NAD %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(NAD.score1)) # FC:Age-Young 1.175182

p1 <- ggpubr::ggboxplot(DB.NAD, x="sample", y="NAD.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NAD score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-NAD
JZ.NAD <- plot.data[plot.data$sub_region == 'JZ',]
JZ.NAD %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(NAD.score1)) # FC:Age-Young 1.035897

p1 <- ggpubr::ggboxplot(JZ.NAD, x="sample", y="NAD.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NAD score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-NAD
LZ.NAD <- plot.data[plot.data$sub_region == 'LZ',]
LZ.NAD %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(NAD.score1)) # FC:Age-Young 0.8743961

p1 <- ggpubr::ggboxplot(JZ.NAD, x="sample", y="NAD.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NAD score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")


## NFKB
colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,28)], Young1.hd.8um@meta.data[,c(19,27)], KO2.hd.8um@meta.data[,c(19,27)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-nfkb
DB.nfkb <- plot.data[plot.data$sub_region == 'DB',]
DB.nfkb %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(nfkb.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.nfkb, x="sample", y="nfkb.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('nfkb score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-nfkb
JZ.nfkb <- plot.data[plot.data$sub_region == 'JZ',]
JZ.nfkb %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(nfkb.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.nfkb, x="sample", y="nfkb.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('nfkb score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-nfkb
LZ.nfkb <- plot.data[plot.data$sub_region == 'LZ',]
LZ.nfkb %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(nfkb.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.nfkb, x="sample", y="nfkb.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('nfkb score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")


# WNT


# 06.DSC signafigure ----
age.dsc <- subset(Age3.jz_db.hd.age.seu.5ct, guess.ct == 'DSCs')
age.dsc <- DietSeurat(age.dsc,assays = 'Spatial.008um')
age.dsc$sample <- 'Age'
young.dsc <- subset(Young1.jz_db.hd.young.seu.5ct, guess.ct == 'DSCs')
young.dsc <- DietSeurat(young.dsc,assays = 'Spatial.008um')
young.dsc$sample <- 'Young'
ko.dsc <- subset(KO2.jz_db.hd.all.seu.5ct, guess.ct == 'DSCs')
ko.dsc <- DietSeurat(ko.dsc,assays = 'Spatial.008um')
ko.dsc$sample <- 'KO'

dsc.8um <- merge(x = age.dsc, y = c(young.dsc,ko.dsc))
dsc.8um[["Spatial.008um"]] <- JoinLayers(dsc.8um[["Spatial.008um"]])
SpatialDimPlot(dsc.8um,crop = F,image.alpha = 0.2)
table(dsc.8um$sample)
dsc.8um$near_samp <- paste0(dsc.8um$sample,'_',dsc.8um$cell_type_mh)

## AUCell
library(AUCell)
expr_mat <- GetAssayData(dsc.8um, assay="Spatial.008um", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

sp.dsc.cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('SASP' = c(sasp.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(sp.dsc.cells_rankings)))
dsc.8um$SASP_AUC <- as.numeric(getAUC(cells_AUC)['SASP',])

sasp.auc <- dsc.8um@meta.data[,c('SASP_AUC','near_samp')] %>% na.omit()
sasp.auc.y_a <- sasp.auc[sasp.auc$sample != 'KO',]
sasp.auc.y_a$sample <- factor(sasp.auc.y_a$sample, levels = c('Young','Age'))
sasp.auc.k_a <- sasp.auc[sasp.auc$sample != 'Young',]

ggplot(sasp.auc.y_a, aes(x=sample, y=SASP_AUC, fill=sample)) +
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'SASP score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  # facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

sasp.auc.y_a$sample <- factor(sasp.auc.y_a$sample, levels = c('Age','Young'))
ggplot(sasp.auc.y_a, aes(x = SASP_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#f69b80','#93aeda')) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

sasp.auc.k_a$sample <- factor(sasp.auc.k_a$sample, levels = c('KO','Age'))
ggplot(sasp.auc.k_a, aes(x = SASP_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80')) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

sasp_near <- sasp.auc[sasp.auc$near_samp == 'Age_Near_mac_DSCs' | sasp.auc$near_samp == 'KO_Near_mac_DSCs' | sasp.auc$near_samp == 'Young_Near_mac_DSCs' , ]
sasp_near$near_samp <- factor(sasp_near$near_samp ,levels = c('KO_Near_mac_DSCs','Age_Near_mac_DSCs','Young_Near_mac_DSCs'))
median(sasp_near[sasp_near$near_samp == 'Age_Near_mac_DSCs',]$SASP_AUC)/median(sasp_near[sasp_near$near_samp == 'Young_Near_mac_DSCs',]$SASP_AUC)
median(sasp_near[sasp_near$near_samp == 'Age_Near_mac_DSCs',]$SASP_AUC)/median(sasp_near[sasp_near$near_samp == 'KO_Near_mac_DSCs',]$SASP_AUC)

ggplot(sasp_near, aes(x = SASP_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.02,0.13) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

sasp_dis <- sasp.auc[sasp.auc$near_samp == 'Young_DSCs' | sasp.auc$near_samp == 'KO_DSCs' | sasp.auc$near_samp == 'Age_DSCs' , ]
sasp_dis$near_samp <- factor(sasp_dis$near_samp ,levels = c('KO_DSCs','Age_DSCs','Young_DSCs'))
median(sasp_dis[sasp_dis$near_samp == 'Age_DSCs',]$SASP_AUC)/median(sasp_dis[sasp_dis$near_samp == 'Young_DSCs',]$SASP_AUC)
median(sasp_dis[sasp_dis$near_samp == 'Age_DSCs',]$SASP_AUC)/median(sasp_dis[sasp_dis$near_samp == 'KO_DSCs',]$SASP_AUC)

ggplot(sasp_dis, aes(x = SASP_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.02,0.13) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")


## NAD
geneset <- list('nad' = c(NAD.meta.gene[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.2 * nrow(sp.dsc.cells_rankings)))
dsc.8um$nad_AUC <- as.numeric(getAUC(cells_AUC)['nad',])

nad.auc <- dsc.8um@meta.data[,c('nad_AUC','near_samp')] %>% na.omit()
nad.auc.y_a <- nad.auc[nad.auc$sample != 'KO',]
nad.auc.y_a$sample <- factor(nad.auc.y_a$sample, levels = c('Young','Age'))
nad.auc.k_a <- nad.auc[nad.auc$sample != 'Young',]

# ggplot(nad.auc.y_a, aes(x=sample, y=nad_AUC, fill=sample)) +
#   geom_boxplot(alpha = .7) + 
#   theme_bw()+
#   labs(x = '', y = 'nad score') + 
#   theme(legend.position = "none") +
#   scale_fill_manual(values = c('#0B996F',"#D6570D")) +
#   # facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
#   ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
#                              method = "wilcox.test")
# 
# nad.auc.y_a$sample <- factor(nad.auc.y_a$sample, levels = c('Age','Young'))
# ggplot(nad.auc.y_a, aes(x = nad_AUC, y = sample, fill = sample)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
#   scale_fill_manual(values = c('#f69b80','#93aeda')) +
#   ggtitle("nad") + xlab('AUC') + ylab("") + 
#   theme_test(base_size = 15) +
#   theme(panel.border=element_rect(linewidth = 1, color = "black"),
#         strip.background = element_rect(linewidth = 1, fill = "white"),
#         strip.text = element_text(size = 18),
#         axis.title.x = element_text(size = 16),
#         axis.text = element_text(size = 16, colour = "black"),
#         axis.line = element_line(size = 0.6), 
#         legend.position ="none")
# 
# nad.auc.k_a$sample <- factor(nad.auc.k_a$sample, levels = c('KO','Age'))
# ggplot(nad.auc.k_a, aes(x = nad_AUC, y = sample, fill = sample)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
#   scale_fill_manual(values = c('#9b93c7','#f69b80')) +
#   ggtitle("nad") + xlab('AUC') + ylab("") + 
#   theme_test(base_size = 15) +
#   theme(panel.border=element_rect(linewidth = 1, color = "black"),
#         strip.background = element_rect(linewidth = 1, fill = "white"),
#         strip.text = element_text(size = 18),
#         axis.title.x = element_text(size = 16),
#         axis.text = element_text(size = 16, colour = "black"),
#         axis.line = element_line(size = 0.6), 
#         legend.position ="none")

nad_near <- nad.auc[nad.auc$near_samp == 'Age_Near_mac_DSCs' | nad.auc$near_samp == 'KO_Near_mac_DSCs' | nad.auc$near_samp == 'Young_Near_mac_DSCs' , ]
nad_near$near_samp <- factor(nad_near$near_samp ,levels = c('KO_Near_mac_DSCs','Age_Near_mac_DSCs','Young_Near_mac_DSCs'))
median(nad_near[nad_near$near_samp == 'Age_Near_mac_DSCs',]$nad_AUC)/median(nad_near[nad_near$near_samp == 'Young_Near_mac_DSCs',]$nad_AUC)
median(nad_near[nad_near$near_samp == 'Age_Near_mac_DSCs',]$nad_AUC)/median(nad_near[nad_near$near_samp == 'KO_Near_mac_DSCs',]$nad_AUC)

ggplot(nad_near, aes(x = nad_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("nad") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0,0.22) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

nad_dis <- nad.auc[nad.auc$near_samp == 'Young_DSCs' | nad.auc$near_samp == 'KO_DSCs' | nad.auc$near_samp == 'Age_DSCs' , ]
nad_dis$near_samp <- factor(nad_dis$near_samp ,levels = c('KO_DSCs','Age_DSCs','Young_DSCs'))
median(nad_dis[nad_dis$near_samp == 'Age_DSCs',]$nad_AUC)/median(nad_dis[nad_dis$near_samp == 'Young_DSCs',]$nad_AUC)
median(nad_dis[nad_dis$near_samp == 'Age_DSCs',]$nad_AUC)/median(nad_dis[nad_dis$near_samp == 'KO_DSCs',]$nad_AUC)

ggplot(nad_dis, aes(x = nad_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("nad") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0,0.22) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")



## INFM
dsc.8um <- AddModuleScore(dsc.8um, features =infm.geneuse, name = 'INFM.score',slot = 'data')

colnames(dsc.8um@meta.data)
plot.data <- dsc.8um@meta.data[,c(30,37)]

plot.data$sample <- factor(plot.data$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(plot.data, x="sample", y="INFM.score1", add = 'boxplot',add.params = list(fill = 'gray90',alpha = 0.6),
                       width = 0.5, 
                       color = "black",#轮廓颜色 
                       fill="sample",#填充
                       palette = npg.cols[c(2,1,3)],
                       xlab = F, #不显示x轴的标签
                       bxp.errorbar=T,#显示误差条
                       bxp.errorbar.width=0.5, #误差条大小
                       size=.5, #箱型图边线的粗细 
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.8) + 
  # geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('Inflammation reponse score')  + 
  theme_bw() + 
  # ylim(0,0.1) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

p2 <- ggpubr::ggviolin(plot.data, x="sample", y="INFM.score1", 
                       # add = 'boxplot',add.params = list(fill = 'gray90',alpha = 0.6),
                       width = 0.5, 
                       color = "black",#轮廓颜色 
                       fill="sample",#填充
                       palette = npg.cols[c(2,1,3)],
                       xlab = F, #不显示x轴的标签
                       bxp.errorbar=T,#显示误差条
                       bxp.errorbar.width=0.5, #误差条大小
                       size=.5, #箱型图边线的粗细 
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.8) + 
  # geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('Inflammation reponse score')  + 
  theme_bw() + 
  # ylim(0,0.1) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p1 + p2

## AUCell

geneset <- list('infm' = c(infm.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(sp.dsc.cells_rankings)))
dsc.8um$infm_AUC <- as.numeric(getAUC(cells_AUC)['infm',])

infm.auc <- dsc.8um@meta.data[,c('infm_AUC','near_samp')] %>% na.omit()

infm_near <- infm.auc[infm.auc$near_samp == 'Age_Near_mac_DSCs' | infm.auc$near_samp == 'KO_Near_mac_DSCs' | infm.auc$near_samp == 'Young_Near_mac_DSCs' , ]
infm_near$near_samp <- factor(infm_near$near_samp ,levels = c('KO_Near_mac_DSCs','Age_Near_mac_DSCs','Young_Near_mac_DSCs'))
median(infm_near[infm_near$near_samp == 'Age_Near_mac_DSCs',]$infm_AUC)/median(infm_near[infm_near$near_samp == 'Young_Near_mac_DSCs',]$infm_AUC)
median(infm_near[infm_near$near_samp == 'Age_Near_mac_DSCs',]$infm_AUC)/median(infm_near[infm_near$near_samp == 'KO_Near_mac_DSCs',]$infm_AUC)

ggplot(infm_near, aes(x = infm_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("infm") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.03,0.1) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

infm_dis <- infm.auc[infm.auc$near_samp == 'Young_DSCs' | infm.auc$near_samp == 'KO_DSCs' | infm.auc$near_samp == 'Age_DSCs' , ]
infm_dis$near_samp <- factor(infm_dis$near_samp ,levels = c('KO_DSCs','Age_DSCs','Young_DSCs'))
median(infm_dis[infm_dis$near_samp == 'Age_DSCs',]$infm_AUC)/median(infm_dis[infm_dis$near_samp == 'Young_DSCs',]$infm_AUC)
median(infm_dis[infm_dis$near_samp == 'Age_DSCs',]$infm_AUC)/median(infm_dis[infm_dis$near_samp == 'KO_DSCs',]$infm_AUC)

ggplot(infm_dis, aes(x = infm_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("infm") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.03,0.1) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")




infm.auc.y_a <- infm.auc[infm.auc$sample != 'KO',]
infm.auc.y_a$sample <- factor(infm.auc.y_a$sample, levels = c('Young','Age'))
infm.auc.k_a <- infm.auc[infm.auc$sample != 'Young',]

ggplot(infm.auc.y_a, aes(x=sample, y=infm_AUC, fill=sample)) +
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'infm score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  # facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

infm.auc.y_a$sample <- factor(infm.auc.y_a$sample, levels = c('Age','Young'))
ggplot(infm.auc.y_a, aes(x = infm_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#f69b80','#93aeda')) +
  ggtitle("Inflammation reponse") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

infm.auc.k_a$sample <- factor(infm.auc.k_a$sample, levels = c('KO','Age'))
ggplot(infm.auc.k_a, aes(x = infm_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80')) +
  ggtitle("Inflammation reponse") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

infm.auc$sample <- factor(infm.auc$sample, levels = c('KO','Age','Young'))
ggplot(infm.auc, aes(x = infm_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("Inflammation reponse") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.045,0.095) + 
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

## NAD
dsc.8um <- AddModuleScore(dsc.8um, features =nad.geneuse, name = 'NAD.score',slot = 'data')

colnames(dsc.8um@meta.data)
plot.data <- dsc.8um@meta.data[,c(30,38)]
plot.data <- plot.data[plot.data$NAD.score1 < 0.19,]
p1 <- ggpubr::ggboxplot(plot.data, x="sample", y="NAD.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  # geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('NAD score')  +
  theme_bw() +
  # ylim(-0.2,0.2) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none')

my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## AUCell

geneset <- list('nad' = c(NAD.meta.gene[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(sp.dsc.cells_rankings)))
dsc.8um$nad_AUC <- as.numeric(getAUC(cells_AUC)['nad',])

nad.auc <- dsc.8um@meta.data[,c('nad_AUC','sample')] %>% na.omit()
nad.auc.y_a <- nad.auc[nad.auc$sample != 'KO',]
nad.auc.y_a$sample <- factor(nad.auc.y_a$sample, levels = c('Young','Age'))
nad.auc.k_a <- nad.auc[nad.auc$sample != 'Young',]

ggplot(nad.auc.y_a, aes(x=sample, y=nad_AUC, fill=sample)) +
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'nad score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  # facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

nad.auc.y_a$sample <- factor(nad.auc.y_a$sample, levels = c('Age','Young'))
ggplot(nad.auc.y_a, aes(x = nad_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#f69b80','#93aeda')) +
  ggtitle("Inflammation reponse") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

nad.auc.k_a$sample <- factor(nad.auc.k_a$sample, levels = c('Age','KO'))
ggplot(nad.auc.k_a, aes(x = nad_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#f69b80','#9b93c7')) +
  ggtitle("Inflammation reponse") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

# zaoxue 
hematop.geneuse

expr_mat <- GetAssayData(merge_hd_8, assay="Spatial.008um", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

sp.all8.cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('hema' = c(hematop.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.all8.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.2 * nrow(sp.all8.cells_rankings)))
merge_hd_8$hema_AUC <- as.numeric(getAUC(cells_AUC)['hema',])

hema.auc <- merge_hd_8@meta.data[,c('hema_AUC','sample')] %>% na.omit()
hema.auc.y_a <- hema.auc[hema.auc$sample != 'KO',]
hema.auc.y_a$sample <- factor(hema.auc.y_a$sample, levels = c('Young','Age'))
hema.auc.k_a <- hema.auc[hema.auc$sample != 'Young',]

ggplot(hema.auc.y_a, aes(x=sample, y=hema_AUC, fill=sample)) +
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'hema score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0B996F',"#D6570D")) +
  # facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

hema.auc.y_a$sample <- factor(hema.auc.y_a$sample, levels = c('Age','Young'))
ggplot(hema.auc.y_a, aes(x = hema_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#f69b80','#93aeda')) +
  ggtitle("Hematopoiesis") + xlab('AUC') + ylab("") +
  theme_test(base_size = 15) +
  xlim(0,0.6) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

hema.auc.k_a$sample <- factor(hema.auc.k_a$sample, levels = c('KO','Age'))
ggplot(hema.auc.k_a, aes(x = hema_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80')) +
  ggtitle("Hematopoiesis") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

hema.auc$sample <- factor(hema.auc$sample, levels = c('KO','Age','Young'))
ggplot(hema.auc, aes(x = hema_AUC, y = sample, fill = sample)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#9b93c7','#f69b80','#93aeda')) +
  ggtitle("Hematopoiesis") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.01,0.5) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none")

## trans noise
library(MASS)
library(ggpubr)
celltypes <- 'DSCs'
Idents(dsc.8um) <- dsc.8um$full_first_type
set.seed(1234)

getEuclideanDistance <- function(celltype, obj, assay, slot, ident1, ident2, group.by, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- subset(obj, cells = WhichCells(obj, idents = celltype))
  
  counts <- GetAssayData(object = tmp, slot = slot, assay = assay)
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) > 0
  expr <- counts[keep_genes, ]
  
  ifelse(min(table(tmp@meta.data[[group.by]])) > 600,
         expr <- expr[,c(rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,]),600)],
                         rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,]),600)])
         ],
         expr <- expr)
  tmp <- subset(tmp,cells = colnames(expr))
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data[[group.by]])[c(ident1,ident2)])
  
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  print(nsample)
  ident2_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident2)], nsample)
  ident1_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident1)], nsample)
  ds_expr_r <- ds_expr[, c(ident1_r, ident2_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, ident1, ident2){
    tmp <- data.matrix(sqrt(matr[genes, ident1]))
    mean <- rowMeans(sqrt(matr[genes, ident1]))
    d_ident1 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident1) <- ident1
    
    tmp <- data.matrix(sqrt(matr[genes, ident2]))
    mean <- rowMeans(sqrt(matr[genes, ident2]))
    d_ident2 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident2) <- ident2
    
    list(ident1 = d_ident1, ident2 = d_ident2)
  }
  ds <- calcEuclDist(matr = ds_expr_r, ident2 = ident2_r, ident1 = ident1_r)
  ds
}

res.dsc <- lapply(celltypes, function(x) getEuclideanDistance(x, 
                                                              obj = dsc.8um,
                                                              assay = 'Spatial.008um',
                                                              slot = 'counts',
                                                              group.by = 'sample',
                                                              ident1 = 'Young',
                                                              ident2 = 'Age',
                                                              lowcv = F))
names(res.dsc) <- celltypes
res.dsc.df <- data.frame(TN.value = unlist(do.call(c, res.dsc))) %>% 
  rownames_to_column(var = 'info') %>% 
  separate(col = info,
           into = c('Celltype','Condition','Sample_cells'),
           sep = '\\.',
           remove = T,
           extra = "merge")
res.dsc.df$Condition <- factor(ifelse(res.dsc.df$Condition == 'ident1','Young','Age'), levels = c('Young','Age'))

my_comparisons <- list(c("Young","Age"))
ggplot(res.dsc.df, aes(x=Condition, y=TN.value, fill=Condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw() +
  labs(x = '', y = 'Transcriptional \nheterogeneity') +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#0b996f',"#d6570d")) +
  facet_wrap(~Celltype, scale="free",nrow = 3) +
  ggtitle('DSCs') + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

ggdensity(res.dsc.df, 
          x = "TN.value",
          add = "median", rug = F,
          color = "Condition", fill = "Condition",
          palette = c('#0B996F',"#D6570D")) +
  ylab("Density") + 
  xlab('Transcriptional \nheterogeneity') + 
  # ggtitle('infm') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) 

## SASP 

dsc.8um <- merge(x = age.dsc, y = c(young.dsc,ko.dsc))
dsc.8um[["Spatial.008um"]] <- JoinLayers(dsc.8um[["Spatial.008um"]])
dsc.8um <- AddModuleScore(dsc.8um, features = sasp.geneuse, name = 'SASP.score',slot = 'data')

colnames(dsc.8um@meta.data)
plot.data <- dsc.8um@meta.data[,c(30,36)]

my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
plot.data$sample <- factor(plot.data$sample, levels = c('Young','Age','KO'))
p1 <- ggpubr::ggviolin(plot.data, x="sample", y="SASP.score1", add = 'boxplot',add.params = list(fill = 'gray90',alpha = 0.6),
                       width = 0.5, 
                       color = "black",#轮廓颜色 
                       fill="sample",#填充
                       palette = npg.cols[c(2,1,3)],
                       xlab = F, #不显示x轴的标签
                       bxp.errorbar=T,#显示误差条
                       bxp.errorbar.width=0.5, #误差条大小
                       size=.5, #箱型图边线的粗细 
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.8) + 
  # geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('SASP score')  + 
  theme_bw() + 
  # ylim(0,0.1) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

p2 <- ggpubr::ggviolin(plot.data, x="sample", y="SASP.score1", 
                       # add = 'boxplot',add.params = list(fill = 'gray90',alpha = 0.6),
                       width = 0.5, 
                       color = "black",#轮廓颜色 
                       fill="sample",#填充
                       palette = npg.cols[c(2,1,3)],
                       xlab = F, #不显示x轴的标签
                       bxp.errorbar=T,#显示误差条
                       bxp.errorbar.width=0.5, #误差条大小
                       size=.5, #箱型图边线的粗细 
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.8) + 
  # geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('SASP score')  + 
  theme_bw() + 
  # ylim(0,0.1) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p1 + p2
# DC geneset ----
DC.geneset <- rio::import_list('./DC.geneset.xlsx')
genage.use <- list(c(str_to_title(str_to_lower(DC.geneset[["GenAge"]][["Symbol"]]))))
infm.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["Inflammatory response genes"]]))))
dna_da.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["DNA damage genes"]]))))
dna_re.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["DNA repair genes"]]))))
dna_jc.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["Cell junction genes"]]))))
ros.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["ROS genes"]]))))
nfkb.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["NF-κB pathway genes"]]))))
wnt.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["WNT pathway genes"]]))))
hed.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["Hedgehog pathway genes"]]))))
bmp.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["BMP pathway genes"]]))))
tgf.use <- list(c(str_to_title(str_to_lower(DC.geneset[["Gene set score analysis"]][["TGF-β pathway genes"]]))))

Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = nfkb.use, name = 'nfkb.use',slot = 'data')
summary(Age3.hd.8um$nfkb.use1)
SpatialFeaturePlot(Age3.hd.8um,features = 'nfkb.use1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.5)

colnames(Age3.hd.8um@meta.data)
plot.data <- Age3.hd.8um@meta.data[,c(15,28)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="nfkb.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('NFKB score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')


Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = nfkb.use, name = 'nfkb.use',slot = 'data')
summary(Young1.hd.8um$genage.score1)
SpatialFeaturePlot(Young1.hd.8um,features = 'nfkb.use1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.5)
colnames(Young1.hd.8um@meta.data)
plot.data <- Young1.hd.8um@meta.data[,c(15,22)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')


KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = nfkb.use, name = 'nfkb.use',slot = 'data')
summary(KO2.hd.8um$nfkb.use1)
SpatialFeaturePlot(KO2.hd.8um,features = 'nfkb.use1',crop = F,image.alpha = 0.2,
                   max.cutoff = 0.5)
colnames(KO2.hd.8um@meta.data)
plot.data <- KO2.hd.8um@meta.data[,c(15,22)]
p1 <- ggpubr::ggboxplot(plot.data, x="sub_region", y="SASP.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sub_region",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('SASP score')  +
  theme_bw() +
  ylim(-0.1,0.4) +
  theme(axis.title.x = element_blank(),legend.position = 'none')



## wnt
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = wnt.use, name = 'wnt.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = wnt.use, name = 'wnt.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = wnt.use, name = 'wnt.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,29)], Young1.hd.8um@meta.data[,c(19,28)], KO2.hd.8um@meta.data[,c(19,28)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-wnt
DB.wnt <- plot.data[plot.data$sub_region == 'DB',]
DB.wnt %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(wnt.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.wnt, x="sample", y="wnt.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('wnt score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-wnt
JZ.wnt <- plot.data[plot.data$sub_region == 'JZ',]
JZ.wnt %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(wnt.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.wnt, x="sample", y="wnt.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('wnt score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-wnt
LZ.wnt <- plot.data[plot.data$sub_region == 'LZ',]
LZ.wnt %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(wnt.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.wnt, x="sample", y="wnt.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('wnt score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## hed
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = hed.use, name = 'hed.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = hed.use, name = 'hed.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = hed.use, name = 'hed.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,30)], Young1.hd.8um@meta.data[,c(19,29)], KO2.hd.8um@meta.data[,c(19,29)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-hed
DB.hed <- plot.data[plot.data$sub_region == 'DB',]
DB.hed %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(hed.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.hed, x="sample", y="hed.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('hed score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-hed
JZ.hed <- plot.data[plot.data$sub_region == 'JZ',]
JZ.hed %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(hed.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.hed, x="sample", y="hed.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('hed score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-hed
LZ.hed <- plot.data[plot.data$sub_region == 'LZ',]
LZ.hed %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(hed.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.hed, x="sample", y="hed.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('hed score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")


## bmp
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = bmp.use, name = 'bmp.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = bmp.use, name = 'bmp.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = bmp.use, name = 'bmp.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,31)], Young1.hd.8um@meta.data[,c(19,30)], KO2.hd.8um@meta.data[,c(19,30)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-bmp
DB.bmp <- plot.data[plot.data$sub_region == 'DB',]
DB.bmp %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(bmp.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.bmp, x="sample", y="bmp.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('bmp score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-bmp
JZ.bmp <- plot.data[plot.data$sub_region == 'JZ',]
JZ.bmp %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(bmp.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.bmp, x="sample", y="bmp.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('bmp score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-bmp
LZ.bmp <- plot.data[plot.data$sub_region == 'LZ',]
LZ.bmp %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(bmp.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.bmp, x="sample", y="bmp.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('bmp score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

## tgf
Age3.hd.8um <- AddModuleScore(Age3.hd.8um, features = tgf.use, name = 'tgf.use',slot = 'data')
Young1.hd.8um <- AddModuleScore(Young1.hd.8um, features = tgf.use, name = 'tgf.use',slot = 'data')
KO2.hd.8um <- AddModuleScore(KO2.hd.8um, features = tgf.use, name = 'tgf.use',slot = 'data')

colnames(Age3.hd.8um@meta.data)
colnames(Young1.hd.8um@meta.data)
colnames( KO2.hd.8um@meta.data)
plot.data <- rbind(Age3.hd.8um@meta.data[,c(15,32)], Young1.hd.8um@meta.data[,c(19,31)], KO2.hd.8um@meta.data[,c(19,31)])
plot.data$sample <- c(Age3.hd.8um$sample, Young1.hd.8um$sample, KO2.hd.8um$sample)

### DB-tgf
DB.tgf <- plot.data[plot.data$sub_region == 'DB',]
DB.tgf %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(tgf.use1)) # FC:Age-Young 1.195652

p1 <- ggpubr::ggboxplot(DB.tgf, x="sample", y="tgf.use1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('tgf score')  +
  ggtitle('DB') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### JZ-tgf
JZ.tgf <- plot.data[plot.data$sub_region == 'JZ',]
JZ.tgf %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(tgf.score1)) # FC:Age-Young 1.030928

p1 <- ggpubr::ggboxplot(DB.tgf, x="sample", y="tgf.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('tgf score')  +
  ggtitle('JZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")

### LZ-tgf
LZ.tgf <- plot.data[plot.data$sub_region == 'LZ',]
LZ.tgf %>% dplyr::group_by(sample) %>% dplyr::summarise(median = median(tgf.score1)) # FC:Age-Young 0.8105727

p1 <- ggpubr::ggboxplot(DB.tgf, x="sample", y="tgf.score1",
                        width = 0.5,
                        color = "black",#轮廓颜色
                        fill="sample",#填充
                        palette = "npg",
                        xlab = F, #不显示x轴的标签
                        bxp.errorbar=T,#显示误差条
                        bxp.errorbar.width=0.5, #误差条大小
                        size=1, #箱型图边线的粗细
                        outlier.shape=NA, #不显示outlier
                        legend = "right",
                        alpha = 0.8) +
  ylab('tgf score')  +
  ggtitle('LZ') +
  theme_bw() +
  # ylim(-0.1,0.1) +
  theme(axis.title.x = element_blank(),legend.position = 'none')
my_comparisons <- list(c("Age", "Young"),c('Age','KO'))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox.test")


library(CellChat)

color.use <- c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b')

Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop 

# Age3
Age3.jz_db_6ct.8um <- subset(Age3.hd.8um,cells = Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop$Barcode)
Age3.jz_db_6ct.8um@meta.data <- Age3.jz_db.hd.age.seu.5ct.filter.na.add.macrop
dim(Age3.jz_db_6ct.8um)
Age3.jz_db_6ct.8um$guess.ct <- factor(Age3.jz_db_6ct.8um$guess.ct, levels = c('DSCs','ECs','GCs','Macrophages.like','S-TGCs','SpTs'))

Idents(Age3.jz_db_6ct.8um) <- Age3.jz_db_6ct.8um$guess.ct
names(color.use) <- levels(Age3.jz_db_6ct.8um)

SpatialDimPlot(Age3.jz_db_6ct.8um,group.by = 'guess.ct',label = F,repel = T,image.alpha = .2,crop = F,alpha = .2, cols = color.use) 

data.input <- GetAssayData(Age3.jz_db_6ct.8um,slot = 'data',assay = 'Spatial.008um')
meta <- data.frame(labels = Idents(Age3.jz_db_6ct.8um),
                   samples = 'Age',
                   row.names = names(Idents(Age3.jz_db_6ct.8um)))
meta$samples <- factor(meta$samples)
spatial.loc <- GetTissueCoordinates(Age3.jz_db_6ct.8um,scale = NULL, cols = c('imagerow','imagecol'))
spatial.loc <- spatial.loc[,-3]
spot.size <- 8
conversion.factor <- 1
spatial.factor <- data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.loc,ratio = spatial.factor$ratio, tol = spatial.factor$tol)
min(d.spatial[d.spatial!=0])

Age.hd.cc <- createCellChat(data.input, meta = meta, group.by = 'labels',
                            datatype = 'spatial',coordinates = spatial.loc,spatial.factors = spatial.factor)

cellchatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(cellchatDB)
Age.hd.cc@DB <- CellChatDB.use

Age.hd.cc <- subsetData(Age.hd.cc)

Age.hd.cc <- Age.hd.cc %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions(variable.both = F)
Age.hd.cc <- Age.hd.cc %>% 
computeCommunProb(type = 'truncatedMean', trim = 0.1, distance.use = T, interaction.range = 250,
                    scale.distance = 0.05,  contact.dependent = T, contact.range = 21) 
Age.hd.cc <- Age.hd.cc %>% 
  filterCommunication(min.cells = 30) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()
  
groupSize <- as.numeric(table(Age.hd.cc@idents))
par(mfrow = c(1,1),xpd = T)
netVisual_circle(Age.hd.cc@net$count, vertex.weight = rowSums(Age.hd.cc@net$count), weight.scale = T,label.edge = T,title.name = 'Number of CCI',color.use = color.use)
netVisual_circle(Age.hd.cc@net$count, vertex.weight = rowSums(Age.hd.cc@net$weight), weight.scale = T,label.edge = T,title.name = 'Weight of CCI',color.use = color.use)

netVisual_heatmap(Age.hd.cc, measure = 'count')
netVisual_heatmap(Age.hd.cc, measure = 'weight')
par(mfrow = c(1,1),xpd = T)
netVisual_aggregate(Age.hd.cc, signaling = 'IGF', layout = 'circle')

netVisual_aggregate(Age.hd.cc, signaling = 'IGF', layout = 'spatial',
                    edge.width.max = .1, vertex.weight.max = .1, vertex.size.max = .1,
                    alpha.image = 0.2, vertex.label.cex = 5,color.use = color.use)
Age.hd.cc <- netAnalysis_computeCentrality(Age.hd.cc, slot.name = 'netP')
netAnalysis_signalingRole_network(Age.hd.cc, signaling = 'IGF', width = 8, height = 3, font.size = 10)

spatialFeaturePlot(Age.hd.cc, features = c('Igf1', 'Igf1r'), point.size = .8, color.heatmap = 'Reds', direction = 1)

## Young 
Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$guess.ct[which(str_detect(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$macrop.like, "^macrop.like"))] <- "Macrophages.like"
table(Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$guess.ct)
Young1.jz_db_6ct.8um <- subset(Young1.hd.8um,cells = Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop$Barcode)
Young1.jz_db_6ct.8um@meta.data <- Young1.jz_db.hd.young.seu.5ct.filter.na.add.macrop
dim(Young1.jz_db_6ct.8um)
Young1.jz_db_6ct.8um$guess.ct <- factor(Young1.jz_db_6ct.8um$guess.ct, levels = c('DSCs','ECs','GCs','Macrophages.like','S-TGCs','SpTs'))

Idents(Young1.jz_db_6ct.8um) <- Young1.jz_db_6ct.8um$guess.ct

SpatialDimPlot(Young1.jz_db_6ct.8um,group.by = 'guess.ct',label = F,repel = T,image.alpha = .2,crop = F,alpha = .2, cols = color.use) 

data.input <- GetAssayData(Young1.jz_db_6ct.8um,slot = 'data',assay = 'Spatial.008um')
meta <- data.frame(labels = Idents(Young1.jz_db_6ct.8um),
                   samples = 'Young',
                   row.names = names(Idents(Young1.jz_db_6ct.8um)))
meta$samples <- factor(meta$samples)
spatial.loc <- GetTissueCoordinates(Young1.jz_db_6ct.8um,scale = NULL, cols = c('imagerow','imagecol'))
spatial.loc <- spatial.loc[,-3]
spot.size <- 8
conversion.factor <- 1
spatial.factor <- data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.loc,ratio = spatial.factor$ratio, tol = spatial.factor$tol)
min(d.spatial[d.spatial!=0])

Young.hd.cc <- createCellChat(data.input, meta = meta, group.by = 'labels',
                            datatype = 'spatial',coordinates = spatial.loc,spatial.factors = spatial.factor)

Young.hd.cc@DB <- CellChatDB.use
Young.hd.cc <- subsetData(Young.hd.cc)

Young.hd.cc <- Young.hd.cc %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions(variable.both = F) %>% 
  computeCommunProb(type = 'truncatedMean', trim = 0.1, distance.use = T, interaction.range = 250,
                    scale.distance = 0.05,  contact.dependent = T, contact.range = 21) %>% 
  filterCommunication(min.cells = 30) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

groupSize <- as.numeric(table(Young.hd.cc@idents))
par(mfrow = c(1,2),xpd = T)
netVisual_circle(Young.hd.cc@net$count, vertex.weight = rowSums(Young.hd.cc@net$count), weight.scale = T,label.edge = T,title.name = 'Number of CCI',color.use = color.use)
netVisual_circle(Young.hd.cc@net$count, vertex.weight = rowSums(Young.hd.cc@net$weight), weight.scale = T,label.edge = F,title.name = 'Weight of CCI',color.use = color.use)

## KO 
KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct[which(str_detect(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$macrop.like, "^macrop.like"))] <- "Macrophages.like"
table(KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$guess.ct)
KO2.jz_db_6ct.8um <- subset(KO2.hd.8um,cells = KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop$Barcode)
KO2.jz_db_6ct.8um@meta.data <- KO2.jz_db.hd.all.seu.5ct.filter.na.add.macrop
dim(KO2.jz_db_6ct.8um)
KO2.jz_db_6ct.8um$guess.ct <- factor(KO2.jz_db_6ct.8um$guess.ct, levels = c('DSCs','ECs','GCs','Macrophages.like','S-TGCs','SpTs'))

Idents(KO2.jz_db_6ct.8um) <- KO2.jz_db_6ct.8um$guess.ct

SpatialDimPlot(KO2.jz_db_6ct.8um,group.by = 'guess.ct',label = F,repel = T,image.alpha = .2,crop = F,alpha = .2, cols = color.use) 

data.input <- GetAssayData(KO2.jz_db_6ct.8um,slot = 'data',assay = 'Spatial.008um')
meta <- data.frame(labels = Idents(KO2.jz_db_6ct.8um),
                   samples = 'KO',
                   row.names = names(Idents(KO2.jz_db_6ct.8um)))
meta$samples <- factor(meta$samples)
spatial.loc <- GetTissueCoordinates(KO2.jz_db_6ct.8um,scale = NULL, cols = c('imagerow','imagecol'))
spatial.loc <- spatial.loc[,-3]
spot.size <- 8
conversion.factor <- 1
spatial.factor <- data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.loc,ratio = spatial.factor$ratio, tol = spatial.factor$tol)
min(d.spatial[d.spatial!=0])

KO.hd.cc <- createCellChat(data.input, meta = meta, group.by = 'labels',
                              datatype = 'spatial',coordinates = spatial.loc,spatial.factors = spatial.factor)

KO.hd.cc@DB <- CellChatDB.use
KO.hd.cc <- subsetData(KO.hd.cc)

KO.hd.cc <- KO.hd.cc %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions(variable.both = F) %>% 
  computeCommunProb(type = 'truncatedMean', trim = 0.1, distance.use = T, interaction.range = 250,
                    scale.distance = 0.05,  contact.dependent = T, contact.range = 21) %>% 
  filterCommunication(min.cells = 30) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

groupSize <- as.numeric(table(KO.hd.cc@idents))
par(mfrow = c(1,2),xpd = T)
netVisual_circle(KO.hd.cc@net$count, vertex.weight = rowSums(KO.hd.cc@net$count), weight.scale = T,label.edge = T,title.name = 'Number of CCI',color.use = color.use)
netVisual_circle(KO.hd.cc@net$count, vertex.weight = rowSums(KO.hd.cc@net$weight), weight.scale = T,label.edge = F,title.name = 'Weight of CCI',color.use = color.use)

# merge age-young
object.list <- list( Young = Young.hd.cc,Age = Age.hd.cc)
age_young.hd.cc <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(age_young.hd.cc, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(age_young.hd.cc, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1), xpd=TRUE)
netVisual_diffInteraction(age_young.hd.cc, weight.scale = T,label.edge = T,sources.use = 'Macrophages.like')
netVisual_diffInteraction(age_young.hd.cc, weight.scale = T, measure = "weight",label.edge = T)

gg1 <- netVisual_heatmap(age_young.hd.cc)
gg2 <- netVisual_heatmap(age_young.hd.cc, measure = "weight")
gg1 + gg2
netVisual_heatmap(Age.hd.cc, color.heatmap = "RdBu",measure = 'weight')

# merge age-ko
object.list <- list( Age = Age.hd.cc,KO = KO.hd.cc)
age_ko.hd.cc <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
gg1 <- compareInteractions(age_ko.hd.cc, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(age_ko.hd.cc, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1), xpd=TRUE)
netVisual_diffInteraction(age_ko.hd.cc, weight.scale = T,label.edge = T)
netVisual_diffInteraction(age_ko.hd.cc, weight.scale = T, measure = "weight",label.edge = T)

gg1 <- netVisual_heatmap(age_ko.hd.cc)
gg2 <- netVisual_heatmap(age_ko.hd.cc, measure = "weight")
gg1 + gg2
gg1 <- rankNet(age_ko.hd.cc, mode = "comparison", measure = "weight", sources.use = 'DSCs', targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(age_ko.hd.cc, mode = "comparison", measure = "weight", sources.use = 'DSCs', targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

p1 <- netVisual_heatmap(age_young.hd.cc, title.name = 'Age vs Young number')
p2 <- netVisual_heatmap(age_ko.hd.cc, title.name = 'Age vs KO number')
p1 + p2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(age_ko.hd.cc, weight.scale = T,label.edge = T,sources.use = 'DSCs',title.name = 'Age vs KO')
netVisual_diffInteraction(age_young.hd.cc, weight.scale = T,label.edge = T,sources.use = 'DSCs',title.name = 'Age vs Young')

netVisual_heatmap(age_ko.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like',measure = 'weight')
netVisual_heatmap(age_young.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like',measure = 'weight')

netVisual_heatmap(age_ko.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like')
netVisual_heatmap(age_young.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like')

plot_heatmap_fixed_range(age_ko.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like')
plot_heatmap_fixed_range(age_young.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like')

plot_heatmap_fixed_range(age_ko.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like',measure = 'weight')
plot_heatmap_fixed_range(age_young.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like',measure = 'weight')

plot_cellchat_heatmap_fixed_0_15(Age.hd.cc, title.name = 'Aged',sources.use = 'Macrophages.like')
plot_cellchat_heatmap_fixed_0_15(Young.hd.cc, title.name = 'Young',sources.use = 'Macrophages.like')
plot_cellchat_heatmap_fixed_0_15(Age.hd.cc, title.name = 'Aged',sources.use = 'Macrophages.like',measure = 'weight')
netVisual_heatmap(Young.hd.cc, title.name = 'Young',sources.use = 'Macrophages.like',measure = 'weight')
netVisual_heatmap(Age.hd.cc, title.name = 'Young',sources.use = 'Macrophages.like',measure = 'weight')
plot_cellchat_heatmap_fixed_0_15(Young.hd.cc, title.name = 'Young',sources.use = 'Macrophages.like',measure = 'weight')


netVisual_heatmap(age_young.hd.cc, title.name = 'Age vs Young number',sources.use = 'Macrophages.like')

age.dsc <- subset(Age3.jz_db.hd.age.seu.5ct, guess.ct == 'DSCs')
age.dsc <- DietSeurat(age.dsc,assays = 'Spatial.008um')
age.dsc$sample <- 'Age'
young.dsc <- subset(Young1.jz_db.hd.young.seu.5ct, guess.ct == 'DSCs')
young.dsc <- DietSeurat(young.dsc,assays = 'Spatial.008um')
young.dsc$sample <- 'Young'
ko.dsc <- subset(KO2.jz_db.hd.all.seu.5ct, guess.ct == 'DSCs')
ko.dsc <- DietSeurat(ko.dsc,assays = 'Spatial.008um')
ko.dsc$sample <- 'KO'

dsc.8um <- merge(x = age.dsc, y = c(young.dsc,ko.dsc))
dsc.8um[["Spatial.008um"]] <- JoinLayers(dsc.8um[["Spatial.008um"]])
SpatialDimPlot(dsc.8um,crop = F,image.alpha = 0.2)
table(dsc.8um$sample)
dsc.8um$near_samp <- paste0(dsc.8um$sample,'_',dsc.8um$cell_type_dis_sq_5bin)
table(dsc.8um$near_samp)

# young ----
library(AUCell)
expr_mat <- GetAssayData(young.dsc, assay="Spatial.008um", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

young.dsc.cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('SASP' = c(sasp.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            young.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(young.dsc.cells_rankings)))
young.dsc$SASP_AUC <- as.numeric(getAUC(cells_AUC)['SASP',])
table(young.dsc$cell_type_dis_sq_5bin)

sasp.auc <- young.dsc@meta.data[,c('SASP_AUC','cell_type_dis_sq_5bin')] %>% na.omit()
table(sasp.auc$cell_type_dis_sq_5bin)
median(sasp.auc[sasp.auc$cell_type_dis_sq_5bin == 'Near_mac_DSCs',]$SASP_AUC)/median(sasp.auc[sasp.auc$cell_type_dis_sq_5bin == 'DSCs',]$SASP_AUC)

ggplot(sasp.auc, aes(x = SASP_AUC, y = cell_type_dis_sq_5bin, fill = cell_type_dis_sq_5bin)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF')) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  xlim(0.03,0.12) + 
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 
# + 
#   ggpubr::stat_compare_means(comparisons = list(c('DSCs','Near_mac_DSCs')),paired = F,
#                              method = "wilcox.test")

geneset <- list('infm' = c(infm.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            young.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.05 * nrow(young.dsc.cells_rankings)))
young.dsc$infm_AUC <- as.numeric(getAUC(cells_AUC)['infm',])

infm.auc <- young.dsc@meta.data[,c('infm_AUC','cell_type_dis_sq_5bin')] %>% na.omit()
table(infm.auc$cell_type_dis_sq_5bin)
median(infm.auc[infm.auc$cell_type_dis_sq_5bin == 'Near_mac_DSCs',]$infm_AUC)/median(infm.auc[infm.auc$cell_type_dis_sq_5bin == 'DSCs',]$infm_AUC)

ggplot(infm.auc, aes(x = infm_AUC, y = cell_type_dis_sq_5bin, fill = cell_type_dis_sq_5bin)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF')) +
  ggtitle("infm") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.02,0.08)+
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

# age -----
library(AUCell)
expr_mat <- GetAssayData(age.dsc, assay="Spatial.008um", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

age.dsc.cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('SASP' = c(sasp.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            age.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(age.dsc.cells_rankings)))
age.dsc$SASP_AUC <- as.numeric(getAUC(cells_AUC)['SASP',])
table(age.dsc$cell_type_dis_sq_5bin)

sasp.auc <- age.dsc@meta.data[,c('SASP_AUC','cell_type_dis_sq_5bin')] %>% na.omit()
table(sasp.auc$cell_type_dis_sq_5bin)
median(sasp.auc[sasp.auc$cell_type_dis_sq_5bin == 'Near_mac_DSCs',]$SASP_AUC)/median(sasp.auc[sasp.auc$cell_type_dis_sq_5bin == 'DSCs',]$SASP_AUC)

ggplot(sasp.auc, aes(x = SASP_AUC, y = cell_type_dis_sq_5bin, fill = cell_type_dis_sq_5bin)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF')) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  xlim(0.03,0.14) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 
# + 
#   ggpubr::stat_compare_means(comparisons = list(c('DSCs','Near_mac_DSCs')),paired = F,
#                              method = "wilcox.test")

geneset <- list('infm' = c(infm.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            age.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.05 * nrow(age.dsc.cells_rankings)))
age.dsc$infm_AUC <- as.numeric(getAUC(cells_AUC)['infm',])

infm.auc <- age.dsc@meta.data[,c('infm_AUC','cell_type_dis_sq_5bin')] %>% na.omit()
table(infm.auc$cell_type_dis_sq_5bin)
median(infm.auc[infm.auc$cell_type_dis_sq_5bin == 'Near_mac_DSCs',]$infm_AUC)/median(infm.auc[infm.auc$cell_type_dis_sq_5bin == 'DSCs',]$infm_AUC)

ggplot(infm.auc, aes(x = infm_AUC, y = cell_type_dis_sq_5bin, fill = cell_type_dis_sq_5bin)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF')) +
  ggtitle("infm") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.025,0.08)+
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

# ko -----
library(AUCell)
expr_mat <- GetAssayData(ko.dsc, assay="Spatial.008um", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

ko.dsc.cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

geneset <- list('SASP' = c(sasp.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            ko.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(ko.dsc.cells_rankings)))
ko.dsc$SASP_AUC <- as.numeric(getAUC(cells_AUC)['SASP',])
table(ko.dsc$cell_type_dis_sq_5bin)

sasp.auc <- ko.dsc@meta.data[,c('SASP_AUC','cell_type_dis_sq_5bin')] %>% na.omit()
table(sasp.auc$cell_type_dis_sq_5bin)
median(sasp.auc[sasp.auc$cell_type_dis_sq_5bin == 'Near_mac_DSCs',]$SASP_AUC)/median(sasp.auc[sasp.auc$cell_type_dis_sq_5bin == 'DSCs',]$SASP_AUC)

ggplot(sasp.auc, aes(x = SASP_AUC, y = cell_type_dis_sq_5bin, fill = cell_type_dis_sq_5bin)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF')) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  # xlim(0.03,0.14) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 
# + 
#   ggpubr::stat_compare_means(comparisons = list(c('DSCs','Near_mac_DSCs')),paired = F,
#                              method = "wilcox.test")

geneset <- list('infm' = c(infm.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            ko.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.05 * nrow(ko.dsc.cells_rankings)))
ko.dsc$infm_AUC <- as.numeric(getAUC(cells_AUC)['infm',])

infm.auc <- ko.dsc@meta.data[,c('infm_AUC','cell_type_dis_sq_5bin')] %>% na.omit()
table(infm.auc$cell_type_dis_sq_5bin)
median(infm.auc[infm.auc$cell_type_dis_sq_5bin == 'Near_mac_DSCs',]$infm_AUC)/median(infm.auc[infm.auc$cell_type_dis_sq_5bin == 'DSCs',]$infm_AUC)

ggplot(infm.auc, aes(x = infm_AUC, y = cell_type_dis_sq_5bin, fill = cell_type_dis_sq_5bin)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#6EE2FFFF','#D0DFE6FF','#F7C530FF')) +
  ggtitle("infm") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.02,0.08)+
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

# all dsc -----

## not integrate
expr_mat <- GetAssayData(dsc.8um, assay="Spatial.008um", layer="data")
rownames(expr_mat)[1:5]
colnames(expr_mat)[1:5]
is.null(rownames(expr_mat))
is.null(colnames(expr_mat))

sp.dsc.cells_rankings <- AUCell_buildRankings(expr_mat,nCores = 30)

## sasp
geneset <- list('SASP' = c(sasp.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.1 * nrow(sp.dsc.cells_rankings)))
dsc.8um$SASP_AUC <- as.numeric(getAUC(cells_AUC)['SASP',])

sasp.auc <- dsc.8um@meta.data[,c('SASP_AUC','near_samp')] %>% na.omit()
table(sasp.auc$near_samp)
sasp.auc_near <- sasp.auc[sasp.auc$near_samp %in% c('Age_Near_mac_DSCs','KO_Near_mac_DSCs','Young_Near_mac_DSCs'),]
table(sasp.auc_near$near_samp)
sasp.auc_near$near_samp <- factor(sasp.auc_near$near_samp, levels = c('KO_Near_mac_DSCs','Age_Near_mac_DSCs','Young_Near_mac_DSCs'))

ggplot(sasp.auc_near, aes(x = SASP_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = (c('#f69b80','red3','#9b93c7'))) +
  ggtitle("SASP") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.02,0.14)+
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

## infm
geneset <- list('infm' = c(infm.geneuse[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.05 * nrow(sp.dsc.cells_rankings)))
dsc.8um$infm_AUC <- as.numeric(getAUC(cells_AUC)['infm',])

infm.auc <- dsc.8um@meta.data[,c('infm_AUC','near_samp')] %>% na.omit()
table(infm.auc$near_samp)
infm.auc_near <- infm.auc[infm.auc$near_samp %in% c('Age_Near_mac_DSCs','KO_Near_mac_DSCs','Young_Near_mac_DSCs'),]
table(infm.auc_near$near_samp)
infm.auc_near$near_samp <- factor(infm.auc_near$near_samp, levels = c('KO_Near_mac_DSCs','Age_Near_mac_DSCs','Young_Near_mac_DSCs'))

ggplot(infm.auc_near, aes(x = infm_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#f69b80','red3','#9b93c7')) +
  ggtitle("infm") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  xlim(0.02,0.08)+
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 

## nad
geneset <- list('nad' = c(NAD.dehy.comp.gene[[1]]))
cells_AUC <- AUCell_calcAUC(geneset,                     
                            sp.dsc.cells_rankings,                 
                            nCores = 30,               
                            aucMaxRank = ceiling(0.05 * nrow(sp.dsc.cells_rankings)))

dsc.8um$nad_AUC <- as.numeric(getAUC(cells_AUC)['nad',])

nad.auc <- dsc.8um@meta.data[,c('nad_AUC','near_samp')] %>% na.omit()
table(nad.auc$near_samp)
nad.auc_near <- nad.auc[nad.auc$near_samp %in% c('Age_Near_mac_DSCs','KO_Near_mac_DSCs','Young_Near_mac_DSCs'),]
table(nad.auc_near$near_samp)
nad.auc_near$near_samp <- factor(nad.auc_near$near_samp, levels = c('KO_Near_mac_DSCs','Age_Near_mac_DSCs','Young_Near_mac_DSCs'))

ggplot(nad.auc_near, aes(x = nad_AUC, y = near_samp, fill = near_samp)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.8) +
  scale_fill_manual(values = c('#f69b80','red3','#9b93c7')) +
  ggtitle("nad") + xlab('AUC') + ylab("") + 
  theme_test(base_size = 15) +
  # xlim(0.02,0.08)+
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6), 
        legend.position ="none") 
