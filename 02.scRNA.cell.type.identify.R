# intergrate data ------
plan(multisession, workers=40)
Clean_sct.inte <- Clean_qc.merge.filtered %>%
 SCTransform(vars.to.regress = c('mitoRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1',split.by = 'condition')
VlnPlot(Clean_sct.inte, group.by = 'SCT_snn_res.1', pt.size = 0,features = c('nFeature_RNA','nCount_RNA','mitoRatio','rpRatio'))
DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.2',split.by = 'condition')

FeaturePlot(Clean_sct.inte,features = c('nFeature_RNA','nCount_RNA','mitoRatio','rpRatio'))
summary(Clean_sct.inte@meta.data[Clean_sct.inte$SCT_snn_res.1 == '5',]$nCount_RNA)
summary(Clean_sct.inte@meta.data[Clean_sct.inte$SCT_snn_res.1 == '5',]$nFeature_RNA)
summary(Clean_sct.inte$nCount_RNA)
summary(Clean_sct.inte$nFeature_RNA)
table(Clean_sct.inte@meta.data[Clean_sct.inte$SCT_snn_res.1 == '5',]$sample)

Idents(Clean_sct.inte) <- Clean_sct.inte$SCT_snn_res.0.2
# Clean_sct.inte <- FindSubCluster(Clean_sct.inte, cluster = '4',graph.name = 'SCT_snn')
# DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'sub.cluster') + scale_color_igv()

FeaturePlot(Clean_sct.inte, features = c('Pparg','Ppara', # 0
                                           'Ghr','Tshr', # 0
                                           'Pdgfra','Fbn1', # 7
                                           'Pecam1','Ptprb', # 9
                                           'Reln','Abcc9', # 9
                                           'Kcnq5','Notch3', # 
                                           'Myh8','Myh11',
                                           'Nkain2','Grid2',
                                           'Hbb-bs','Hba-a1',
                                           'Adgre1','Apoe',  # 1 11 8
                                           'Dock2','Pax5', # B 7 12 
                                           'Krt17','Krt14', # Keratinocytes 14
                                           'Skap1','Ms4a4b', # T memory 4 
                                           'Upk3b','Nkain4','Msln', # Mesothelial  2
                                           'C1qc', 'Mrc1', # marcop 1
                                           'H2-Eb1','Ciita', # DCs 5
                                           'Sirpb1c', # mono 8
                                           'S100a9', # neut 11
                                           'Fscn1' # DC pro 15
))

# cell types identify----
Idents(Clean_sct.inte) <- Clean_sct.inte$SCT_snn_res.1
DefaultAssay(Clean_sct.inte) <- 'SCT'
Clean_sct.inte <- PrepSCTFindMarkers(Clean_sct.inte)
deg.res.1 <- FindAllMarkers(Clean_sct.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1) %>% 
  dplyr::filter(p_val_adj < 0.05)
deg.res.1.top10 <- deg.res.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

FeaturePlot(Clean_sct.inte, 
            features = c('Lrp1b'),cols = c('gray90','red3'),split.by = 'condition',
            order = T)
## major cell type
# endo 20  LEC 25
FeaturePlot(Clean_sct.inte, 
            features = c('Pecam1','Kdr','Tek','Vwf','Cldn5'),
            order = T)
# Adipocyte 1 2 12 13 or 10
FeaturePlot(Clean_sct.inte, 
            features = c('Pparg','Ppara', 
                         'Ghr','Tshr'),
            order = T)
# FAP 7
FeaturePlot(Clean_sct.inte, 
            features = c('Pdgfra','Fbn1'),
            order = T)
# PC/SMC 28
FeaturePlot(Clean_sct.inte, 
            features = c( 'Kcnq5','Notch3','Abcc9','Myh8','Mybpc1','Rgs5'),
            order = T)

# epdi cell 21 22
FeaturePlot(Clean_sct.inte, 
            features = c('Ank3','Abcb5','Pax2'),
            order = T)
# sch cell 10
FeaturePlot(Clean_sct.inte, 
            features = c('Chl1','Grid2','Nrxn3'),
            order = T)
# meso cell 0 6
FeaturePlot(Clean_sct.inte, 
            features = c('Cdh11','Adamtsl1','Upk3b'),
            order = T)

# macrop 3458
FeaturePlot(Clean_sct.inte, 
            features = c('Ptprc','Adgre1','C1qa','Cd80','Mrc1','Il6'),
            order = T)
# neut NA
FeaturePlot(Clean_sct.inte, 
            features = c('S100a9','S100a8'),
            order = T)
# DC 11 26
FeaturePlot(Clean_sct.inte, 
            features = c('H2-Aa','H2-Eb1'),
            order = T)
# mono 17
FeaturePlot(Clean_sct.inte, 
            features = c('Adgre4','Pou2f2'),
            order = T)
# B 7 18 19
FeaturePlot(Clean_sct.inte, 
            features = c('Bank1','Cd79b','Cd38','Top2a','Jchain'),
            order = T)
# T 9 16 23 27
FeaturePlot(Clean_sct.inte, 
            features = c('Ms4a4b','Skap1'),
            order = T)
# Epididymal ductal cells 24
FeaturePlot(Clean_sct.inte, 
            features = c('Acta2','Krt14','Hydin','Spef2'),
            order = T)

Clean_sct.inte <- RenameIdents(Clean_sct.inte,
                               '1' = 'Adipocytes',
                               '2' = 'Adipocytes',
                               '12' = 'Adipocytes',
                               '13' = 'Adipocytes',
                               '14' = 'FAPs',
                               '15' = 'FAPs',
                               '0' = 'MCs',
                               '6' = 'MCs',
                               '20' = 'ECs',
                               '25' = 'LECs',
                               '28' = 'Pericytes',
                               '10' = 'Schwann cells',
                               '21' = 'Epididymal cells',
                               '22' = 'Epididymal cells',
                               '24' = 'Epididymal ductal cells',
                               '3' = 'Macrophages',
                               '4' = 'Macrophages',
                               '5' = 'Macrophages',
                               '8' = 'Macrophages',
                               '11' = 'DCs',
                               '26' = 'DCs',
                               '17' = 'Monocytes',
                               '7' = 'B cells',
                               '18' = 'B cells',
                               '19' = 'B cells',
                               '9' = 'T cells',
                               '16' = 'T cells',
                               '23' = 'T cells',
                               '27' = 'T cells'
                               )
Clean_sct.inte$cell.type<- Idents(Clean_sct.inte)
DimPlot(Clean_sct.inte, group.by = 'cell.type')

# remove Epididymal ductal cells -----
Clean_sct.inte.rm.edc <- subset(Clean_sct.inte, cells = 
                                  rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type == 'Epididymal ductal cells',]), invert = T)
Clean_sct.inte.rm.edc <- RunUMAP(Clean_sct.inte.rm.edc,dims = 1:30)
DimPlot(Clean_sct.inte.rm.edc,group.by = 'cell.type',label = T, cols = c(npg.cols, use.cols))

Clean_sct.inte.rm.edc <- subset(Clean_qc.merge.filtered, cells = rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type != 'Epididymal ductal cells',]))
Clean_sct.inte.rm.edc <- Clean_sct.inte.rm.edc %>%
  SCTransform(vars.to.regress = c('mitoRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

DimPlot(Clean_sct.inte.rm.edc, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.5') + scale_color_igv()
Clean_sct.inte.rm.edc$cell.type.raw <- Clean_sct.inte@meta.data[match(colnames(Clean_sct.inte.rm.edc), colnames(Clean_sct.inte)),]$cell.type
DimPlot(Clean_sct.inte.rm.edc, reduction = "umap", label = T,group.by = 'cell.type.raw')
FeaturePlot(Clean_sct.inte.rm.edc, 
            features = c('Adgre4','Pou2f2','Chl1','Nrxn3'),
            order = T)
Clean_sct.inte.rm.edc$cell.type.percise.new <- Clean_sct.inte.rm.edc$cell.type.raw
