
########################################################################## bulk RNA-seq section ############################################################

# 01.data import -----

# gene info 
mus.gene.info <- read.table('../mouse.110.gtf',header = F,sep = '\t')
colnames(mus.gene.info) <- c('geneid','symbol','biotype')
mus.gene.info[mus.gene.info$symbol == 'ensembl',]$symbol <- mus.gene.info[mus.gene.info$symbol == 'ensembl',]$geneid
mus.gene.info[mus.gene.info$symbol == 'ensembl_havana',]$symbol <- mus.gene.info[mus.gene.info$symbol == 'ensembl_havana',]$geneid
# counts
rna.counts <- read.table('../04.rnaseq/02.quant/all_counts.txt',header = T,sep = '\t') %>% rownames_to_column(var = 'geneid')
colnames(rna.counts)

samp.name <- c(paste0('Young',c(10,seq(1,9))), 
               paste0('Aged+78c',c(10,seq(1,9))),
               paste0('Aged+NR',c(10,seq(1,9))),
               paste0('Aged',c(10,seq(1,9))))

colnames(rna.counts) <- c('geneid',samp.name)
rna.counts.pcg <- rna.counts %>% 
  left_join(y = mus.gene.info, by = 'geneid') %>% 
  na.omit() %>% 
  dplyr::filter(biotype == 'protein_coding') %>% 
  dplyr::group_by(`symbol`) %>% slice_sample(n = 1) %>% 
  column_to_rownames(var = 'symbol') %>% 
  dplyr::select(!c('geneid','biotype')) %>% 
  round(0)

# TPM 
rna.tpm <- read.table('../04.rnaseq/02.quant/all_quant_TPM.txt',header = T,sep = '\t') %>% rownames_to_column(var = 'geneid')
colnames(rna.tpm) <- c('geneid',samp.name)
rna.tpm.pcg <- rna.tpm %>% 
  left_join(y = mus.gene.info, by = 'geneid') %>% 
  na.omit() %>% 
  dplyr::filter(biotype == 'protein_coding') %>% 
  dplyr::group_by(`symbol`) %>% slice_sample(n = 1) %>% 
  column_to_rownames(var = 'symbol') %>% 
  dplyr::select(!c('geneid','biotype')) %>% round(2)


my.samp <- c(paste0('Aged',c(1,2,3,5,6,7)),
               paste0('Young',c(1,4,5,8,9,10)),
               paste0('Aged+NR',c(1,3,5,7,8,10)),
               paste0('Aged+78c',c(1,2,3,4,6,9)))

rna.counts.clean <- rna.counts.pcg[pick.genes,]
rna.counts.clean <- rna.counts.clean[,which(names(rna.counts.clean) %in% c(pick.samp))]
rna.tpm.clean <- rna.tpm.clean[,which(names(rna.counts.clean) %in% c(pick.samp))]

colData = data.frame(row.names = colnames(rna.counts.clean),group_list = c(rep('Young',6),
                                                                           rep('Aged+78c',6),
                                                                           rep('Aged+NR',6),
                                                                           rep('Aged',6)))
dds <- DESeqDataSetFromMatrix(countData = rna.counts.clean,colData = colData,design = ~group_list)

rld_rlog <- rlog(dds, blind=TRUE,fitType = "glmGamPoi")
rlog_mat <- assay(rld_rlog) %>% as.data.frame()

colData$group_list <- factor(colData$group_list, levels = levels(group_list_large))
d_p(rlog_mat,group_list = colData$group_list,'') + 
  labs(x = 'PC1 (21.9%)', y = 'PC2 (13.9%)') + 
  theme_bw() + 
  theme(axis.title = element_text(size = 15))

rlog_cor <- rlog_mat %>%
  cor(method = 'spearman',use = 'pairwise.complete.ob')
ComplexHeatmap::pheatmap(rlog_cor,
                         name = "Spearman's coeff",
                         border_color = NA, 
                         clustering_method = 'average',
                         show_colnames = F,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),
                                                                alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45')

rna.counts.clean %>% rownames_to_column(var = 'Symbol') %>% rio::export( file = '25.5.20/rna.clean.counts.xlsx')

# 04.Euclidean Distance ----
d_p_pca <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
}
rna.tpm.clean <- rna.tpm.pcg[pick.genes,]
rna.tpm.clean <- rna.tpm.clean[,colnames(rna.counts.clean)]
tpm.pca <- d_p_pca(rna.tpm.clean,colData$group_list)

dist.pca <- tpm.pca$ind$coord[,1:2] %>% as.data.frame()
dist.pca.euc <- dist(dist.pca,method = 'euclidean') %>% as.matrix()
colnames(dist.pca.euc)
dist.young.pbs_age.pbs <- dist.pca.euc[c(1:6,19:24), c(1:6,19:24)] %>% as.data.frame()
dist.young.pbs_age.nr <- dist.pca.euc[c(1:6,13:18), c(1:6,13:18)] %>% as.data.frame()
dist.young.pbs_age.cd38i <- dist.pca.euc[1:12, 1:12] %>% as.data.frame()

dist.young.pbs_age.pbs.long <- as.numeric(dist.young.pbs_age.pbs[upper.tri(dist.young.pbs_age.pbs)])
dist.young.pbs_age.nr.long <- as.numeric(dist.young.pbs_age.nr[upper.tri(dist.young.pbs_age.nr)])
dist.young.pbs_age.cd38i.long <- as.numeric(dist.young.pbs_age.cd38i[upper.tri(dist.young.pbs_age.cd38i)])

dist.group <- data.frame(dist = c(dist.young.pbs_age.pbs.long, dist.young.pbs_age.nr.long, dist.young.pbs_age.cd38i.long), group = c(rep('Young_Aged',66),rep('Young_Aged+NR',66), rep('Young_Age+78c',66)))

my_comparisons <- list(c("Young_Aged", "Young_Aged+NR"),c('Young_Aged','Young_Age+78c'))
ggpubr::ggboxplot(dist.group, x="group", y="dist", width = 0.4, 
                       # add = 'boxplot',add.params = list(fill = 'gray90',alpha = 0.6),
                       color = "black",#轮廓颜色 
                       fill="group",#填充
                       # palette = "npg",
                       xlab = F, #不显示x轴的标签
                       bxp.errorbar=T,#显示误差条
                       bxp.errorbar.width=0.2, #误差条大小
                       size=0.5, #箱型图边线的粗细 
                       outlier.shape=NA, #不显示outlier
                       legend = "right",
                       alpha = 0.6) + 
  geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('Euclidean distance') +
  scale_fill_manual(values = c('#d62828','#eae2b7','#00afb9')) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),legend.position = 'none') +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")


# 05.DEGs ----

options(scipen = 10)

library(IHW)
library(ashr)
library(tidyverse)

# young vs age
young_age <- rna.counts.clean[,c(1:6,19:24)]
colData_young_age <- data.frame(group_list = colData[c(1:6,19:24),])
rownames(colData_young_age) <- colnames(young_age)

dds <- DESeqDataSetFromMatrix(countData = young_age,colData = colData_young_age,design = ~group_list)
dds <- DESeq(dds,minReplicatesForReplace = 3)
plotDispEsts(dds)
resultsNames(dds) 

young_age_ashr <- lfcShrink(dds,type = 'youyihmezhzhuzhukan',
                            coef='group_list_Aged_vs_Young') # adjust LFC
young_age_ashr_res <- young_age_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
young_age_ihw <- results(dds, contrast = c('group_list','Aged','Young'), filterFun = ihw) # adjust P value
young_age_ihw_res <- young_age_ihw %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
young_age_ihw_ashr_res <- young_age_ihw_res %>%
  left_join(y = young_age_ashr_res[,c(1,3)], by = 'geneid')
young_age_ihw_ashr_res$change <- factor(ifelse(young_age_ihw_ashr_res$padj < 0.05 & young_age_ihw_ashr_res$log2FoldChange.y > log2(1.5),'Aged up', ifelse(young_age_ihw_ashr_res$padj < 0.05 & young_age_ihw_ashr_res$log2FoldChange.y < -log2(1.5) ,'Young up','No sig')))
table(young_age_ihw_ashr_res$change)

pheatmap::pheatmap(
  na.omit(rna.tpm.clean[NAD.dehy.comp.gene[[1]],]),
  scale = 'row',
  show_rownames = F,
  cluster_cols = F
)
pheatmap::pheatmap(
  na.omit(rna.tpm.clean[sasp.geneuse[[1]],]),
  scale = 'row',
  show_rownames = F,
  cluster_cols = F
)
pheatmap::pheatmap(
  na.omit(rna.tpm.clean[nad.geneuse[[1]],c(1:6,19:24)]),
  scale = 'row',
  show_rownames = F,
  cluster_cols = F
)
pheatmap::pheatmap(
  na.omit(rna.tpm.clean[infm.geneuse[[1]],]),
  scale = 'row',
  clustering_method = 'ward.D',
  show_rownames = F,
  cluster_cols = F
)


# nr vs age
nr_age <- rna.counts.clean[,c(13:24)]
colnames(nr_age)
colData_nr_age <- data.frame(group_list = colData[c(13:24),])
rownames(colData_nr_age) <- colnames(nr_age)

dds <- DESeqDataSetFromMatrix(countData = nr_age,colData = colData_nr_age,design = ~group_list)
dds <- DESeq(dds,minReplicatesForReplace = 3)
plotDispEsts(dds)
resultsNames(dds) 

nr_age_ashr <- lfcShrink(dds,type = 'apeglm',
                         coef='group_list_Aged.NR_vs_Aged') # adjust LFC
nr_age_ashr_res <- nr_age_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
nr_age_ihw <- results(dds, contrast = c('group_list','Aged','Aged+NR'), filterFun = ihw) # adjust P value
nr_age_ihw_res <- nr_age_ihw %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
nr_age_ihw_ashr_res <- nr_age_ihw_res %>%
  left_join(y = nr_age_ashr_res[,c(1,3)], by = 'geneid')
nr_age_ihw_ashr_res$log2FoldChange.y <- -nr_age_ihw_ashr_res$log2FoldChange.y

nr_age['Eif1ad12',]

nr_age_ihw_ashr_res$change <- factor(ifelse(nr_age_ihw_ashr_res$padj < 0.05 & nr_age_ihw_ashr_res$log2FoldChange.y > log2(1.5),'Aged up', ifelse(nr_age_ihw_ashr_res$padj < 0.05 & nr_age_ihw_ashr_res$log2FoldChange.y < -log2(1.5) ,'NR up','No sig')))
table(nr_age_ihw_ashr_res$change)

# nr vs age
c78_age <- rna.counts.clean[,c(7:12,19:24)]
colnames(c78_age)
colData_c78_age <- data.frame(group_list = colData[c(7:12,19:24),])
rownames(colData_c78_age) <- colnames(c78_age)

dds <- DESeqDataSetFromMatrix(countData = c78_age,colData = colData_c78_age,design = ~group_list)
dds <- DESeq(dds,minReplicatesForReplace = 3)
plotDispEsts(dds)
resultsNames(dds) 

c78_age_ashr <- lfcShrink(dds,type = 'apeglm',
                          coef='group_list_Aged.78c_vs_Aged') # adjust LFC
c78_age_ashr_res <- c78_age_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
c78_age_ihw <- results(dds, contrast = c('group_list','Aged','Aged+78c'), filterFun = ihw) # adjust P value
c78_age_ihw_res <- c78_age_ihw %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')
c78_age_ihw_ashr_res <- c78_age_ihw_res %>%
  left_join(y = c78_age_ashr_res[,c(1,3)], by = 'geneid')
c78_age_ihw_ashr_res$log2FoldChange.y <- -c78_age_ihw_ashr_res$log2FoldChange.y
c78_age['Wfdc10',]

c78_age_ihw_ashr_res$change <- factor(ifelse(c78_age_ihw_ashr_res$padj < 0.05 & c78_age_ihw_ashr_res$log2FoldChange.y > log2(1.5),'Aged up', ifelse(c78_age_ihw_ashr_res$padj < 0.05 & c78_age_ihw_ashr_res$log2FoldChange.y < -log2(1.5) ,'C78 up','No sig')))
table(c78_age_ihw_ashr_res$change)

rnaseq.deg <- list(young_age = young_age_ihw_ashr_res, nr_age = nr_age_ihw_ashr_res, c78_age = c78_age_ihw_ashr_res)
rio::export(rnaseq.deg, file = './25.5.20/rnaseq.deg.xlsx')


# ciber
ciber.fc <- read.csv('decon/ciber_function0523.csv')
ciber.fc$Mixture <- pick.samp
ciber.fc <- ciber.fc[,-c(17:19)] %>% column_to_rownames(var = 'Mixture') %>% scale() %>%  t() %>% as.data.frame() %>% na.omit()
# ciber.fc <- ciber.fc %>% t() %>% as.data.frame()

split= factor(c(rep('Young',6),rep('Aged+78c',6),rep('Aged+NR',6),rep('Aged',6)), levels = c('Young','Aged+78c','Aged+NR','Aged'))
ciber.fc <- ciber.fc[,c(7:12,1:6,13:24)]

ComplexHeatmap::pheatmap(ciber.fc,
                         name = "Spearman's coeff",
                         # scale = 'row',
                         border_color = NA, 
                         cluster_cols = F,cluster_rows = F,
                         column_split = split,
                         show_colnames = T,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),
                                                                alpha=T,bias=1)(256),alpha = 1),
                         top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c('#00afb9','#d62828','#eae2b7','#f77f00')),labels = c("Young","Aged","78C","NR"),labels_gp = gpar(col = "white"))),
                         angle_col = '45')

ciber.ic <- read.csv('decon/ciber_immune0523.csv')
ciber.ic$Mixture <- pick.samp
ciber.ic <- ciber.ic %>% column_to_rownames(var = 'Mixture') %>% scale() %>%  t() %>% as.data.frame() %>% na.omit()
split= c(rep(1,6),rep(2,6),rep(3,6),rep(4,6))
ciber.ic <- ciber.ic[,c(7:12,1:6,13:24)]

ComplexHeatmap::pheatmap(ciber.ic,
                         name = "Spearman's coeff",
                         scale = 'row',
                         border_color = NA, 
                         cluster_cols = F,cluster_rows = F,
                         column_split = split,
                         show_colnames = T,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),
                                                                alpha=T,bias=1)(256),alpha = 1),
                         top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c('#00afb9','#d62828','#eae2b7','#f77f00')),labels = c("Young","Aged","78C","NR"),labels_gp = gpar(col = "white"))),
                         angle_col = '45')

ciber.all <- rbind(ciber.fc,ciber.ic)
ComplexHeatmap::pheatmap(ciber.all,
                         name = " ",
                         scale = 'row',
                         border_color = NA, 
                         cluster_cols = F,cluster_rows = F,
                         column_split = split,
                         show_colnames = F,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),
                                                                alpha=T,bias=1)(256),alpha = 1),
                         top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c('#00afb9','#d62828','#eae2b7','#f77f00')),labels = c("Young","Aged","78C","NR"),labels_gp = gpar(col = "white"))),
                         
                         angle_col = '45')

## hl pick pathway
sc.pick <- rio::import('25.5.20/hd.pick.path.xlsx')

sc.pick$group = rep(c("Young","Aged"), each = 8)
colnames(sc.pick)
go_up$LogP <- -go_up$LogP
go_up <- sc.pick %>%
  group_by(Description) %>%
  # 如果同一通路在同一组里可能有多条记录，可取最大值
  # summarize(LogP = max(LogP), group = first(group)) %>%
  ungroup() %>%
  # 重新设因子水平，最小 LogP 在上，最大在下
  mutate(Pathway = fct_reorder(Description, LogP))
go_up$Description <- factor(go_up$Description, levels = rev(go_up$Description))

ggplot(go_up, aes(x = LogP, y = Description, fill = group)) +
  geom_col(width = 0.7) +
  # 按组做左右两列分面；free_x 保证每列横轴独立缩放
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  # 美化横轴
  scale_x_continuous(
    labels = number_format(accuracy = 1),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  # 自定义配色（可改）
  scale_fill_manual(
    values = c("Aged" = "#d62828", "Young" = "#00afb9")
  ) +
  labs(
    x    = expression(-log[10](Pvalue)),
    y    = NULL,
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.text         = element_text(face = "bold"),
    axis.text.y        = element_text(size = 10)
  )

sc.pick <- rio::import('25.5.20/sc.pick.path.xlsx')
sc.pick$LogP <- -sc.pick$LogP
sc.pick$group <- rep(c("Young","Aged"), each = 10)
sc.pick.young <- sc.pick[sc.pick$group == 'Young',]
sc.pick.young <- sc.pick.young[order(sc.pick.young$LogP),]
sc.pick.age <- sc.pick[sc.pick$group == 'Aged',]
sc.pick.age <- sc.pick.age[order(sc.pick.age$LogP),]
sc.pick <- rbind(sc.pick.young,sc.pick.age)
# sc.pick$LogP <- -sc.pick$LogP
sc.pick$Description <- factor(sc.pick$Description, levels = c(sc.pick$Description))

ggplot(sc.pick, aes(x = LogP, y = Description, fill = group)) +
  geom_col(width = 0.7) +
  # 按组做左右两列分面；free_x 保证每列横轴独立缩放
  # facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  # 美化横轴
  scale_x_continuous(
    labels = number_format(accuracy = 1),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  # 自定义配色（可改）
  scale_fill_manual(
    values = c("Aged" = "#d62828", "Young" = "#00afb9")
  ) +
  labs(
    x    = expression(-log[10](Pvalue)),
    y    = NULL,
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.text         = element_text(face = "bold"),
    axis.text.y        = element_text(size = 10)
  )


# 10.pick pathway heatmap ----

# raw
ha_left.col <- list(group_list = c('#00afb9','#eae2b7','#f77f00','#d62828'))
names(ha_left.col$group_list) <- c('Young','Aged+78c','Aged+NR','Aged')

split= factor(c(rep('Young',6),rep('Aged+78c',6),rep('Aged+NR',6),rep('Aged',6)), levels = c('Young','Aged','Aged+78c','Aged+NR'))
ComplexHeatmap::pheatmap(na.omit(rlog_mat[sasp.geneuse[[1]],]),
                         scale = 'row',
                         name = "SASP",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         cluster_cols = F,
                         show_colnames = F,
                         show_rownames = F,
                         column_split = split,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),
                                                                alpha=T,bias=1)(256),alpha = 1),
                         annotation_col = colData,
                         annotation_colors = ha_left.col,
                         angle_col = '45')

ComplexHeatmap::pheatmap(na.omit(rlog_mat[infm.geneuse[[1]],]),
                         scale = 'row',
                         name = "INFM",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         cluster_cols = F,
                         show_colnames = F,
                         show_rownames = F,
                         column_split = split,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),
                                                                alpha=T,bias=1)(256),alpha = 1),
                         annotation_col = colData,
                         annotation_colors = ha_left.col,
                         angle_col = '45')


ComplexHeatmap::pheatmap(na.omit(rlog_mat[NAD.dehy.comp.gene[[1]],]),
                         scale = 'row',
                         name = "NADH",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         cluster_cols = F,
                         show_colnames = F,
                         show_rownames = F,
                         column_split = split,
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),
                                                                alpha=T,bias=1)(256),alpha = 1),
                         annotation_col = colData,
                         annotation_colors = ha_left.col,
                         angle_col = '45')

sasp.gene.inset <- intersect(sasp.geneuse[[1]],rownames(rlog_mat)) %>% unlist()
infm.gene.inset <- intersect(infm.geneuse[[1]], rownames(rlog_mat)) %>% unlist()
nadh.gene.inset <- intersect(NAD.dehy.comp.gene[[1]], rownames(rlog_mat)) %>% unlist()
gene_group <- factor(c(rep('SASP',length(sasp.gene.inset)),rep('INFM',length(infm.gene.inset)),rep('NADH',length(nadh.gene.inset))),levels = c('SASP','INFM','NADH'))

pick.rlog <- rlog_mat[c(sasp.gene.inset,infm.gene.inset,nadh.gene.inset), ]

Heatmap(
  t(scale(t(pick.rlog))),
  name              = "Expr",           
  col               = colorRamp2(c(-3,0,3), c("navy","white","firebrick")),
  cluster_rows      = TRUE,          
  cluster_columns   = FALSE,          
  row_split         = gene_group,       
  column_split = split,
  show_row_names    = FALSE,
  show_column_names = F,
  # row_title         = "Gene group",
  # row_title_side    = "left",
  show_row_dend = F,
  layer_fun = function(nr, nc, r_index, c_index, ...) {
    grid.rect(
      x      = unit(0.5, "npc"),
      y      = unit(0.5, "npc"),
      width  = unit(1, "npc"),
      height = unit(1, "npc"),
      gp     = gpar(col = "black", fill = NA, lwd = 1)
    )
  }
  # heatmap_legend_param = list(title_position="top")
)
