
############################################################################## scRNA-seq section ########################################################################################

# umap plot ----

## raw umap
Clean_sct.inte$cell.type.percise <- factor(Clean_sct.inte$cell.type.percise, levels = c('GCs','Syns','S-TGCs','SpTs','LpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells'))
p1 <- DimPlot(Clean_sct.inte, group.by = 'cell.type.percise', cols = c( npg.cols,use.cols), label = F) + ggtitle("") + NoAxes()
p2 <- DimPlot(Clean_sct.inte, group.by = 'cell.type.percise', cols = c(npg.cols,use.cols), label = T) + ggtitle("") + NoAxes() + NoLegend()
p3 <- DimPlot(Clean_sct.inte, group.by = 'cell.type.percise', cols = c(npg.cols,use.cols), label = T,repel = T) + ggtitle("") + NoAxes()
p1 + p2 + p3

Clean_sct.inte$cell.type.percise.rm.syn <- Clean_sct.inte$cell.type.percise
Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.percise == 'Syns',]$cell.type.percise.rm.syn <- 'GCs'
p1 <- DimPlot(Clean_sct.inte, group.by = 'cell.type.percise.rm.syn', cols = c(npg.cols,use.cols), label = F) + ggtitle("") + NoAxes()
p2 <- DimPlot(Clean_sct.inte, group.by = 'cell.type.percise.rm.syn', cols = c(npg.cols,use.cols), label = T) + ggtitle("") + NoAxes() + NoLegend()
p3 <- DimPlot(Clean_sct.inte, group.by = 'cell.type.percise.rm.syn', cols = c(npg.cols,use.cols), label = T,repel = T) + ggtitle("") + NoAxes()
p1 + p2 + p3

## remove LpTs
p1 <- DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c( use.cols,npg.cols), label = F) + ggtitle("") + NoAxes()
p2 <- DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols), label = T,repel = T) + ggtitle("") + NoAxes() + NoLegend()
p3 <- DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols), label = T,repel = T) + ggtitle("") + NoAxes()
p1 +p3
DimPlot(Clean_sct.inte.rm.lpt, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols), label = T,repel = T, split.by = 'condition') + ggtitle("") + NoAxes()

table(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Young',]$cell.type.percise.new)/nrow(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Young',])
table(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Age',]$cell.type.percise.new)/nrow(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Age',])


# miloR:composition of neighbourhood of cells ----
library(miloR)
library(scales)
library(SingleCellExperiment)
library(ggbeeswarm)

milo.all.seu <- subset(Clean_sct.inte.rm.lpt, cells = c(sample( 
  rownames(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Age',]), 21867),
  rownames(Clean_sct.inte.rm.lpt@meta.data[Clean_sct.inte.rm.lpt$condition == 'Young',] )))
milo.all.seu <- as.SingleCellExperiment(milo.all.seu,assay = 'SCT') %>%
  Milo() %>%
  miloR::buildGraph(k = 30, d = 50) %>% 
  makeNhoods(
             prop = 0.2,                
             k = 30,        
             d=50,                   
             refined = T)

milo.all.seu <- countCells(milo.all.seu, 
                           meta.data = data.frame(colData(milo.all.seu)),                     
                          sample="orig.ident")
milo.traj_design <- data.frame(colData(milo.all.seu))[,c("orig.ident", "condition")]
milo.traj_design$orig.ident <- as.factor(milo.traj_design$orig.ident)
milo.traj_design <- distinct(milo.traj_design)
rownames(milo.traj_design) <- milo.traj_design$orig.ident

milo.all.seu <- calcNhoodDistance(milo.all.seu, d=50)
milo.da_results <- testNhoods(milo.all.seu,                 
                         design = ~ condition,            
                         design.df = milo.traj_design)
milo.all.seu <- buildNhoodGraph(milo.all.seu)
ggplot(milo.da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(milo.da_results, aes(logFC, -log10(SpatialFDR))) +  
  geom_point() + 
  geom_hline(yintercept = 1)
scater::plotReducedDim(milo.all.seu, dimred = "UMAP", colour_by="cell.type.percise.new", text_size = 3, point_size=0.1)

milo.da_results$logFC <- -milo.da_results$logFC
miloR::plotNhoodGraphDA(milo.all.seu, milo.da_results, alpha=0.1,) +  
  scale_fill_gradient2(low="#F2F2F2",         
                       mid="#F2F2F2",
                       high="#FF5E5B",                  
                       name="log2FC",                   
                       limits=c(-5,5),                  
                       oob=squish)

milo.da_results <- annotateNhoods(milo.all.seu, milo.da_results, coldata_col = "cell.type.percise.new")
milo.da_results$cell.type.percise.new <- factor(milo.da_results$cell.type.percise.new, levels = rev(c('GCs','S-TGCs','SpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','Monocytes','DCs','GMPs','Mast cells','T cells','B cells')))

plotbee <- function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL) 
{
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
           group.by, ")?")
    }
    if (is.numeric(da.res[, group.by])) {
    }
    da.res <- mutate(da.res, group_by = da.res[, group.by])
  }
  else {
    da.res <- mutate(da.res, group_by = "g1")
  }
  if (!is.factor(da.res[, "group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, 
                                               levels = unique(group_by)))
  }
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }
  beeswarm_pos <- ggplot_build(da.res %>% mutate(is_signif = ifelse(SpatialFDR < 
                                                                      alpha, 1, 0)) %>% arrange(group_by) %>% ggplot(aes(group_by, 
                                                                                                                         logFC)) + geom_quasirandom())
  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y
  n_groups <- unique(da.res$group_by) %>% length()
  da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 
                                       1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 
                                                                                1, logFC, NA)) %>% arrange(group_by) %>% mutate(Nhood = factor(Nhood, 
                                                                                                                                               levels = unique(Nhood))) %>% mutate(pos_x = pos_x, pos_y = pos_y) %>% 
    ggplot(aes(pos_x, pos_y, color = logFC_color)) + scale_color_gradient2() + 
    guides(color = "none") + xlab(group.by) + ylab("Log Fold Change") + 
    scale_x_continuous(breaks = seq(1, n_groups), labels = setNames(levels(da.res$group_by), 
                                                                    seq(1, n_groups))) + geom_point(size = .5) + coord_flip() + 
    theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0))
}

plotbee (milo.da_results, alpha = 0.1,
               group.by = "cell.type.percise.new") +  
  scale_color_gradient2(low="#4cc9f0",           
                        mid="#F2F2F2",                  
                        high="#FF5E5B",                  
                        limits=c(-5,5),                  
                        oob=squish) + 
  labs(x="", y="Log2 Fold Change") +  
  theme_bw(base_size=10)+  
  theme(axis.text = element_text(colour = 'black')) + 
  scale_size_continuous(range = .1)

# ltsr composition of cell types ----
library(magrittr)
library(lme4)
library(numDeriv)
## data perpare
Clean_sct.inte.rm.lpt$orig.ident <- factor(Clean_sct.inte.rm.lpt$orig.ident, levels = c(paste0('Young',seq(1,3)), paste0('Age',seq(1,3))))
Clean_sct.inte.rm.lpt$condition <- factor(Clean_sct.inte.rm.lpt$condition, levels = c('Young','Age'))
cell.number <- FetchData(Clean_sct.inte.rm.lpt, 
                         vars = c("orig.ident", "cell.type.percise.new")) %>%
  dplyr::count(orig.ident, cell.type.percise.new) %>% 
  tidyr::spread(orig.ident, n) 

sample_ids <- colnames(cell.number)[-1]
cell_types <- cell.number$cell.type.percise.new
n_cells_per_sample <- colSums(cell.number[,-1])
n_var_cats <- 2 

sample_cats <- tibble(
  Sample_ID = sample_ids,
  Treatment = c(rep('Young',3),rep('Age',3)),
  Rep = c(rep(c('one','two','three'),2)),
  cell.num = n_cells_per_sample
)

sample_num1_values<- rep(1,6)
obs_tbl <- data.frame(
  Sample_ID = rep(sample_ids, c(n_cells_per_sample)),
  Treatment = rep(sample_cats$Treatment, c(n_cells_per_sample)),
  Rep = rep(sample_cats$Rep, c(n_cells_per_sample)),
  Var_Num1 = rep(sample_num1_values, c(n_cells_per_sample))
)

obs_tbl$Cell_type <- c(rep(cell.number$cell.type.percise.new,c(cell.number$Young1)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$Young2)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$Young3)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$Age1)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$Age2)),
                       rep(cell.number$cell.type.percise.new,c(cell.number$Age3))
)

## RUN LTSR 
source('LTSR.raw.code.R')
results <- CellTypeCompositionAnalysis(obs_tbl, "Sample_ID", "Cell_type", c("Treatment",'Rep'), "Var_Num1")
ranef_tbl <- results$ranef
sdse_tbl <- results$sdse

vars1 <- list(Treatment = c('Young','Age'))

plot_ranef.new <- function(ranef_tbl, vars, celltypes = NULL, celltype_order = "hclust", references = NULL,
         maxFC = 3, LTSR2p = F, highlightLtsr = 0.0, filterLtsr = 0.0, swap_axes = F) {
  ranef_tbl <- .getCondValLtsr(ranef_tbl, vars, celltypes = celltypes, references = references)
  save(ranef_tbl, file = "ranef_tbl.RData")
  condval_mat <- ranef_tbl %>%
    dplyr::select(
      Celltype, grpval, condval
    ) %>%
    spread(
      "grpval", "condval"
    ) %>%
    column_to_rownames(
      var = "Celltype"
    ) %>%
    as.matrix()
  if (length(celltype_order) == 1 && celltype_order == "hclust") {
    dendy <- hclust(dist(condval_mat))
    ordered_celltype <- rownames(condval_mat)[dendy$ord]
  } else if (!is.null(celltype_order) && length(celltype_order) == dim(condval_mat)[1]) {
    ordered_celltype <- celltype_order
  }
  
  ranef_tbl <- ranef_tbl %>% mutate(
    Celltype = factor(Celltype, levels = ordered_celltype),
    condval = condval %>% pmin(log(maxFC)) %>% pmax(log(1 / maxFC)),
    ltsr = ltsr %>% pmin(0.9999) %>% pmax(0.5)
  )
  
  if (swap_axes) {
    ranef_tbl$Celltype <- factor(ranef_tbl$Celltype, levels = rev(levels(ranef_tbl$Celltype)))
    ranef_tbl$grpval <- factor(ranef_tbl$grpval, levels = rev(levels(ranef_tbl$grpval)))
  }
  
  if (filterLtsr > 0) {
    filtered_celltypes <- ranef_tbl %>%
      group_by(Celltype) %>%
      summarise(maxLtsr = max(ltsr)) %>%
      dplyr::filter(maxLtsr >= filterLtsr) %>%
      dplyr::select(Celltype) %>%
      unlist(use.names = F)
    ranef_tbl <- ranef_tbl %>% dplyr::filter(Celltype %in% filtered_celltypes)
  }
  
  geom_dots <- geom_point(
    aes(
      fill = log2(exp(condval)),
      size = -log10(1 - ltsr)
    ),
    color = "white",
    shape = 21
  )
  
  if (swap_axes) {
    p <- (
      ggplot(ranef_tbl, aes(y = grpval, x = Celltype)) +
        facet_grid(grpvar ~ ., scales = "free_y", space = "free_y", switch = "x") +
        geom_dots
    )
  } else {
    p <- (
      ggplot(ranef_tbl, aes(x = grpval, y = Celltype)) +
        facet_grid(. ~ grpvar, scales = "free_x", space = "free_x", switch = "x") +
        geom_dots
    )
  }
  
  p <- (
    p + scale_fill_distiller(
      palette = "RdBu",
      limits = log2(c(1 / maxFC, maxFC)),
      breaks = log2(c(1 / maxFC, maxFC)),
      labels = c(paste0("1/", maxFC), maxFC),
      oob = squish,
      guide = guide_colorbar(
        title = "Fold change", title.position = "top", direction = "horizontal",
        barwidth = 5, barheight = 0.75, raster = F, order = 1)
    )
    + scale_size(
      limits = -log10(1 - c(0.5, 0.9999)),
      breaks = -log10(1 - c(0.5, 0.7, 0.9, 0.99)),
      range = c(0.5, 9),
      labels = ifelse(
        rep(LTSR2p, 4),
        c("0.5", "0.3", "0.1", "<0.01"),
        c("0.5", "0.7", "0.9", ">0.99")
      ),
      guide = guide_legend(
        title = ifelse(LTSR2p, "p", "LTSR"), reverse = T, order = 2,
        override.aes = list(fill = "black", color = "white")
      )
    )
    + theme_bw()
    + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      legend.spacing.y = unit(0.5, "line")
    )
  )
  
  if (highlightLtsr > 0) {
    p <- (
      p + geom_point(
        aes(
          color = ltsr > highlightLtsr,
          alpha = ltsr > highlightLtsr
        ),
        shape = 21, size = 9
      )
      + scale_color_manual(
        label = c(
          "", ifelse(LTSR2p, paste0("p < ", 1 - highlightLtsr), paste0("LTSR > ", highlightLtsr))
        ),
        values = c("white", "red"),
        guide = guide_legend(title = NULL, override.aes = list(size = 9), reverse = T, order = 3)
      )
      + scale_alpha_manual(values = c(0, 1), guide = F)
    )
  }
  
  p
}


ranef_plot <- plot_ranef.new(ranef_tbl, vars = vars1, celltypes = cell_types, celltype_order = rev(cell_types),
                         maxFC = 1.5, LTSR2p = FALSE) + xlab('Condition')
sdse_plot <- plot_sdse(sdse_tbl, "Sample_ID", ci = 0.95, xlim = c(0, 1))

ranef_plot + sdse_plot

## sankey
library(ggalluvial)

plot.data.num <- FetchData(Clean_sct.inte.rm.lpt, 
                           vars = c("condition", "cell.type.percise.new")) %>%
  dplyr::count(condition, cell.type.percise.new) %>% 
  tidyr::spread(condition, n) 

plot.data.num[is.na(plot.data.num)] <- 0

plot.data <- FetchData(Clean_sct.inte.rm.lpt, 
                       vars = c("condition", "cell.type.percise.new")) %>%
  dplyr::count(condition, cell.type.percise.new) %>% 
  group_by(condition) %>% 
  summarise(Prop = n / sum(n))

plot.data$large_ct <- factor(rep(c(levels(Clean_sct.inte.rm.lpt$cell.type.percise.new)), 2), levels = levels(Clean_sct.inte.rm.lpt$cell.type.percise.new))

plot.data$levels <- c(rep(1:15,2))
plot.data$Prop <- plot.data$Prop*100
ggplot(plot.data,
       aes(x = condition, stratum = large_ct, alluvium = levels,
           y = Prop,
           fill = large_ct, label = large_ct)) +
  # scale_fill_manual(values = c(lancent.col[c(2,1,3:6)])) +
  scale_fill_manual(values = c(use.cols,npg.cols)) +  
  geom_flow(stat = "alluvium", 
            lode.guidance = "frontback",
            width = 0.3,
            color = "darkgray") +
  geom_stratum(alpha = .8) +
  guides(fill = guide_legend(title = 'Cell type')) +
  ylab('Proportion (%)')+
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15))

# DEGs in cell types ----
library(ClusterGVis)

ht.data <- prepareDataFromscRNA(Clean_sct.inte.rm.lpt,
                                diffData = rm.lpt.ct.deg,
                                showAverage = T,
                                assays = 'SCT',slot = 'data',
                                group.by = 'cell.type.percise.new',keep.uniqGene = F,
                                scale.data = T)

enrich.go <- enrichCluster(object = ht.data,
                           OrgDb = org.Mm.eg.db,
                           type = "BP",
                           organism = "mmu",
                           pvalueCutoff = 0.01,
                           topn = 5,
                           seed = 1234)
enrich.go$ratio <- -log(enrich.go$pvalue)
rio::export(rm.lpt.ct.deg, file = 'R.hd/table1.scRNAseq cell types DEGs.xlsx')

pdf('./add.analysis/rm.lpts/24.20.newplot/05.cell.type.degs.heatmap.pdf',height = 12,width = 14,onefile = F)

visCluster(object = ht.data,
           ht.col.list = list(col_range = c(-4, 0, 4)),
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           # markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich.go,
           line.side = "left",
           cluster.order = c(1:15),
           go.col = rep(c(use.cols,npg.cols[1:2]),each = 5),
           sample.col = c(use.cols,npg.cols[1:2]),
           ctAnno.col = c(use.cols,npg.cols[1:2]),
           go.size = 8,
           add.bar = T)
dev.off()
# dotplot of key genes -----
library(plot1cell)

complex_dotplot_single <- function (seu_obj, feature, celltypes = NULL, groups, splitby = NULL, 
                                    color.palette = NULL, font.size = 12, strip.color = NULL, 
                                    do.scale = T, scale.by = "radius") 
{
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1", 
                                        "indianred1", "darkred"))(255)
  }
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (is.null(celltypes)) {
    celltypes <- levels(seu_obj)
  }
  if (length(groups) == 1) {
    groups_level <- levels(seu_obj@meta.data[, groups])
    if (is.null(groups_level)) {
      seu_obj@meta.data[, groups] <- factor(seu_obj@meta.data[, 
                                                              groups], levels = names(table(seu_obj@meta.data[, 
                                                                                                              groups])))
      groups_level <- levels(seu_obj@meta.data[, groups])
    }
    if (!is.null(splitby)) {
      if (is.null(levels(seu_obj@meta.data[, splitby]))) {
        seu_obj@meta.data[, splitby] <- factor(seu_obj@meta.data[, 
                                                                 splitby], levels = names(table(seu_obj@meta.data[, 
                                                                                                                  splitby])))
      }
      splitby_level <- levels(seu_obj@meta.data[, splitby])
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = c(groups, 
                                                                             splitby))
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], count_df[, splitby], 
                                  sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$splitby <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                            split = "___"), FUN = function(x) {
                                                              x[[3]]
                                                            }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$splitby <- factor(data_plot$splitby, levels = splitby_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    else {
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = groups)
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
      geom_tile(fill = "white", color = "white") + geom_point(aes(colour = avg.exp, 
                                                                  size = pct.exp)) + scale_color_gradientn(colours = color.palette) + 
      theme(panel.background = element_rect(fill = "white", 
                                            colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                          hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                            2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
            legend.text = element_text(size = (font.size - 
                                                 2)), legend.title = element_text(size = (font.size)), 
            strip.text = element_text(size = font.size), 
            legend.position = "right") + ylab("") + xlab("") + 
      ggtitle(feature)
    if (do.scale) {
      p = p + scale_size(range = c(0, 10))
    }
    else {
      if (max(data_plot$pct.exp) >= 20) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                        20))
      }
    }
    if (!is.null(splitby)) {
      p <- p + facet_wrap(~splitby, scales = "free_x")
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid.draw(g))
    }
    else {
      p
    }
  }
  else {
    gene_count <- extract_gene_count(seu_obj = seu_obj, 
                                     features = feature, cell.types = celltypes, meta.groups = c(groups, 
                                                                                                 splitby))
    allgroups <- c(groups, splitby)
    for (i in 1:length(allgroups)) {
      if (is.null(levels(seu_obj@meta.data[, allgroups[i]]))) {
        seu_obj@meta.data[, allgroups[i]] <- factor(seu_obj@meta.data[, 
                                                                      allgroups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                            allgroups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, allgroups[i]])
      gene_count[, allgroups[i]] <- factor(gene_count[, 
                                                      allgroups[i]], levels = group_level)
    }
    gene_count$celltype <- factor(gene_count$celltype, levels = celltypes)
    all_levels <- list()
    for (i in 1:length(groups)) {
      if (is.null(levels(seu_obj@meta.data[, groups[i]]))) {
        seu_obj@meta.data[, groups[i]] <- factor(seu_obj@meta.data[, 
                                                                   groups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                      groups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, groups[i]])
      all_levels[[i]] <- group_level
    }
    all_levels <- as.character(unlist(all_levels))
    data_plot <- list()
    for (i in 1:length(groups)) {
      count_df <- gene_count
      count_df$new_group <- paste(gene_count[, groups[i]], 
                                  gene_count[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      df1 <- merge(exp_df, pct_df, by = "new_group")
      df1$groupID <- groups[i]
      data_plot[[i]] <- df1
    }
    data_plot <- do.call("rbind", data_plot)
    data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                         split = "___"), FUN = function(x) {
                                                           x[[1]]
                                                         }))
    data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[2]]
                                                           }))
    data_plot$groups <- factor(data_plot$groups, levels = all_levels)
    data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    data_plot$groupID <- factor(data_plot$groupID, levels = groups)
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    if (is.null(splitby)) {
      p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_wrap(~groupID, scales = "free_x")
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid::grid.draw(g))
    }
    else {
      df2 <- reshape2::melt(gene_count[, c(groups, splitby)], 
                            measure.vars = groups)
      df2 <- df2[!duplicated(df2$value), ]
      colnames(df2)[colnames(df2) == "value"] <- "groups"
      data_plot2 <- list()
      for (i in 1:length(groups)) {
        df3 <- data_plot[data_plot$groupID == groups[i], 
        ]
        df4 <- df2[df2$variable == groups[i], c("groups", 
                                                splitby[i])]
        colnames(df4)[2] <- "split"
        df5 <- merge(df3, df4, by = "groups")
        data_plot2[[i]] <- df5
      }
      data_plot2 <- do.call("rbind", data_plot2)
      fill_x1 <- grDevices::rainbow(length(groups), alpha = 0.5)
      fill_x2 <- list()
      for (i in 1:length(splitby)) {
        n_col <- unique(gene_count[, splitby[i]])
        fill_x2[[i]] <- (scales::hue_pal(l = 90))(length(n_col))
      }
      fill_x2 <- as.character(unlist(fill_x2))
      fill_x <- c(fill_x1, fill_x2)
      p <- ggplot(data_plot2, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_nested(~groupID + split, 
                                        scales = "free_x", strip = strip_nested(background_x = elem_list_rect(fill = fill_x)))
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      p
    }
  }
}
com.dot.new <- com.dot.new <- function (seu_obj, features, celltypes = NULL, groups, color.palette = NULL, 
                                        strip.color = NULL) 
{
  pb <- progress_bar$new(format = "  Ploting [:bar] :percent eta: :eta", 
                         clear = FALSE, total = length(features), width = 100)
  plot_list <- list()
  for (i in 1:length(features)) {
    pp <- invisible(complex_dotplot_single(seu_obj = seu_obj, 
                                           feature = features[i], groups = groups, celltypes = celltypes))
    pp <- pp$data
    pp$gene <- features[i]
    plot_list[[i]] <- pp
    pb$tick()
    Sys.sleep(1/length(features))
  }
  all_data <- do.call("rbind", plot_list)
  all_data$gene <- factor(all_data$gene, levels = rev(features))
  all_data$celltype <- factor(all_data$celltype, levels = levels(seu_obj))
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1",
                                        "indianred1", "darkred"))(255)
  }
  p <- invisible(ggplot(all_data, aes(x = groups, y = gene)) + 
                   geom_tile(fill = "white", color = "white") + 
                   geom_point(aes(colour = avg.exp, size = pct.exp), alpha = 0.9) + 
                   scale_color_gradientn(colours = color.palette) + 
                   scale_size(range = c(0, 5)) +
                   theme(
                     panel.background = element_rect(fill = "white", colour = "black"), 
                     axis.text.x = element_text(angle = 45,hjust = 1),
                     axis.text.y = element_text(face = 'italic'),
                     plot.title = element_text(size = 10, hjust = 0.5,face = "bold"), 
                     axis.text = element_text(size = 12), 
                     axis.title = element_text(size = 8), 
                     legend.text = element_text(size = 8), 
                     legend.title = element_text(size = 12),
                     legend.position = "right", 
                     strip.text = element_text(size = 8, colour = "black",face = "bold")) + 
                   ylab("") + xlab("") + ggtitle("") + 
                   facet_wrap(~celltype, ncol = length(levels(seu_obj))))
  g <- change_strip_background(p, type = "top", strip.color = strip.color)
  print(grid.draw(g))
}
extract_gene_count <- function (seu_obj, features, cell.types = NULL, data.type = "data", 
                                meta.groups = NULL) 
{
  if (is.null(cell.types)) {
    cell.types = levels(seu_obj)
  }
  seu_obj@meta.data$celltype <- as.character(seu_obj@active.ident)
  if (is.null(meta.groups)) {
    meta.groups = colnames(seu_obj@meta.data)
  }
  if (!is.null(cell.types)) {
    new_seu <- subset(seu_obj, idents = cell.types)
  }
  feature_count <- Seurat::FetchData(new_seu, slot = data.type, 
                                     vars = c(features, meta.groups, "celltype"))
  umap_data <- data.frame(new_seu[["umap"]]@cell.embeddings)
  feature_count$UMAP1 <- umap_data$UMAP_1
  feature_count$UMAP2 <- umap_data$UMAP_2
  feature_count
}
change_strip_background <- function (ggplt_obj, type = "top", strip.color = NULL) 
{
  g <- ggplot_gtable(ggplot_build(ggplt_obj))
  if (type == "top") {
    strip_both <- which(grepl("strip-t", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else if (type == "right") {
    strip_both <- which(grepl("strip-r", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else {
    strip_t <- which(grepl("strip-t", g$layout$name))
    strip_r <- which(grepl("strip-r", g$layout$name))
    strip_both <- c(strip_t, strip_r)
    fills <- strip.color
    if (is.null(fills)) {
      fills <- c((scales::hue_pal(l = 90))(length(strip_t)), 
                 (scales::hue_pal(l = 90))(length(strip_r)))
    }
  }
  k <- 1
  for (i in strip_both) {
    j <- which(grepl("rect", g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  g
}

Idents(Clean_sct.inte.rm.lpt) <- Clean_sct.inte.rm.lpt$cell.type.percise.new
Clean_sct.inte.rm.lpt <- PrepSCTFindMarkers(Clean_sct.inte.rm.lpt)
rm.lpt.ct.deg <- FindAllMarkers(Clean_sct.inte.rm.lpt, only.pos = T,logfc.threshold = 0.25, min.pct = 0.2) %>% dplyr::filter(p_val_adj < 0.05)
rm.lpt.ct.deg.top10 <- rm.lpt.ct.deg %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Clean_sct.inte.rm.lpt$cell.type.percise.new <- factor(Clean_sct.inte.rm.lpt$cell.type.percise.new, levels = c('GCs','S-TGCs','SpTs','ECs','DSCs','MSCs','EryCs','Macrophages','Neutrophils','GMPs','Monocytes','DCs','Mast cells','T cells','B cells'))
Idents(Clean_sct.inte.rm.lpt) <- Clean_sct.inte.rm.lpt$cell.type.percise.new
com.dot.new(Clean_sct.inte.rm.lpt, feature = c('Kcnk5','Ncam1','Prl2a1',
                                               'Lepr','Ctsq','Stra6',
                                               'Prl8a9','Prl3a1','Ctsm',
                                               'Pecam1','Ptprb','Kdr',
                                               'Pbx1','Pgr','Cryab',
                                               'Col1a1','Mfap4','Myl9',
                                               'Hba-a1','Hba-a2','Hbb-bs',
                                               'C1qc','Apoe','Ms4a7',
                                               'S100a9','G0s2','S100a8',
                                               'Ngp','Stfa1','Stfa2',
                                               'Plcb1','Myo1g','Adgre4',
                                               'H2-Aa','H2-Eb1','Cd74',
                                               'Mcpt8','Cpa3','Fcer1a',
                                               'Cd3d','Cd3g','Cd3e',
                                               'Cd79a','Cd79b','Ms4a1')
            ,groups = "condition",strip.color = c(use.cols,npg.cols)[1:15])

# umap of key genes ----
p1 <- FeaturePlot(Clean_sct.inte.rm.lpt, features = c('Cd38'),split.by = 'condition',order = T ,combine = F, cols = c('gray90','red3'))
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 2)

macrop.seu <- subset(Clean_sct.inte.rm.lpt, subset = (cell.type.percise.new == 'Macrophages'))
VlnPlot(macrop.seu,features =  c('Irf5','Cd38'),
        # split.by = 'cell.type.percise.new',
        split.plot = F,
        group.by = 'condition',
        assay = 'SCT',layer = 'data',
        pt.size = 0,cols = c('#0B996F',"#D6570D"))

DimPlot(macrop.seu, group.by = 'condition',split.by = 'condition')
macrop.seu[['RNA']] <- split(macrop.seu[['RNA']], f = macrop.seu$sample)
macrop.seu.inte <- macrop.seu %>%
  # SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = .5) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)
p1 <- DimPlot(macrop.seu.inte, group.by = 'condition',split.by = 'condition',cols = c(npg.cols[c(2,1)])) + 
  ggtitle('Macrophages integrated') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p2 <- FeaturePlot(macrop.seu.inte, features = 'Cd38',split.by = 'condition', order = T,min.cutoff = 'q1',
                  cols =c("grey90", "red3"),) &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p1 / p2

macrop.seu.inte <- macrop.seu.inte %>% FindNeighbors( reduction = "integrated.dr", 
               dims = 1:30) %>% 
  FindClusters(resolution = .3) 

DimPlot(macrop.seu.inte, split.by = 'condition',group.by = 'SCT_snn_res.0.3',label = T) + 
  scale_color_simpsons() +
  ggtitle('Macrophages integrated') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 

macrop.seu <- RunUMAP(macrop.seu,reduction = "integrated.dr", 
                      dims = 1:30)
DimPlot(macrop.seu, group.by = 'condition',split.by = 'condition',cols = c(npg.cols[c(2,1)]))
FeaturePlot(macrop.seu, features = 'Cd38',split.by = 'condition')
macrop.seu.split <- SplitObject(macrop.seu,split.by = 'sample')
macrop.seu.merge <- merge(x = macrop.seu.split$Young1, y = c(macrop.seu.split$Young2, 
                                                             macrop.seu.split$Young3,
                                                             macrop.seu.split$Age1,
                                                             macrop.seu.split$Age2,
                                                             macrop.seu.split$Age3))
macrop.seu.merge <- macrop.seu.merge %>% 
  SCTransform(vars.to.regress = c('mitoRatio','rpRatio','G2M.Score','S.Score')) %>%
  RunPCA() %>% 
  RunUMAP(reduction = "pca", 
                      dims = 1:30)
macrop.seu.merge$condition <- factor(macrop.seu.merge$condition, levels = c('Young','Age'))
p1 <- DimPlot(macrop.seu.merge, group.by = 'condition',split.by = 'condition') + 
  scale_color_simpsons() +
  ggtitle('Macrophages only merged') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p2 <- FeaturePlot(macrop.seu.merge, features = 'Cd38',split.by = 'condition', order = T,
                  # min.cutoff = 'q1',
                  cols =c("grey90", "red3"),) &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p1 / p2

FeaturePlot(macrop.seu.merge, features = 'Cd38',
            split.by = 'SCT_snn_res.0.3', 
            order = T,min.cutoff = 'q1',
            ncol = 4,
            cols =c("grey90", "red3"))


macrop.seu.merge <- macrop.seu.merge %>% FindNeighbors( reduction = "pca", 
                                                      dims = 1:30) %>% 
  FindClusters(resolution = .3) 

DimPlot(macrop.seu.merge, split.by = 'condition',group.by = 'SCT_snn_res.0.3',label = T) + 
  scale_color_simpsons() +
  ggtitle('Macrophages merged') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 

Idents(macrop.seu.merge) <- 'SCT_snn_res.0.3'
macro.ct.deg <- FindAllMarkers(macrop.seu.merge,only.pos = T,logfc.threshold = 0.25,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
rio::export(macro.ct.deg, file = 'R.hd/macrophages.merged.cluster.deg.xlsx')

macrop.seu.merge <- macrop.seu.merge %>% 
  FindClusters(resolution = .8) 

p1 <- DimPlot(macrop.seu.merge, group.by = 'SCT_snn_res.0.8',split.by = 'condition',label = T) +
  scale_color_simpsons() +
  ggtitle('Macrophages only merged') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p2 <- FeaturePlot(macrop.seu.merge, features = 'Cd38',split.by = 'condition', order = T,
                  # min.cutoff = 'q1',
                  cols =c("grey90", "red3"),) &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p1 / p2

Idents(macrop.seu.merge) <- 'SCT_snn_res.0.8'
macrop.seu.merge <- RenameIdents(macrop.seu.merge,
                                 '0' = 'Cd38 +',
                                 '1' = 'Cd38 +',
                                 '2' = 'Cd38 +',
                                 '3' = 'Cd38 +',
                                 '4' = 'Cd38 +',
                                 '5' = 'Cd38 +'
                                 )
macrop.seu.merge$macrop.type <- as.character(Idents(macrop.seu.merge))
macrop.seu.merge$macrop.type[macrop.seu.merge$macrop.type != 'Cd38 +' ] <- 'Cd38 -'
macrop.seu.merge$macrop.type <- factor(macrop.seu.merge$macrop.type,levels = c('Cd38 +','Cd38 -'))

p1 <- DimPlot(macrop.seu.merge, group.by = 'macrop.type',label = T) +
  scale_color_simpsons() +
  ggtitle('Macrophages subtypes') &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p2 <- FeaturePlot(macrop.seu.merge, features = 'Cd38', order = T,
                  # min.cutoff = '0.2',
                  cols =c("grey90", "red3"),) &
  NoAxes() & 
  theme(panel.background=element_rect(fill='transparent', color='black'),title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.height=unit(1,"line")) 
p1 + p2

macrop.subtype.deg <- FindMarkers(macrop.seu.merge, group.by = 'macrop.type',min.pct = 0.1,logfc.threshold = 0.25,
                                  ident.1 = 'Cd38 +', ident.2 = 'Cd38 -', recorrect_umi = F) %>% 
  dplyr::filter( p_val_adj < 0.05) %>% rownames_to_column(var = 'genes')
macrop.subtype.deg$type <- factor(ifelse(macrop.subtype.deg$avg_log2FC > 0,'Cd38 + up', 'Cd38 - up'))
rio::export(macrop.subtype.deg, file = 'R.hd/macrop.subtype.deg.xlsx')

library(ensembldb)
library(org.Mm.eg.db)
hub<-AnnotationHub::AnnotationHub()
annolist <- AnnotationHub::query(hub, "Mus musculus")
ensdb110 <- hub[["AH113713"]]
plot.list <- list()
macrop.term <- list()

library(clusterProfiler)
for (ct in levels(macrop.subtype.deg$type)) {
  tmp <- macrop.subtype.deg[macrop.subtype.deg$type == ct,]
  gene2entrzid <- bitr(tmp$genes, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,' up'))
  erich.name <- paste0('macro_',ct,'_plot')
  macrop.term[[erich.name]] <- erich.go.BP
  # assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}
patchwork::wrap_plots(c(plot.list),nrow = 1)
rio::export(macrop.term, file = 'R.hd/macrop.term.xlsx')
write.csv(macrop.term$`macro_Cd38 + up_plot`, file = 'R.hd/macrop.cd38+_up.term.csv')
write.csv(macrop.term$`macro_Cd38 - up_plot`, file = 'R.hd/macrop.cd38-_up.term.csv')

# SASP ----
sasp.genenset <- rio::import('SASP.geneset.xlsx')
sasp.geneuse <- sasp.genenset$`Gene ID` %>% tolower() %>% stringr::str_to_title() 
sasp.geneuse <- list(c(sasp.geneuse))

# run AUCells

# inflammatory reponse----
infm.genenset <- rio::import('inflammatory_response.xlsx')
infm.geneuse <- infm.genenset$inflammatory.response
infm.geneuse <- list(c(infm.geneuse))

# run AUCells

# DEGs cross condition -----
comlist <- t(combn(unique(Clean_sct.inte.rm.lpt$condition), 2))
for (ct in unique(Clean_sct.inte.rm.lpt$cell.type.percise.new)) {
  print(ct)
  Idents(Clean_sct.inte.rm.lpt) <- 'cell.type.percise.new'
  sub.seu.1 <- subset(Clean_sct.inte.rm.lpt, idents = ct)
  DefaultAssay(sub.seu.1) <- 'SCT'
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    age1 <- comlist[gp,1]
    age2 <- comlist[gp,2]
    DEGs <- FindMarkers(sub.seu.1,
                        ident.1 = age1,
                        ident.2 = age2,
                        logfc.threshold = 0.2,
                        group.by = 'condition',min.pct = 0.2,
                        recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
    DEGs.name <- paste('sc.deg.pct02_sct',ct,age1,age2,sep = '_')
    assign(DEGs.name,DEGs)
  }
}

hub<-AnnotationHub::AnnotationHub()
annolist <- AnnotationHub::query(hub, "Mus musculus")
ensdb110 <- hub[["AH113713"]]

plot.list <- list()
library(clusterProfiler)
for (ct in levels(Clean_sct.inte.rm.lpt$cell.type.percise.new)) {
  tmp <- base::get(paste('sc.deg',ct,age1,age2,sep = '_')) %>% rownames_to_column(var = 'gene')
  gene2entrzid <- bitr(tmp[tmp$avg_log2FC < 0,]$gene, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,' Age up'))
  erich.name <- paste0(ct,'_plot')
  # assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}
patchwork::wrap_plots(c(plot.list),nrow = 3)

plot.list <- list()
for (ct in levels(Clean_sct.inte.rm.lpt$cell.type.percise.new)) {
  tmp <- base::get(paste('sc.deg',ct,age1,age2,sep = '_')) %>% rownames_to_column(var = 'gene')
  gene2entrzid <- bitr(tmp[tmp$avg_log2FC < 0,]$gene, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  er.plot <- dotplot(erich.go.BP,showCategory = 15) + ggtitle(paste0(ct,' Young up'))
  erich.name <- paste0(ct,'_plot')
  # assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}
patchwork::wrap_plots(c(plot.list),nrow = 3)

# pseudobulk
library(DESeq2)

dsc.seu <- subset(Clean_sct.inte.rm.lpt, idents = 'DSCs')
dsc.seu.bulk <- AggregateExpression(dsc.seu, assays = 'RNA',group.by = 'sample') %>% as.data.frame()
colData = data.frame(row.names = colnames(dsc.seu.bulk),group_list = c(rep('Age',3),rep('Young',3)))
dds <- DESeqDataSetFromMatrix(countData = dsc.seu.bulk,colData = colData,design = ~group_list)
dds <- DESeq(dds)
dsc_ihw <- DESeq2::results(dds, contrast = c('group_list','Age','Young'), filterFun = ihw) 
dds_ashr <- lfcShrink(dds, contrast = c('group_list','Age','Young'),type = 'normal') # adjust LFC
summary(dds_ashr)
dsc_deseq<- dds_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  na.omit() %>% 
  arrange(log2FoldChange)
dsc_deseq$change <- factor(ifelse(dsc_deseq$padj < 0.05 & dsc_deseq$log2FoldChange > 0.75,'Age up',
                                  ifelse(dsc_deseq$padj < 0.05 & dsc_deseq$log2FoldChange < -0.75, 'Young up','not')
                                  ))
table(dsc_deseq$change)
FeaturePlot(dsc.seu, features = 'Fcgr4',split.by = 'condition',order = T)
rio::export(dsc_deseq, file = '25.3/dsc_deseq.lfc0.75.xlsx')

# NAD pathway -----
library(KEGGREST)
nad.pathway <- keggGet('mmu00760')
nad.pathway <- sapply(seq(2,86,2), function(x){
  gns <- unlist(strsplit(nad.pathway[[1]][['GENE']][x],";"))[1]
})

nad.geneuse <- list(c(nad.pathway))


# cellchat -----
library(CellChat)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "Macrophages", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "DSCs", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netVisual_bubble(cellchat.merge, sources.use = 8, targets.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Age", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat.merge, sources.use = 8, targets.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Age", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- rankNet(cellchat.merge, mode = "comparison", measure = "weight", sources.use = 8, targets.use = 5, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.merge, mode = "comparison", measure = "weight", sources.use = 8, targets.use = 5, stacked = F, do.stat = TRUE)

gg1 + gg2


pos.dataset = "Age"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

cellchat.merge <- identifyOverExpressedGenes(cellchat.merge, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

net <- netMappingDEG(cellchat.merge, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat.merge, net = net, datasets = "Age",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Young",ligand.logFC = -0.05, receptor.logFC = NULL)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 8, targets.use = 5, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 8, targets.use = 5, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

## new 
netVisual_heatmap(Age.cc,  color.heatmap = "Reds")
netVisual_heatmap(Young.cc,  color.heatmap = "Reds")

netVisual_heatmap(Age.cc,  color.heatmap = "Reds",measure = 'weight')
netVisual_heatmap(Young.cc,  color.heatmap = "Reds",measure = 'weight')

gg1 <- netVisual_heatmap(cellchat.merge,)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat.merge, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
