plot_heatmap_fixed_range <- function (
    object,
    comparison = c(1, 2),
    measure = c("count", "weight"),
    signaling = NULL,
    slot.name = c("netP", "net"),
    color.use = NULL,
    color.heatmap = NULL,
    title.name = NULL,
    width = NULL,
    height = NULL,
    font.size = 8,
    font.size.title = 10,
    cluster.rows = FALSE,
    cluster.cols = FALSE,
    sources.use = NULL,
    targets.use = NULL,
    remove.isolate = FALSE,
    row.show = NULL,
    col.show = NULL
) {
  measure <- match.arg(measure)
  slot.name <- match.arg(slot.name)
  
  # 获取 net.diff 矩阵
  if (is.list(object@net[[1]])) {
    if (is.null(color.heatmap)) color.heatmap <- c("#2166ac", "#b2182b")
    obj1     <- object@net[[comparison[1]]][[measure]]
    obj2     <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    legend.name <- if (measure=="count") "Differential number of interactions" else "Differential interaction strength"
    title.name  <- title.name %||% legend.name
  } else {
    if (is.null(color.heatmap)) color.heatmap <- "Reds"
    if (!is.null(signaling)) {
      net.diff   <- slot(object, slot.name)$prob[,,signaling]
      title.name <- title.name %||% paste0(signaling, " signaling network")
      legend.name <- "Communication Prob."
    } else {
      net.diff   <- object@net[[measure]]
      legend.name <- title.name %||% if (measure=="count") "Number of interactions" else "Interaction strength"
    }
  }
  
  net <- net.diff
  
  # 可选筛选 sources/targets
  if (!is.null(sources.use) || !is.null(targets.use)) {
    df.net <- melt(net, value.name="value")
    colnames(df.net)[1:2] <- c("source","target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) sources.use <- rownames(net)[sources.use]
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) targets.use <- colnames(net)[targets.use]
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$value[is.na(df.net$value)] <- 0
    net <- with(df.net, tapply(value, list(source, target), sum))
  }
  
  mat <- net
  mat[is.na(mat)] <- 0
  
  # 默认颜色
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  
  color.use.row <- color.use
  color.use.col <- color.use
  
  # 去除孤立节点
  if (remove.isolate) {
    zr <- which(rowSums(mat)==0)
    zc <- which(colSums(mat)==0)
    if (length(zr)>0) {
      mat <- mat[-zr, , drop=FALSE]
      color.use.row <- color.use.row[-zr]
    }
    if (length(zc)>0) {
      mat <- mat[, -zc, drop=FALSE]
      color.use.col <- color.use.col[-zc]
    }
  }
  
  # 根据用户指定行/列选取
  if (!is.null(row.show)) {
    mat <- mat[row.show, , drop=FALSE]
    color.use.row <- color.use.row[row.show]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show, drop=FALSE]
    color.use.col <- color.use.col[col.show]
  }
  
  # === 核心修改：固定热图颜色映射范围 ===
  # 无论实际数据如何，都使用 -15 到 15
  # 负正双色
  color.heatmap.use <- colorRamp3(
    c(-0.02, 0, 0.015),
    c(color.heatmap[1], "#f7f7f7", color.heatmap[2])
  )
  colorbar.break <- c(-0.02, 0, 0.015)
  # =======================================
  
  # 构建行/列注释
  df.col <- data.frame(group = colnames(mat), row.names=colnames(mat))
  df.row <- data.frame(group = rownames(mat), row.names=rownames(mat))
  col_anno <- HeatmapAnnotation(
    df = df.col, col=list(group=color.use.col),
    which="column", show_legend=FALSE, show_annotation_name=FALSE,
    simple_anno_size=unit(0.2,"cm")
  )
  row_anno <- HeatmapAnnotation(
    df = df.row, col=list(group=color.use.row),
    which="row", show_legend=FALSE, show_annotation_name=FALSE,
    simple_anno_size=unit(0.2,"cm")
  )
  ha_top <- HeatmapAnnotation(
    Strength = anno_barplot(colSums(abs(mat)),
                            border=FALSE, gp=gpar(fill=color.use.col, col=color.use.col)
    ), show_annotation_name=FALSE
  )
  ha_left <- rowAnnotation(
    Strength = anno_barplot(rowSums(abs(mat)),
                            border=FALSE, gp=gpar(fill=color.use.row, col=color.use.row)
    ), show_annotation_name=FALSE
  )
  
  # 绘图
  ht <- Heatmap(
    mat,
    name = legend.name,
    col = color.heatmap.use,
    na_col = "white",
    bottom_annotation = col_anno,
    left_annotation   = row_anno,
    top_annotation    = ha_top,
    right_annotation  = ha_left,
    cluster_rows      = cluster.rows,
    cluster_columns   = cluster.cols,
    row_names_side    = "left",
    row_names_rot     = 0,
    row_names_gp      = gpar(fontsize=font.size),
    column_names_gp   = gpar(fontsize=font.size),
    column_title      = title.name,
    column_title_gp   = gpar(fontsize=font.size.title),
    column_names_rot  = 90,
    row_title         = "Source",
    row_title_rot     = 90,
    row_title_gp      = gpar(fontsize=font.size.title),
    heatmap_legend_param = list(
      title_gp        = gpar(fontsize=font.size, fontface="plain"),
      title_position  = "leftcenter-rot",
      border          = NA,
      legend_height   = unit(20, "mm"),
      labels_gp       = gpar(fontsize=font.size),
      grid_width      = unit(2, "mm"),
      at              = colorbar.break
    )
  )
  
  return(ht)
}


plot_cellchat_heatmap_fixed_0_15 <- function (
    object,
    comparison = c(1, 2),
    measure = c("count", "weight"),
    signaling = NULL,
    slot.name = c("netP", "net"),
    color.use = NULL,
    color.heatmap = c("#2166ac", "#f7f7f7", "#b2182b"),
    title.name = NULL,
    width = NULL,
    height = NULL,
    font.size = 8,
    font.size.title = 10,
    cluster.rows = FALSE,
    cluster.cols = FALSE,
    sources.use = NULL,
    targets.use = NULL,
    remove.isolate = FALSE,
    row.show = NULL,
    col.show = NULL
) {
  # match args
  measure   <- match.arg(measure)
  slot.name <- match.arg(slot.name)
  
  # extract differential matrix
  if (is.list(object@net[[1]])) {
    message("Using merged comparisons\n")
    if (length(color.heatmap) != 3) stop("color.heatmap must be length 3 for fixed 0–15 range")
    obj1     <- object@net[[comparison[1]]][[measure]]
    obj2     <- object@net[[comparison[2]]][[measure]]
    mat      <- obj2 - obj1
    legend.name <- if (measure=="count") "Δ interaction count" else "Δ interaction strength"
    title.name  <- title.name %||% legend.name
  } else {
    message("Using single CellChat object\n")
    mat <- switch(
      TRUE,
      !is.null(signaling)     -> slot(object, slot.name)$prob[,,signaling],
      !is.null(measure)       -> object@net[[measure]],
      stop("Either signaling or measure must be provided")
    )
    title.name  <- title.name %||%
      if (!is.null(signaling)) paste0(signaling, " signaling") else
        if (measure=="count") "Interaction count" else "Interaction strength"
    legend.name <- title.name
  }
  
  # optional filter
  if (!is.null(sources.use) || !is.null(targets.use)) {
    df <- reshape2::melt(mat, value.name="value")
    colnames(df)[1:2] <- c("source","target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) sources.use <- rownames(mat)[sources.use]
      df <- subset(df, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) targets.use <- colnames(mat)[targets.use]
      df <- subset(df, target %in% targets.use)
    }
    df$value[is.na(df$value)] <- 0
    mat <- with(df, tapply(value, list(source,target), sum))
  }
  
  mat[is.na(mat)] <- 0
  
  # clamp to 0–15
  mat[mat <  0]  <- 0
  mat[mat > 0.02]  <- 0.02
  
  # default color.use
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  color.use.row <- color.use
  color.use.col <- color.use
  
  # remove isolates
  if (remove.isolate) {
    zr <- which(rowSums(mat)==0)
    zc <- which(colSums(mat)==0)
    if (length(zr)>0) {
      mat <- mat[-zr,,drop=FALSE]
      color.use.row <- color.use.row[-zr]
    }
    if (length(zc)>0) {
      mat <- mat[,-zc,drop=FALSE]
      color.use.col <- color.use.col[-zc]
    }
  }
  
  # subsetting rows/cols
  if (!is.null(row.show)) {
    mat <- mat[row.show,,drop=FALSE]
    color.use.row <- color.use.row[row.show]
  }
  if (!is.null(col.show)) {
    mat <- mat[,col.show,drop=FALSE]
    color.use.col <- color.use.col[col.show]
  }
  
  # create color function for 0–15
  color.heatmap.fun <- colorRamp3(
    c(0, 0.01, 0.02),
    color.heatmap
  )
  
  # annotations
  df.col <- data.frame(group=colnames(mat), row.names=colnames(mat))
  df.row <- data.frame(group=rownames(mat), row.names=rownames(mat))
  col_anno <- HeatmapAnnotation(
    df=df.col, col=list(group=color.use.col),
    which="column", show_legend=FALSE, show_annotation_name=FALSE,
    simple_anno_size=unit(0.2,"cm")
  )
  row_anno <- HeatmapAnnotation(
    df=df.row, col=list(group=color.use.row),
    which="row", show_legend=FALSE, show_annotation_name=FALSE,
    simple_anno_size=unit(0.2,"cm")
  )
  top_anno <- HeatmapAnnotation(
    Strength = anno_barplot(colSums(mat),
                            border=FALSE, gp=gpar(fill=color.use.col, col=color.use.col)
    ),
    show_annotation_name=FALSE
  )
  left_anno <- rowAnnotation(
    Strength = anno_barplot(rowSums(mat),
                            border=FALSE, gp=gpar(fill=color.use.row, col=color.use.row)
    ),
    show_annotation_name=FALSE
  )
  
  # draw heatmap
  ht <- Heatmap(
    mat,
    name                = legend.name,
    col                 = color.heatmap.fun,
    na_col              = "white",
    bottom_annotation   = col_anno,
    left_annotation     = row_anno,
    top_annotation      = top_anno,
    right_annotation    = left_anno,
    cluster_rows        = cluster.rows,
    cluster_columns     = cluster.cols,
    row_names_side      = "left",
    row_names_rot       = 0,
    row_names_gp        = gpar(fontsize=font.size),
    column_names_gp     = gpar(fontsize=font.size),
    column_title        = title.name,
    column_title_gp     = gpar(fontsize=font.size.title),
    column_names_rot    = 90,
    row_title           = "Source",
    row_title_rot       = 90,
    row_title_gp        = gpar(fontsize=font.size.title),
    heatmap_legend_param = list(
      at            = c(0,0.01,0.02),
      labels        = c("0","0.01","0.02"),
      title_gp      = gpar(fontsize=font.size, fontface="plain"),
      title_position= "leftcenter-rot",
      border        = NA,
      legend_height = unit(20,"mm"),
      labels_gp     = gpar(fontsize=font.size),
      grid_width    = unit(2,"mm")
    ),
    width  = width,
    height = height
  )
  
  return(ht)
}
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(grid)

plot_cellchat_heatmap_fixed_0_15 <- function (
    object,
    comparison = c(1, 2),
    measure = c("count", "weight"),
    signaling = NULL,
    slot.name = c("netP", "net"),
    color.use = NULL,
    color.heatmap = c("#2166ac", "#f7f7f7", "#b2182b"),
    title.name = NULL,
    width = NULL,
    height = NULL,
    font.size = 8,
    font.size.title = 10,
    cluster.rows = FALSE,
    cluster.cols = FALSE,
    sources.use = NULL,
    targets.use = NULL,
    remove.isolate = FALSE,
    row.show = NULL,
    col.show = NULL
) {
  # 参数匹配
  measure   <- match.arg(measure)
  slot.name <- match.arg(slot.name)
  
  # 提取矩阵
  if (is.list(object@net[[1]])) {
    message("Using merged comparisons\n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    mat  <- obj2 - obj1
    legend.name <- if (measure=="count") "Δ interaction count" else "Δ interaction strength"
    title.name  <- title.name %||% legend.name
  } else {
    message("Using single CellChat object\n")
    if (!is.null(signaling)) {
      mat <- slot(object, slot.name)$prob[,,signaling]
      title.name <- title.name %||% paste0(signaling, " signaling")
      legend.name <- "Communication Prob."
    } else {
      mat <- object@net[[measure]]
      legend.name <- title.name %||%
        if (measure=="count") "Interaction count" else "Interaction strength"
    }
  }
  
  # 筛选 sources/targets
  if (!is.null(sources.use) || !is.null(targets.use)) {
    df <- reshape2::melt(mat, value.name="value")
    colnames(df)[1:2] <- c("source","target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) sources.use <- rownames(mat)[sources.use]
      df <- subset(df, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) targets.use <- colnames(mat)[targets.use]
      df <- subset(df, target %in% targets.use)
    }
    df$value[is.na(df$value)] <- 0
    mat <- with(df, tapply(value, list(source,target), sum))
  }
  
  mat[is.na(mat)] <- 0
  
  # 截断到 0–15 范围
  mat[mat <  0]  <- 0
  mat[mat > 0.015]  <- 0.015
  
  # 处理 color.use 长度
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  if (length(color.use) == 1) {
    color.use <- rep(color.use, length.out = ncol(mat))
  } else if (length(color.use) != ncol(mat)) {
    stop(sprintf(
      "color.use 长度为 %d，不等于列数 %d，请提供 %d 个颜色。",
      length(color.use), ncol(mat), ncol(mat)
    ))
  }
  names(color.use) <- colnames(mat)
  color.use.row <- color.use
  color.use.col <- color.use
  
  # 去除孤立节点
  if (remove.isolate) {
    zr <- which(rowSums(mat)==0)
    zc <- which(colSums(mat)==0)
    if (length(zr)>0) {
      mat <- mat[-zr,,drop=FALSE]
      color.use.row <- color.use.row[-zr]
    }
    if (length(zc)>0) {
      mat <- mat[,-zc,drop=FALSE]
      color.use.col <- color.use.col[-zc]
    }
  }
  
  # 子集行/列
  if (!is.null(row.show)) {
    mat <- mat[row.show,,drop=FALSE]
    color.use.row <- color.use.row[row.show]
  }
  if (!is.null(col.show)) {
    mat <- mat[,col.show,drop=FALSE]
    color.use.col <- color.use.col[col.show]
  }
  
  # 创建 0–15 三色渐变
  if (length(color.heatmap) != 3) {
    stop("color.heatmap 必须是长度为 3 的颜色向量，例如 c('#2166ac','#f7f7f7','#b2182b')")
  }
  color.heatmap.fun <- colorRamp3(
    c(0, 0.005, 0.015),
    color.heatmap
  )
  
  # 注释条
  df.col <- data.frame(group=colnames(mat), row.names=colnames(mat))
  df.row <- data.frame(group=rownames(mat), row.names=rownames(mat))
  col_anno <- HeatmapAnnotation(
    df=df.col, col=list(group=color.use.col),
    which="column", show_legend=FALSE, show_annotation_name=FALSE,
    simple_anno_size=unit(0.2,"cm")
  )
  row_anno <- HeatmapAnnotation(
    df=df.row, col=list(group=color.use.row),
    which="row", show_legend=FALSE, show_annotation_name=FALSE,
    simple_anno_size=unit(0.2,"cm")
  )
  top_anno <- HeatmapAnnotation(
    Strength = anno_barplot(colSums(mat),
                            border=FALSE, gp=gpar(fill=color.use.col, col=color.use.col)
    ),
    show_annotation_name=FALSE
  )
  left_anno <- rowAnnotation(
    Strength = anno_barplot(rowSums(mat),
                            border=FALSE, gp=gpar(fill=color.use.row, col=color.use.row)
    ),
    show_annotation_name=FALSE
  )
  
  # 绘制热图
  ht <- Heatmap(
    mat,
    name                = legend.name,
    col                 = color.heatmap.fun,
    na_col              = "white",
    bottom_annotation   = col_anno,
    left_annotation     = row_anno,
    top_annotation      = top_anno,
    right_annotation    = left_anno,
    cluster_rows        = cluster.rows,
    cluster_columns     = cluster.cols,
    row_names_side      = "left",
    row_names_rot       = 0,
    row_names_gp        = gpar(fontsize=font.size),
    column_names_gp     = gpar(fontsize=font.size),
    column_title        = title.name,
    column_title_gp     = gpar(fontsize=font.size.title),
    column_names_rot    = 90,
    row_title           = "Source",
    row_title_rot       = 90,
    row_title_gp        = gpar(fontsize=font.size.title),
    heatmap_legend_param = list(
      at            = c(0,0.005,0.015),
      labels        = c("0","0.005","0.015"),
      title_gp      = gpar(fontsize=font.size, fontface="plain"),
      title_position= "leftcenter-rot",
      border        = NA,
      legend_height = unit(20,"mm"),
      labels_gp     = gpar(fontsize=font.size),
      grid_width    = unit(2,"mm")
    ),
    width  = width,
    height = height
  )
  
  return(ht)
}
