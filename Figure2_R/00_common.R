# ==============================================================================
# 00_common.R — 公共配置与工具函数
# 复现论文: Allen W.E., Blosser T.R., Sullivan Z.A., Dulac C., Zhuang X.
#   "Molecular and spatial signatures of mouse brain aging at single-cell
#    resolution." Cell 186(1):194-208.e18, 2023. DOI: 10.1016/j.cell.2022.12.010
# 复现内容: Figure 2 (A) 细胞类型密度随年龄变化 (B) 细胞状态比例随年龄变化
#                     (C) 各细胞状态示例基因表达小提琴图
# 数据: CELLxGENE 导出的 MERFISH 数据 (378,918 cells x 374 genes)
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------- 路径设置 ---------------------------------------
# 如果换机器运行, 只需修改这几个路径
PROJECT_DIR <- "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure2_R"
H5AD_FILE   <- "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/BrainAgingSpatialAtlas_MERFISH.h5ad"
# snRNA-seq (10x) 数据: 含全转录组 20,929 基因, 用于 Fig 2C 中不在 MERFISH panel
# 内的 15 个基因 (原文由 snRNA-seq 插补; 本项目用同状态 snRNA 细胞实测表达近似)
SNRNA_FILE  <- "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/BrainAgingSpatialAtlas_snRNAseq.h5ad"

CACHE_DIR <- file.path(PROJECT_DIR, "cache")
FIG_DIR   <- file.path(PROJECT_DIR, "figures")
TAB_DIR   <- file.path(PROJECT_DIR, "tables")
# 高分辨率图 (供 DOCX / PPTX 排版; 普通 figures/ 仍是 300 dpi)
HI_DIR    <- file.path(PROJECT_DIR, "figures_hires")
for (d in c(CACHE_DIR, FIG_DIR, TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

CACHE_FILE      <- file.path(CACHE_DIR, "merfish_cache.rds")
SNRNA_CACHE_FILE <- file.path(CACHE_DIR, "snrna_states_cache.rds")

# ----------------------------- 全局常量 ---------------------------------------
AGE_LEVELS <- c("4wk", "24wk", "90wk")
# 年龄配色 —— 取自作者原始代码 (权威来源):
#   analyze_merfish_celltyping.ipynb:  age_colors = ['cornflowerblue','thistle','lightcoral']
# 这三个是 matplotlib 的 X11 命名色, 对应下面三个 hex。
# 注: 之前用 tools/sample_paper_colors.py 从 PDF 图上采样得到 #7296D8/#D2C0D2/#DE8A8A,
#     那是 JPEG 压缩 + alpha 混合后冲淡的值, 应以代码里的原始值为准。
AGE_COLORS <- c("4wk" = "#6495ED",   # cornflowerblue
                "24wk" = "#D8BFD8",  # thistle
                "90wk" = "#F08080")  # lightcoral
AGE_LABELS <- c("4wk" = "Juvenile (4 wk)", "24wk" = "Young adult (24 wk)",
                "90wk" = "Old (90 wk)")
AGE_LABELS_SHORT <- c("4wk" = "Juvenile", "24wk" = "Young adult", "90wk" = "Old")

# 主要细胞类型排列顺序 (与原文 Fig 2A 一致)
CELLTYPE_ORDER <- c("ExN", "InN", "MSN", "Olig", "Astro", "Endo", "Micro",
                    "OPC", "Peri", "VLMC", "Epen", "Macro", "T cell")
# 低丰度细胞类型 (Fig 2A 一级 inset)
CELLTYPE_INSET <- c("VLMC", "Epen", "Macro", "T cell")
# 极低丰度 (Fig 2A 二级 inset)
CELLTYPE_INSET2 <- "T cell"

# Fig 2B 涉及的细胞类型及其状态 (子集群)
STATES_OF <- list(
  Olig  = c("Olig-1", "Olig-2", "Olig-3"),
  Micro = c("Micro-1", "Micro-2", "Micro-3"),
  Endo  = c("Endo-1", "Endo-2", "Endo-3"),
  Astro = c("Astro-1", "Astro-2")
)
# Fig 2B/C 状态配色: 同一大类内深浅渐变
#   Oligo 黑→深灰→浅灰; Micro 深紫红→亮粉; Endo 橄榄→黄; Astro 深绿→浅绿
# Fig 2B/C 状态配色: 逐块从原文 Figure 2 采样得到 (tools/sample_paper_colors.py)。
#   Olig  黑 → 中灰 → 浅灰
#   Micro 深紫红 → 品红 → 亮粉
#   Endo  橄榄 → 暗黄 → 亮黄
#   Astro 深绿 → 中绿
STATE_COLORS <- c(
  "Olig-1"  = "#211F1F", "Olig-2"  = "#919396", "Olig-3"  = "#CCCCD2",
  "Micro-1" = "#96185A", "Micro-2" = "#C91778", "Micro-3" = "#FC1191",
  "Endo-1"  = "#8A8454", "Endo-2"  = "#BAB46C", "Endo-3"  = "#EAE48A",
  "Astro-1" = "#246642", "Astro-2" = "#2A8A54"
)
ALL_STATES <- unlist(STATES_OF, use.names = FALSE)

# Fig 2C 示例基因 (与原文一致, 每类 8 个共 32 个)
# 不在 MERFISH 374 基因 panel 内的基因, 用 snRNA-seq 同状态细胞实测表达替代
# (原文为基于 snRNA-seq 的插补; 图中这些基因名后带 * 标注)
GENES_TO_SHOW <- list(
  Olig  = c("Ank2", "Trim2", "Robo1", "Spock1", "C4b", "Neat1", "Dgki", "Il33"),
  Micro = c("Dpp10", "Meg3", "Ctss", "Chrm3", "Trem2", "Trim30a", "Itgb2", "B2m"),
  Endo  = c("Ly6a", "Apoe", "B2m", "Rbm20", "Tmem209", "Xdh", "Fmo2", "Bgn"),
  Astro = c("Lsamp", "Gpc5", "Slc1a2", "Luzp2", "Trpm3", "Brinp3", "Gfap", "C4b")
)
# 缓存计数矩阵时要额外保留的基因 (供后续扩展分析; 修改后需删除 cache 重跑)
EXTRA_GENES <- character(0)

# ----------------------------- 工具函数 ---------------------------------------
#' 读取 h5ad 中 obs 的分类(categorical)列, 返回解码后的字符串向量
.read_obs_cat <- function(h5, col) {
  g <- h5[[paste0("obs/", col)]]
  codes <- g[["codes"]]$read()
  cats  <- g[["categories"]]$read()
  out <- rep(NA_character_, length(codes))
  keep <- codes >= 0
  out[keep] <- cats[codes[keep] + 1L]
  out
}

#' 从 h5ad 读取细胞元信息 + 指定基因的原始计数 (内存友好的分块读取)
#'
#' 说明: h5ad 的 raw/X 是 CSR 稀疏矩阵 (行=细胞 378,918, 列=基因 374,
#'       元素=原始计数)。为控制内存占用, 只读取 genes_keep 中指定的基因列。
#'       meta 的行顺序与 counts 完全一致 (即 h5ad 中的细胞顺序)。
read_h5ad_merfish <- function(path, genes_keep = NULL) {
  message("Reading: ", path)
  h5 <- hdf5r::H5File$new(path, "r")
  on.exit(h5$close_all())

  n_cells <- h5[["obs/_index"]]$dims

  meta <- data.frame(
    age         = factor(.read_obs_cat(h5, "age"), levels = AGE_LEVELS),
    # 用 cell_type (CELLxGENE 原始 13 类) 再映射到原文简写;
    # 不能用 cell_type_annot —— 那个字段把 Epen 归到 "cell", 会丢掉 Epen。
    cell_type   = normalise_celltype(.read_obs_cat(h5, "cell_type")),
    clust_annot = .read_obs_cat(h5, "clust_annot"),
    donor_id    = .read_obs_cat(h5, "donor_id"),
    slice       = .read_obs_cat(h5, "slice"),
    stringsAsFactors = FALSE
  )
  # Epen: CELLxGENE 把室管膜细胞的 cell_type 记为空 ("cell"), 需按 clust_annot 补回。
  # 原文也是从 Astro 里把 clust_annot=='Epen' 的细胞拆出来的。
  meta$cell_type[meta$clust_annot == "Epen"] <- "Epen"

  # 空间坐标 (obsm/spatial_coords: 2 x n_cells, 单位: 微米)
  xy <- h5[["obsm/spatial_coords"]]$read()
  meta$x <- xy[1, ]
  meta$y <- xy[2, ]
  rm(xy); invisible(gc(verbose = FALSE))

  # 基因符号 (var/feature_name 为分类编码)
  fn <- h5[["var/feature_name"]]
  gene_sym <- fn[["categories"]]$read()[fn[["codes"]]$read() + 1L]

  # ---- raw/X: CSR 稀疏矩阵, 分块读取并只保留目标基因列 ------------------------
  Xg <- h5[["raw/X"]]
  indptr <- as.double(Xg[["indptr"]]$read())
  stopifnot(length(indptr) == n_cells + 1L)

  wanted <- which(gene_sym %in% genes_keep)
  message("  cells: ", n_cells, " | panel genes: ", length(gene_sym),
          " | genes to keep: ", length(wanted),
          " (", paste(gene_sym[wanted], collapse = ", "), ")")

  counts <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = double(0),
                                 dims = c(n_cells, length(wanted)))
  if (length(wanted) > 0) {
    CH <- 100000L                       # 每块 10 万个细胞
    rows_all <- list(); cols_all <- list(); vals_all <- list(); k <- 0
    for (s in seq(1L, n_cells, by = CH)) {
      e <- min(s + CH - 1L, n_cells)
      a <- indptr[s] + 1L; b <- indptr[e + 1L]
      if (b >= a) {
        idx <- as.integer(Xg[["indices"]][a:b]) + 1L   # 基因列号 (1-based)
        keep <- idx %in% wanted
        if (any(keep)) {
          k <- k + 1
          vals_all[[k]] <- as.numeric(Xg[["data"]][a:b][keep])
          cols_all[[k]] <- match(idx[keep], wanted)
          rows_all[[k]] <- rep(seq(s, e), diff(indptr[s:(e + 1L)]))[keep]
        }
      }
      if (e %% 500000L == 0L) message("    ... ", e, " / ", n_cells, " cells")
    }
    counts <- Matrix::sparseMatrix(
      i = unlist(rows_all), j = unlist(cols_all), x = unlist(vals_all),
      dims = c(n_cells, length(wanted)),
      dimnames = list(NULL, gene_sym[wanted]))
    rm(rows_all, cols_all, vals_all); invisible(gc(verbose = FALSE))
  }

  list(meta = meta, counts = counts, genes = gene_sym,
       genes_kept = gene_sym[wanted])
}

#' 加载数据 (有缓存且缓存覆盖所需基因则直接读; 否则重建)
load_data <- function() {
  genes_keep <- unique(c(unlist(GENES_TO_SHOW, use.names = FALSE), EXTRA_GENES))
  if (file.exists(CACHE_FILE)) {
    obj <- readRDS(CACHE_FILE)
    # 缓存有效性: 所需基因中凡在 panel 内的, 都必须已在缓存里
    missing <- setdiff(intersect(genes_keep, obj$genes), obj$genes_kept)
    if (length(missing) == 0) {
      message("Loading cache: ", CACHE_FILE)
      return(obj)
    }
    message("Cache missing genes: ", paste(missing, collapse = ", "),
            " -> rebuilding ...")
  }
  obj <- read_h5ad_merfish(H5AD_FILE, genes_keep = genes_keep)
  saveRDS(obj, CACHE_FILE, compress = "gzip")
  message("Cache saved: ", CACHE_FILE)
  obj
}

#' 从 snRNA-seq (10x) h5ad 读取指定状态细胞的指定基因原始计数 (分块, 内存友好)
#'
#' snRNA-seq 与 MERFISH 使用同一套联合聚类注释 (clust_annot, 43 个状态),
#' 因此 Fig 2C 中 MERFISH panel 外的基因用同状态 snRNA 细胞的实测表达近似
#' (原文为 snRNA-seq 插补)。返回 list(meta, counts, genes):
#'   meta$clust / meta$age / meta$total_counts (用于 CP10K 归一化)
#'   counts: dgCMatrix (目标细胞 x 指定基因, 原始计数)
load_snrna_states <- function(genes_needed, states = ALL_STATES) {
  if (file.exists(SNRNA_CACHE_FILE)) {
    message("Loading snRNA cache: ", SNRNA_CACHE_FILE)
    obj <- readRDS(SNRNA_CACHE_FILE)
    if (!all(genes_needed %in% obj$genes)) {
      stop("snRNA 缓存缺少基因: ",
           paste(setdiff(genes_needed, obj$genes), collapse = ", "),
           " — 请删除 ", SNRNA_CACHE_FILE, " 后重跑")
    }
    return(obj)
  }
  if (!file.exists(SNRNA_FILE)) {
    stop("未找到 snRNA-seq 文件: ", SNRNA_FILE,
         "\n请从 CELLxGENE 下载: 752ed29d-bbb0-45bd-9c57-19edf68ef73e.h5ad")
  }
  message("Reading snRNA-seq: ", SNRNA_FILE)
  h5 <- hdf5r::H5File$new(SNRNA_FILE, "r")
  on.exit(h5$close_all())

  n_cells <- h5[["obs/_index"]]$dims
  clust <- .read_obs_cat(h5, "clust_annot")
  age   <- .read_obs_cat(h5, "age")
  total_counts <- h5[["obs/total_counts"]]$read()

  fn <- h5[["var/feature_name"]]
  gene_sym <- fn[["categories"]]$read()[fn[["codes"]]$read() + 1L]
  wanted <- which(gene_sym %in% genes_needed)
  message("  cells: ", n_cells, " | genes: ", length(gene_sym),
          " | genes to keep: ", length(wanted))

  Xg <- h5[["raw/X"]]          # CSR: 行=细胞, 列=基因
  indptr <- as.double(Xg[["indptr"]]$read())
  is_target <- clust %in% states

  rows_all <- list(); cols_all <- list(); vals_all <- list(); k <- 0
  CH <- 2000L                  # 每块 2000 个细胞, 控制单块内存 <100MB
  for (s in seq(1L, n_cells, by = CH)) {
    e <- min(s + CH - 1L, n_cells)
    if (!any(is_target[s:e])) next          # 跳过无目标细胞的块
    a <- indptr[s] + 1L; b <- indptr[e + 1L]
    idx <- as.integer(Xg[["indices"]][a:b]) + 1L
    row_of_entry <- rep(s:e, diff(indptr[s:(e + 1L)]))
    keep <- idx %in% wanted & is_target[row_of_entry]
    if (any(keep)) {
      k <- k + 1
      vals_all[[k]] <- as.numeric(Xg[["data"]][a:b][keep])
      cols_all[[k]] <- match(idx[keep], wanted)
      rows_all[[k]] <- row_of_entry[keep]
    }
    if (e %% 20000L == 0L) message("    ... ", e, " / ", n_cells)
  }
  counts <- Matrix::sparseMatrix(
    i = unlist(rows_all), j = unlist(cols_all), x = unlist(vals_all),
    dims = c(n_cells, length(wanted)),
    dimnames = list(NULL, gene_sym[wanted]))
  rm(rows_all, cols_all, vals_all); invisible(gc(verbose = FALSE))

  keep_rows <- which(is_target)
  obj <- list(
    meta = data.frame(clust = clust[keep_rows], age = age[keep_rows],
                      total_counts = total_counts[keep_rows],
                      stringsAsFactors = FALSE),
    counts = counts[keep_rows, , drop = FALSE],
    genes = gene_sym[wanted]
  )
  saveRDS(obj, SNRNA_CACHE_FILE, compress = "gzip")
  message("snRNA cache saved: ", SNRNA_CACHE_FILE)
  obj
}

#' 多边形面积 (shoelace 公式)
poly_area <- function(x, y) {
  n <- length(x)
  abs(sum(x * y[c(2:n, 1)] - x[c(2:n, 1)] * y)) / 2
}

#' 多边形周长
poly_perimeter <- function(x, y) {
  n <- length(x)
  p <- cbind(x, y); p <- rbind(p, p[1, ])
  sum(sqrt(diff(p[, 1])^2 + diff(p[, 2])^2))
}

#' 均值 ± 95% CI (t 分布)
mean_ci <- function(x) {
  n <- length(x); m <- mean(x); s <- sd(x)
  half <- if (n > 1) qt(0.975, n - 1) * s / sqrt(n) else 0
  c(mean = m, lower = m - half, upper = m + half, n = n)
}

#' 两组独立样本 t 检验 (equal variance, 与 scipy.stats.ttest_ind 默认一致)
ttest_p <- function(v1, v2) {
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  if (length(v1) < 2 || length(v2) < 2) return(NA_real_)
  tryCatch(t.test(v1, v2, var.equal = TRUE)$p.value, error = function(e) NA_real_)
}

# ==============================================================================
# 二维核密度小提琴 (diamond violin)
#
# 原文 Figure 2C 的形状不是 geom_violin 那种"宽度 ∝ 一维密度"的常规小提琴,
# 而是在 (表达量, 抖动后的水平坐标) 两个维度上做**二维高斯核密度估计**,
# 再取最外层填充等高线得到的菱形/纺锤形轮廓 (即 seaborn 的
# `sns.kdeplot(x=expr, y=jitter, fill=True)` 画法).
#
# 本机 R 环境没有 ks / KernSmooth 之外的二维 KDE 包, 且作者希望流程图完全基于
# 基础 R + ggplot2, 因此这里用 binning + 可分离卷积实现二维高斯 KDE,
# 再用 marching-squares 提取最外层等高线并归一化宽度 —— 计算量与细胞数近似线性。
# ==============================================================================

#' 极值归一化到 [0, 1]
.norm01 <- function(v) {
  r <- range(v, finite = TRUE)
  if (!is.finite(r[1]) || diff(r) <= 0) return(rep(0.5, length(v)))
  (v - r[1]) / diff(r)
}

#' 1D 高斯核密度估计 (用于确定网格范围与带宽)
.kde1d <- function(v, n = 512L, bw = NULL) {
  v <- v[is.finite(v)]
  if (length(v) < 2L) v <- c(v, v + 1e-6)
  if (is.null(bw)) bw <- stats::bw.nrd0(v)
  if (!is.finite(bw) || bw <= 0) bw <- max(stats::sd(v), 1e-6) * 0.5
  r <- range(v)
  pad <- 3 * bw
  g <- seq(r[1] - pad, r[2] + pad, length.out = n)
  d <- vapply(v, function(vi) stats::dnorm(g, vi, bw), numeric(n))
  list(x = g, y = if (is.null(dim(d))) d else rowMeans(d), bw = bw)
}

#' [内部] 可分离高斯平滑: 先按行, 再按列
.sep_smooth <- function(A, bwx, bwy) {
  lx <- max(1L, ceiling(4 * bwx)); kx <- exp(-0.5 * (seq(-lx, lx) / bwx)^2); kx <- kx / sum(kx)
  ly <- max(1L, ceiling(4 * bwy)); ky <- exp(-0.5 * (seq(-ly, ly) / bwy)^2); ky <- ky / sum(ky)
  n_row <- nrow(A); n_col <- ncol(A)
  Ar <- A
  for (s in seq_len(length(kx))) {
    sh <- s - lx - 1L
    if (sh == 0L) next
    if (sh > 0L) Ar[, (1 + sh):n_col] <- Ar[, (1 + sh):n_col] + kx[s] * A[, 1:(n_col - sh)]
    else         Ar[, 1:(n_col + sh)] <- Ar[, 1:(n_col + sh)] + kx[s] * A[, (1 - sh):n_col]
  }
  At <- Ar
  for (s in seq_len(length(ky))) {
    sh <- s - ly - 1L
    if (sh == 0L) next
    if (sh > 0L) At[(1 + sh):n_row, ] <- At[(1 + sh):n_row, ] + ky[s] * Ar[1:(n_row - sh), ]
    else         At[1:(n_row + sh), ] <- At[1:(n_row + sh), ] + ky[s] * Ar[(1 - sh):n_row, ]
  }
  At
}

#' 二维核密度估计 (binning + 可分离卷积), 返回规则网格上的密度
kde2d_binned <- function(x, y, n = c(140L, 180L), bw_mult = 1,
                         xlim = NULL, ylim = NULL, lev_q = 0.10, lev_adj = 1) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 5L) return(NULL)

  bw_x <- stats::bw.nrd0(x) * bw_mult
  bw_y <- stats::bw.nrd0(y) * bw_mult
  if (!is.finite(bw_x) || bw_x <= 0) bw_x <- max(stats::sd(x), 1e-6) * 0.3
  if (!is.finite(bw_y) || bw_y <= 0) bw_y <- max(stats::sd(y), 1e-6) * 0.3

  if (is.null(xlim)) xlim <- range(x) + c(-4, 4) * bw_x
  if (is.null(ylim)) ylim <- range(y) + c(-4, 4) * bw_y
  gx <- seq(xlim[1], xlim[2], length.out = n[2])
  gy <- seq(ylim[1], ylim[2], length.out = n[1])
  dx <- gx[2] - gx[1]; dy <- gy[2] - gy[1]

  ix <- pmin(n[2], pmax(1L, floor((x - gx[1]) / dx) + 1L))
  iy <- pmin(n[1], pmax(1L, floor((y - gy[1]) / dy) + 1L))
  M  <- matrix(0, nrow = n[1], ncol = n[2])
  idx <- (iy - 1L) * n[2] + ix
  M[idx] <- M[idx] + 1
  M <- .sep_smooth(M, bw_x / dx, bw_y / dy)
  s <- sum(M)
  if (s > 0) M <- M / s

  # 等高线水平: 取"在数据点处的二维密度"的低分位数。
  # 这样网格范围 (以及由此决定的 y 方向长宽比) 不会影响轮廓形状,
  # 得到的等高线恰好包含约 (1 - lev_q) 的数据点, 与原文手动 jitter + kdeplot
  # 得到的纺锤形轮廓等价。
  zpts <- M[cbind(iy, ix)]
  lv <- as.numeric(stats::quantile(zpts, lev_q, names = FALSE, type = 7))
  if (is.finite(lv) && lv > 0) {
    lv <- lv * lev_adj
  } else {
    lv <- max(M) * 0.05
  }
  list(x = gx, y = gy, z = M, dx = dx, dy = dy, level = lv)
}

#' 单变量核密度估计 (用于绘制小提琴宽度)
.uni_kde <- function(v, n = 128L, bw = NULL, from = NULL, to = NULL) {
  v <- v[is.finite(v)]
  if (length(v) < 2L) return(NULL)
  if (is.null(bw)) bw <- stats::bw.nrd0(v)
  if (!is.finite(bw) || bw <= 0) bw <- max(stats::sd(v), 1e-6) * 0.5
  if (is.null(from)) from <- min(v) - 2 * bw
  if (is.null(to))   to   <- max(v) + 2 * bw
  g <- seq(from, to, length.out = n)
  d <- rowMeans(vapply(v, function(vi) stats::dnorm(g, vi, bw), numeric(n)))
  data.frame(y = g, w = d / max(d))
}

#' [内部] marching squares: 逐"情形(case)"向量化提取 z = level 的等值线段
#'
#' 单元 (i, j) 四角: LT = z[i, j] (bit1), RT = z[i, j+1] (bit2),
#'                   RB = z[i+1, j+1] (bit4), LB = z[i+1, j] (bit8)。
#' 每条等值线段的端点落在单元边上, 位置由线性插值确定。
.ms_lines <- function(z, level) {
  nr <- nrow(z); nc <- ncol(z)
  LT <- z[-nr, -nc]; RT <- z[-nr, -1]; RB <- z[-1, -1]; LB <- z[-1, -nc]
  case <- (LT >= level) * 1L + (RT >= level) * 2L +
          (RB >= level) * 4L + (LB >= level) * 8L
  # 四条边上的插值参数 (0 = 上/右角, 1 = 下/左角)
  tT <- .interp(level, LT, RT)     # 上边 LT -> RT
  bB <- .interp(level, LB, RB)     # 下边 LB -> RB
  rR <- .interp(level, RT, RB)     # 右边 RT -> RB
  lL <- .interp(level, LT, LB)     # 左边 LT -> LB
  rc <- which(case > 0L, arr.ind = TRUE)
  if (nrow(rc) == 0L) return(NULL)
  i <- rc[, 1]; j <- rc[, 2]
  cd <- case[rc]
  pT <- data.frame(x = j + .at(tT, i, j),       y = i)
  pB <- data.frame(x = j + .at(bB, i, j),       y = i + 1)
  pR <- data.frame(x = j + 1,                   y = i + .at(rR, i, j))
  pL <- data.frame(x = j,                       y = i + .at(lL, i, j))
  sel <- function(a, b)
    data.frame(x0 = a[, 1], y0 = a[, 2], x1 = b[, 1], y1 = b[, 2])
  # 情形 5 / 10 是鞍点, 拆成两条线段
  parts <- list(
    c1  = sel(pL, pT), c2  = sel(pT, pR), c3  = sel(pL, pR), c4  = sel(pR, pB),
    c6  = sel(pT, pB), c7  = sel(pL, pB), c8  = sel(pB, pL), c9  = sel(pB, pT),
    c11 = sel(pB, pR), c12 = sel(pR, pL), c13 = sel(pT, pR), c14 = sel(pT, pL))
  c5  <- rbind(sel(pL, pT), sel(pR, pB))
  c10 <- rbind(sel(pT, pR), sel(pB, pL))
  cl <- c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8",
          "c9", "c10", "c11", "c12", "c13", "c14")
  out <- lapply(seq_along(cl), function(k) {
    hit <- cd == k
    if (!any(hit)) return(NULL)
    d <- switch(cl[k],
      c1 = parts$c1, c2 = parts$c2, c3 = parts$c3, c4 = parts$c4, c5 = c5,
      c6 = parts$c6, c7 = parts$c7, c8 = parts$c8, c9 = parts$c9, c10 = c10,
      c11 = parts$c11, c12 = parts$c12, c13 = parts$c13, c14 = parts$c14)
    d[hit, , drop = FALSE]
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0L) return(NULL)
  out
}

#' [内部] 线性插值参数 (a、b 同形)
.interp <- function(level, a, b) {
  p <- (level - a) / (b - a)
  p[!is.finite(p)] <- 0.5
  pmin(1, pmax(0, p))
}

#' [内部] 用 (i, j) 索引矩阵
.at <- function(M, i, j) M[cbind(i, j)]

#' [内部] 沿字段方向把线段串成闭合多边形, 返回最长的一条
#'
#' 关键: marching-squares 的相邻线段共享顶点在浮点上**完全相等**
#' (同一个插值表达式), 但直接用 `==` 比较仍会因累加顺序不同而失败,
#' 因此先把顶点坐标量化到 1e-6 的整数格点再建邻接表。
.ms_ring <- function(segs, gx, gy) {
  n <- nrow(segs)
  if (is.null(n) || n < 4L) return(NULL)
  f <- 1e6
  key <- function(x, y) paste0(round(x * f), "_", round(y * f))
  k0 <- key(segs$x0, segs$y0)
  k1 <- key(segs$x1, segs$y1)
  # 顶点 -> 线段编号 (每个内部顶点至多 2 条)
  nb <- new.env(hash = TRUE, parent = emptyenv())
  for (s in seq_len(n)) {
    for (k in c(k0[s], k1[s])) {
      cur <- nb[[k]]
      nb[[k]] <- c(cur, s)
    }
  }
  vert_of <- function(k) {
    p <- as.numeric(strsplit(k, "_", fixed = TRUE)[[1]]) / f
    c(p[1], p[2])
  }
  next_pt <- function(seg_i, kfrom) {
    if (identical(k0[seg_i], kfrom)) k1[seg_i] else k0[seg_i]
  }
  used <- logical(n)
  best <- NULL; best_len <- 0L
  for (s0 in seq_len(n)) {
    if (used[s0]) next
    used[s0] <- TRUE
    ks <- c(k0[s0], k1[s0])
    guard <- 0L
    repeat {
      guard <- guard + 1L
      if (guard > n + 5L) break
      kcur <- ks[length(ks)]
      if (length(ks) > 2L && identical(kcur, ks[1])) break
      cand <- nb[[kcur]]
      cand <- cand[!used[cand]]
      if (length(cand) == 0L) break
      # 优先选真正与原线段共享端点的线段
      pick <- cand[1]
      used[pick] <- TRUE
      ks <- c(ks, next_pt(pick, kcur))
    }
    if (length(ks) - 1L > best_len) {
      best_len <- length(ks) - 1L
      best <- vapply(ks, vert_of, numeric(2))
      best <- t(best)
    }
  }
  if (is.null(best) || best_len < 4L) return(NULL)
  # 线段端点保存的是**分数网格坐标**, 这里线性映射回数据坐标
  # (不能用 round() 取整, 否则会把几十个网格行压成同一行, 轮廓退化)
  xa <- gx[1] + (best[, 1] - 1) * (gx[2] - gx[1])
  ya <- gy[1] + (best[, 2] - 1) * (gy[2] - gy[1])
  keep <- !duplicated(cbind(round(xa, 9), round(ya, 9)))
  data.frame(x = xa[keep], y = ya[keep])
}

#' 二维核密度 "菱形小提琴" 轮廓
#'
#' 复刻原文 Figure 2C 的画法: 把每个基因的表达值沿水平方向加抖动作出二维散点,
#' 再做二维高斯核密度估计, 取一条最外层填充等高线作为形状轮廓, 宽度按该轮廓
#' 自身的最大宽度归一化到 [0, 1] (因此每个 gene × state 的形状都是紧凑的纺锤/
#' 菱形, 而不是跨越整个 jitter 区间的长条)。
#'
#' @param x,y 原始数据点 (x = 表达量, y = 抖动后的水平坐标)
#' @param lev_q 等高线在数据点密度上的分位数 (0.10 = 约含 90% 的数据点)
#' @param bw_mult 核带宽倍数 (越大越平滑)
#' @param n 网格尺寸 c(行 = y 方向, 列 = x 方向)
#' @return data.frame(y = 表达量, w = 归一化宽度); 失败返回 NULL
diamond_outline <- function(x, y, lev_q = 0.10, bw_mult = 1, n = c(140L, 180L),
                            lev_adj = 1) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 20L) return(NULL)
  d <- kde2d_binned(x, y, n = n, xlim = range(x), ylim = range(y),
                    bw_mult = bw_mult, lev_q = lev_q, lev_adj = lev_adj)
  if (is.null(d) || !is.finite(d$level)) return(NULL)
  segs <- .ms_lines(d$z, d$level)
  if (is.null(segs) || nrow(segs) < 4L) return(NULL)
  ring <- .ms_ring(segs, d$x, d$y)
  if (is.null(ring) || nrow(ring) < 3L) return(NULL)
  # 同一表达高度可能出现两个宽度 (轮廓回折); 取每侧最大宽度以获得凸的纺锤形
  agg <- stats::aggregate(w ~ ybin, data = data.frame(ybin = ring$y, w = ring$x),
                          FUN = max)
  agg <- agg[order(agg$ybin), ]
  if (nrow(agg) < 3L) return(NULL)
  data.frame(y = agg$ybin, w = agg$w / max(agg$w))
}

#' 1D 核密度小提琴轮廓: 返回 data.frame(y, w), w ∈ [0, 1] 已按峰值归一化
#'
#' 说明 (为什么最终用 1D KDE 而不是 2D KDE 轮廓):
#' 原文 Figure 2C 用 `sns.kdeplot` 在 (表达量, 水平抖动) 二维空间上做密度填充。
#' 忠实复刻时发现: MERFISH 表达矩阵极度稀疏 (例如 Olig-3 中 C4b 仅约 4% 细胞
#' 为阳性), y 方向 (抖动) 上绝大部分概率质量来自 "表达量 = 0" 的点质量, 它沿 y
#' 方向近似均匀分布, 因此**任何**能包含足够概率质量的二维等高线都会被拉成横跨整个
#' jitter 区间的长条, 无法得到原文那种紧凑的纺锤形 (见 README "与原文的差异" 第 6 条)。
#' 因此这里改为: 对表达量做 1D 高斯核密度估计, 再上下对称地以宽度 = 密度 画填充
#' 多边形 —— 即 geom_violin 的构造原理, 形状与原文 Figure 2C 一致, 但不受
#' 零膨胀点质量的影响。
violin_outline <- function(x, n = 160L, bw_mult = 1, from = NULL, to = NULL,
                           trim = 0.005, pad_frac = 0.15) {
  x <- x[is.finite(x)]
  if (length(x) < 20L) return(NULL)
  bw <- stats::bw.nrd0(x) * bw_mult
  if (!is.finite(bw) || bw <= 0) bw <- max(stats::sd(x), 1e-6) * 0.5
  r <- range(x)
  # 求值区间不能低于数据最小值: 表达量 >= 0, 若把区间向下延伸到负值, KDE 会因
  # 边界外的核质量而在 0 处形成一个虚假的"平底", 使小提琴看起来坐在地上。
  if (is.null(from)) from <- r[1]
  if (is.null(to))   to   <- r[2] + pad_frac * diff(r) + 1e-9
  g <- seq(from, to, length.out = n)
  d <- rowMeans(vapply(x, function(xi) stats::dnorm(g, xi, bw), numeric(n)))
  # 去掉密度极低的尾部, 避免多边形出现超长尖刺
  keep <- d >= max(d) * trim
  if (sum(keep) < 3L) keep <- d >= max(d) * 1e-6
  data.frame(y = g[keep], w = d[keep] / max(d))
}

#' 由 1D KDE 轮廓生成对称填充多边形
#'
#' @param x 表达向量
#' @param from,to KDE 求值区间 (通常取同一大类同一基因所有状态的公共范围)
violin_shape <- function(x, n = 160L, bw_mult = 1, from = NULL, to = NULL) {
  o <- violin_outline(x, n = n, bw_mult = bw_mult, from = from, to = to)
  if (is.null(o)) return(NULL)
  o
}

#' 把 violin_shape() 的轮廓展开成对称多边形
#'
#' @param o violin_shape() 的返回值
#' @param centre 该格在 **面板内** 的 x 位置 (相对序号 1..k, 注意不是全局状态序号)
#' @param half_width 小提琴半宽 (相对列间距)
violin_poly <- function(o, centre, half_width = 0.42) {
  if (is.null(o)) return(NULL)
  centre <- as.numeric(centre)[1]
  pos <- data.frame(x = centre + o$w * half_width, y = o$y)
  neg <- data.frame(x = centre - rev(o$w) * half_width, y = rev(o$y))
  rbind(pos, neg)
}

#' 统一细胞类型命名: CELLxGENE 的 cell_type 原始标签 -> 原文简写
#'
#' CELLxGENE 导出的 MERFISH 数据 `cell_type` 用的是本体名 (neuron / microglial
#' cell / endothelial cell / vascular leptomeningeal cell ...), 原文代码与
#' 本项目其它地方用的是 ExN / Micro / Endo / VLMC 这类简写。
#' 不做映射会导致所有细胞类型都对不上、柱子全空。
normalise_celltype <- function(x) {
  map <- c(
    "neuron"                        = "ExN",
    "inhibitory interneuron"        = "InN",
    "medium spiny neuron"           = "MSN",
    "oligodendrocyte"               = "Olig",
    "astrocyte"                     = "Astro",
    "endothelial cell"              = "Endo",
    "microglial cell"               = "Micro",
    "oligodendrocyte precursor cell"= "OPC",
    "pericyte"                      = "Peri",
    "vascular leptomeningeal cell"  = "VLMC",   # 原文代码写作 "Vlmc"
    "ependymal cell"                = "Epen",
    "macrophage"                    = "Macro",
    "T cell"                        = "T cell",
    "Vlmc"                          = "VLMC"    # 兼容旧缓存
  )
  out <- unname(map[as.character(x)])
  out[is.na(out) & as.character(x) %in% CELLTYPE_ORDER] <- as.character(x)[
    is.na(out) & as.character(x) %in% CELLTYPE_ORDER]
  out
}

#' 供体均值 + **bootstrap 95% CI** (对齐 seaborn barplot 的默认误差线)
#'
#' 原文用 `sns.barplot(...)`, 其默认误差线是 95% **bootstrap** 置信区间
#' (n_boot = 1000), 不是 t 分布 CI。n 只有 4/3/5 时两者差别明显, 会直接
#' 改变显著性判断, 因此必须按 bootstrap 复现。
#'
#' @param x 该组各供体的观测值 (如 cells/mm²)
#' @param n_boot 重采样次数 (seaborn 默认 1000)
#' @param seed 随机种子 (原文未设, 这里固定以保证可重复)
boot_ci <- function(x, n_boot = 1000L, seed = 20230105) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(low = NA_real_, high = NA_real_))
  if (length(x) == 1) return(c(low = x, high = x))
  set.seed(seed)
  b <- vapply(seq_len(n_boot),
              function(i) mean(sample(x, length(x), replace = TRUE)),
              numeric(1))
  q <- stats::quantile(b, c(0.025, 0.975), names = FALSE, type = 7)
  c(low = q[1], high = q[2])
}

#' ggplot 主题
theme_paper <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_fig <- function(plot, name, w, h, dpi = 300, cairo = FALSE,
                     hi_dpi = 900) {
  # cairo = TRUE 时用 cairo_pdf 设备导出 PDF (嵌入字体、跨平台一致;
  # 本项目要求 Arial 无衬线字体, cairo_pdf 能正确嵌入)。
  if (cairo) {
    grDevices::cairo_pdf(file.path(FIG_DIR, paste0(name, ".pdf")),
                         width = w, height = h)
    print(plot)
    grDevices::dev.off()
  } else {
    ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = w, height = h)
  }
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = w, height = h,
         dpi = dpi)
  # 另存一份高分辨率 PNG 供报告/幻灯片排版使用。
  # 默认 300 dpi 的 PNG 在 Word 里被放大到 16 cm 宽时有效分辨率只有约
  # 570 dpi, 进一步放大会明显发虚; 这里额外导出 hi_dpi 版本。
  # 注意 ggsave 会限制超高 dpi, 用 ragg/png 设备直接指定像素尺寸更稳妥。
  if (!is.null(hi_dpi) && hi_dpi > dpi) {
    dir.create(HI_DIR, showWarnings = FALSE, recursive = TRUE)
    px_w <- round(w * hi_dpi); px_h <- round(h * hi_dpi)
    grDevices::png(file.path(HI_DIR, paste0(name, ".png")),
                   width = px_w, height = px_h, res = hi_dpi, type = "cairo")
    print(plot)
    grDevices::dev.off()
  }
  message("Saved figure: ", name, "  (", dpi, " dpi",
          if (!is.null(hi_dpi) && hi_dpi > dpi) paste0(" + ", hi_dpi, " dpi hi-res") else "",
          ")")
}
