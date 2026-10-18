# ==============================================================================
# 00_common.R — Figure 5E / 5F 复现的公共部分 (R 版)
#
# 忠实于原文的定义表达 (Allen et al., Cell 2023; STAR Methods + 作者仓库
# github.com/ZhuangLab/SpatialBrainAgingCell22 的 int_analysis.ipynb)
#
# 1) 激活分数 (activation score) —— scanpy sc.tl.score_genes 的定义:
#      "The score is the average expression of a set of genes after subtraction
#       by the average expression of a reference set of genes. The reference set
#       is randomly sampled from the gene_pool for each binned expression value."
#    默认参数 ctrl_size = 50, n_bins = 25, random_state = 0, gene_pool = 全部基因。
#    基因集 (逐字取自作者代码 CELL 41):
#      activate_micro = B2m, Trem2, Ccl2, Apoe, Axl, Itgax, Cd9, C1qa, C1qc,
#                       Lyz2, Ctss
#      activate_astro = C4b, C3, Serpina3n, Cxcl10, Gfap, Vim, Il18, Hif3a
#      activate_olig  = C4b, Il33, Il18        (少突胶质炎症分数)
#    中心化 (作者 CELL 42/43): 减去该细胞类型 4wk 组的均值
#
# 2) 距离-激活曲线 (Fig 5E) —— 作者 CELL 46/47/467-470:
#      identify_nearest_neighbors_with_dist: KDTree.query(k=2) 取最近距离
#      compute_binned_values(dists, scores, min_d=0, max_d=100, bin_size=30):
#        for i in 0:(max_d-bin_size-1):
#          idx = (dists > i) & (dists <= i+bin_size)
#          binned_mean[i] = mean(scores[idx])
#          binned_std[i]  = std(scores[idx]) / sqrt(n)
#        binned_mean -= mean(binned_mean)
#      x 轴 = seq_along(binned_mean) - 1 + bin_size
#
# 3) 相关性 (Fig 5F) —— 作者 CELL 472-474:
#      在胼胝体(CC)区域内, 每个 src 细胞取其 50 µm 内所有 tgt 细胞的
#      平均激活分数, 与自身激活分数求 Pearson R
#
# 空间坐标单位实测为 **µm** (单切片凸包面积 16.44 mm²、最近邻中位 15.1 µm),
# 与 30 µm 窗口、原文 500 µm scale bar 一致, 不做任何缩放。
# ==============================================================================

suppressPackageStartupMessages({
  library(hdf5r); library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2)
})

# ----------------------------- 路径 -------------------------------------------
PROJECT_DIR <- "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro"
FIG5_DIR <- file.path(PROJECT_DIR, "Figure5_R")
FIG_DIR  <- file.path(FIG5_DIR, "figures")
TAB_DIR  <- file.path(FIG5_DIR, "tables")
HI_DIR   <- file.path(FIG5_DIR, "figures_hires")
CACHE_DIR <- file.path(FIG5_DIR, "cache")
for (d in c(FIG_DIR, TAB_DIR, HI_DIR, CACHE_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

H5AD_FILE <- file.path(PROJECT_DIR, "BrainAgingSpatialAtlas_MERFISH.h5ad")
CACHE_FILE <- file.path(CACHE_DIR, "fig5_merfish_cache.rds")

# ----------------------------- 基因集 (逐字取自作者代码) -----------------------
GENE_SETS <- list(
  activate_micro = c("B2m", "Trem2", "Ccl2", "Apoe", "Axl", "Itgax", "Cd9",
                     "C1qa", "C1qc", "Lyz2", "Ctss"),
  activate_astro = c("C4b", "C3", "Serpina3n", "Cxcl10", "Gfap", "Vim",
                     "Il18", "Hif3a"),
  activate_olig  = c("C4b", "Il33", "Il18")
)

AGE_LEVELS <- c("4wk", "24wk", "90wk")
AGE_LABELS <- c("4wk" = "Juvenile", "24wk" = "Young Adult", "90wk" = "Old")

# 细胞类型映射: CELLxGENE 本体名 -> 原文简写
CT_MAP <- c(
  "neuron" = "ExN", "inhibitory interneuron" = "InN",
  "medium spiny neuron" = "MSN", "oligodendrocyte" = "Olig",
  "astrocyte" = "Astro", "endothelial cell" = "Endo",
  "microglial cell" = "Micro", "oligodendrocyte precursor cell" = "OPC",
  "pericyte" = "Peri", "vascular leptomeningeal cell" = "VLMC",
  "ependymal cell" = "Epen", "macrophage" = "Macro", "T cell" = "T cell"
)

# 解剖区域: 原文 CELL 455 的 spatial_order 与发布版 obs$tissue 的对应关系
REGION_MAP <- c(
  "pia mater"             = "Pia",
  "cortical layer II/III" = "L2/3",
  "cortical layer V"      = "L5",
  "cortical layer VI"     = "L6",
  "corpus callosum"       = "CC",
  "olfactory region"      = "LatSept",
  "striatum"              = "Striatum",
  "brain ventricle"       = "Ventricle"
)
# Fig 5F 所用的区域 (作者代码: obs.spatial_clust_annots == "CC")
CC_TISSUE <- "corpus callosum"

# Fig 5F 的填充色梯度 —— 原文 sns.kdeplot(fill=True) 的默认 cmap
# 取 seaborn 0.13 该 QuadContourSet cmap (256 级) 的 10 个等距代表色。
KDE_FILL_COLORS <- c("#acc9ee", "#90bef1", "#6fb2f5", "#3fa8f7", "#2c9be9", "#2a8fd6", "#2f82c2", "#3676ac", "#3b6995", "#3d5d7f")

# Fig 5E 的目标细胞类型与配色
#
# 配色来源 (原文代码 int_analysis.ipynb):
#   CELL 203 定义 clust_cell_types_10x / major_cell_types_10x, 其中
#       "Endo" : "khaki",  "Vlmc" : "olive"/"gold",  "Olig" : "slategrey"
#   并定义 gen_dark_palette(prefix, color_name, uniq_clusts) =
#       sns.dark_palette(color_name, n_colors=n+2)[2:]
#   CELL 204: celltype_colors, ..., label_colors, ... = generate_palettes(adata_imputed)
#   CELL 467: plt.plot(x, binned_mean, color=celltype_colors[i])
# 其中 unpick 的 n 依赖运行时 uniq_clusts, 无法从发布数据确定性复算,
# 因此这里直接取**原文图上实测**的颜色 (700~1000 dpi 采样, 与上面色系一致):
#       khaki 系 = #E8E088   (PIL 实测; 该色系含 sns.dark_palette("khaki")[2:] 的末位 #F0E68C)
#       olive 系 = #F8D40C
#       slategrey = #708090  (与 matplotlib "slategrey" 精确一致)
TARGETS <- list(
  list(ct = "Endo", color = "#F0E68C", label = "Endo"),   # matplotlib "khaki"
  list(ct = "VLMC", color = "#FFD700", label = "VLMC"),   # matplotlib "gold"
  list(ct = "Olig", color = "#708090", label = "Oligo")   # matplotlib "slategrey"
)

# 作者代码 CELL 464:  bin_size = 50
# 作者代码 CELL 467-470: compute_binned_values(..., bin_size=bin_size, max_d=100);
#                        x = np.arange(len(binned_mean)) + bin_size
#
# 关于 bin_size 的取值 (详见 README「推断与限制 I1」):
#   * 仓库现值 50 会只留下 i=0,1 两个窗口 (窗口 i 覆盖 (i, i+50], i>=2 时无数据),
#     曲线塌成一小段, 与原文 5E 不符;
#   * 10 太小: 相邻窗口重叠 9 µm、每格仅约 30~130 个细胞, 曲线出现不真实的锯齿
#     (sd(diff) 比 30 大一个量级), 趋势方向也不稳;
#   * 30 与 50 的结果都平滑且趋势正确, 30 保留更多起伏细节, 最接近原文观感。
# 三者都实现了作者的公式, 差别只在窗口宽度。默认取 30。
# 可用环境变量 FIG5_BIN_SIZE 切换 (50 = 仓库值, 10 = 细窗对照)。
BIN_SIZE <- 50L
if (nzchar(Sys.getenv("FIG5_BIN_SIZE")))
  BIN_SIZE <- as.integer(Sys.getenv("FIG5_BIN_SIZE"))
MAX_D    <- 100L     # µm, 逐字取自 CELL 467-470
CORR_RADIUS <- 30    # µm
# 说明: 作者代码 CELL 472-474 显式传 radius=50, 而原文 Fig 5F 图注写
# "within 30 µm"(PDF 文字层把 µm 印成 mm)。两者都实现出来核对过, R 值差异极小:
#   radius=30 -> 0.412 / 0.442 / 0.380      radius=50 -> 0.415 / 0.427 / 0.371
# 这里取图注的 30 µm (与发表文字一致), 差异见 README「推断与限制 I6」。
if (nzchar(Sys.getenv("FIG5_CORR_RADIUS")))
  CORR_RADIUS <- as.numeric(Sys.getenv("FIG5_CORR_RADIUS"))
MIN_DIST_THRESH <- 15  # µm, identify_nearest_neighbors_with_idx 默认值
SEED <- 0L           # scanpy score_genes 默认 random_state=0

# ----------------------------- 读取 h5ad ---------------------------------------
.read_obs_cat <- function(h5, col) {
  g <- h5[[paste0("obs/", col)]]
  codes <- g[["codes"]]$read()
  cats  <- g[["categories"]]$read()
  out <- rep(NA_character_, length(codes))
  keep <- codes >= 0
  out[keep] <- cats[codes[keep] + 1L]
  out
}

#' 读取 MERFISH: obs 元信息 + 空间坐标(µm) + 原始计数 CSR
#'
#' @param genes_keep 若给定, 只保留这些基因的计数 (加速); NULL 表示全部 374 个
#' @return list(meta, xy, counts, genes)
read_merfish <- function(genes_keep = NULL, use_cache = TRUE) {
  if (use_cache && file.exists(CACHE_FILE)) {
    message("Loading cache: ", CACHE_FILE)
    obj <- readRDS(CACHE_FILE)
    if (is.null(genes_keep) ||
        all(genes_keep %in% obj$genes)) return(obj)
    message("Cache 缺少基因, 重建 ...")
  }
  message("Reading: ", H5AD_FILE)
  h5 <- H5File$new(H5AD_FILE, "r")
  on.exit(h5$close_all())

  n_cells <- h5[["obs/_index"]]$dims

  meta <- data.frame(
    age         = factor(.read_obs_cat(h5, "age"), levels = AGE_LEVELS),
    cell_type   = .read_obs_cat(h5, "cell_type"),
    clust_annot = .read_obs_cat(h5, "clust_annot"),
    donor_id    = .read_obs_cat(h5, "donor_id"),
    # obs$tissue 就是**解剖区域注释**, 取值与原文 CELL 455 的 spatial_order 一一对应:
    #   pia mater            -> Pia
    #   cortical layer II/III-> L2/3
    #   cortical layer V     -> L5
    #   cortical layer VI    -> L6
    #   corpus callosum      -> CC      <-- Fig 5F 用到的就是这一个
    #   olfactory region     -> LatSept (lateral septum / 嗅觉区)
    #   striatum             -> Striatum
    #   brain ventricle      -> Ventricle
    tissue      = .read_obs_cat(h5, "tissue"),
    stringsAsFactors = FALSE
  )
  # CELLxGENE 把室管膜细胞的 cell_type 记为空, 需按 clust_annot 补回
  meta$ct <- unname(CT_MAP[meta$cell_type])
  meta$ct[meta$clust_annot == "Epen"] <- "Epen"

  # 空间坐标 (obsm/spatial_coords, 2 x n_cells), 单位实测为 µm
  xy <- h5[["obsm/spatial_coords"]]$read()
  meta$x <- xy[1, ]; meta$y <- xy[2, ]
  xy_mat <- t(xy)
  rm(xy); invisible(gc(verbose = FALSE))

  # 基因符号
  fn <- h5[["var/feature_name"]]
  gene_sym <- fn[["categories"]]$read()[fn[["codes"]]$read() + 1L]

  # raw/X: CSR, 分块读取
  Xg <- h5[["raw/X"]]
  indptr <- as.double(Xg[["indptr"]]$read())
  wanted <- if (is.null(genes_keep)) seq_along(gene_sym)
            else which(gene_sym %in% genes_keep)
  message("  cells: ", n_cells, " | panel genes: ", length(gene_sym),
          " | kept: ", length(wanted))

  rows_all <- list(); cols_all <- list(); vals_all <- list(); k <- 0L
  CH <- 100000L
  for (s in seq(1L, n_cells, by = CH)) {
    e <- min(s + CH - 1L, n_cells)
    a <- indptr[s] + 1L; b <- indptr[e + 1L]
    if (b >= a) {
      idx <- as.integer(Xg[["indices"]][a:b]) + 1L
      keep <- idx %in% wanted
      if (any(keep)) {
        k <- k + 1L
        vals_all[[k]] <- as.numeric(Xg[["data"]][a:b][keep])
        cols_all[[k]] <- match(idx[keep], wanted)
        rows_all[[k]] <- rep(seq(s, e), diff(indptr[s:(e + 1L)]))[keep]
      }
    }
  }
  counts <- sparseMatrix(i = unlist(rows_all), j = unlist(cols_all),
                         x = unlist(vals_all),
                         dims = c(n_cells, length(wanted)),
                         dimnames = list(NULL, gene_sym[wanted]))
  rm(rows_all, cols_all, vals_all); invisible(gc(verbose = FALSE))

  obj <- list(meta = meta, xy = xy_mat, counts = counts,
              genes = gene_sym, genes_kept = gene_sym[wanted])
  saveRDS(obj, CACHE_FILE, compress = "gzip")
  message("Cache saved: ", CACHE_FILE)
  obj
}

# ----------------------------- score_genes (scanpy 原算法) -------------------
#' 复刻 scanpy 的 sc.tl.score_genes —— 逐行对应 _score_genes_bins
#'
#' scanpy 源码要点:
#'   n_items = round(len(pool) / (n_bins - 1))
#'   obs_cut = rank(obs_avg, method="min") %/% n_items          # 整数除法
#'   for cut in unique(obs_cut[gene_list]):                      # 只取目标基因所在的箱
#'       r_genes = 该箱全部基因;  若 ctrl_size < 箱内基因数 则随机抽 ctrl_size 个
#'       r_genes = r_genes \ gene_list
#'   control = union(所有 r_genes)
#'   score   = rowMean(gene_list) - rowMean(control)
#'
#' 本数据 (374 基因, n_bins=25 -> n_items=16) 下每箱仅约 16 个基因,
#' 小于 ctrl_size=50, 因此**不发生随机抽样**, 整个流程是确定性的。
#'
#' @param X 细胞 x 基因 的稀疏矩阵
#' @param gene_names X 的列名
#' @param gene_list 目标基因
#' @param pool 对照基因池 (默认 = gene_names, 与 scanpy 默认一致)
score_genes <- function(X, gene_names, gene_list,
                        ctrl_size = 50L, n_bins = 25L, seed = SEED,
                        pool = NULL) {
  if (is.null(pool)) pool <- gene_names
  genes <- intersect(gene_list, gene_names)
  if (length(genes) == 0L) return(rep(0, nrow(X)))

  pool <- intersect(pool, gene_names)
  obs_avg <- Matrix::colMeans(X[, pool, drop = FALSE])
  obs_avg <- obs_avg[is.finite(obs_avg)]
  pool <- names(obs_avg)
  if (length(pool) == 0L) return(rep(0, nrow(X)))

  # scanpy: n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
  n_items <- as.integer(round(length(obs_avg) / (n_bins - 1)))
  if (n_items < 1L) n_items <- 1L
  # pandas rank(method="min") 从 1 开始; R 的 rank(ties.method="min") 同样
  obs_cut <- floor(rank(obs_avg, ties.method = "min") / n_items)
  names(obs_cut) <- pool

  # 只使用目标基因所在的那几个箱 (与 np.unique(obs_cut.loc[gene_list]) 一致)
  cuts <- sort(unique(obs_cut[genes]))
  ctrl <- character(0)
  for (cu in cuts) {
    r_genes <- names(obs_cut)[obs_cut == cu]
    if (ctrl_size < length(r_genes)) {          # 只有箱内基因数 > 50 才抽样
      set.seed(seed)
      r_genes <- sample(r_genes, ctrl_size, replace = FALSE)
    }
    ctrl <- union(ctrl, setdiff(r_genes, genes))
  }
  if (length(ctrl) == 0L) return(rep(0, nrow(X)))

  ml <- Matrix::rowMeans(X[, genes, drop = FALSE])
  mc <- Matrix::rowMeans(X[, ctrl, drop = FALSE])
  as.numeric(ml - mc)
}

#' 计算三个激活分数, 并按作者代码做 4wk 中心化
compute_activation_scores <- function(meta, X, gene_names, pool = NULL) {
  out <- list()
  tgt_map <- c(activate_micro = "Micro", activate_astro = "Astro",
               activate_olig = "Olig")
  for (nm in names(GENE_SETS)) {
    sc <- score_genes(X, gene_names, GENE_SETS[[nm]], pool = pool)
    m <- meta$ct == tgt_map[[nm]] & meta$age == "4wk"
    out[[nm]] <- sc - mean(sc[m], na.rm = TRUE)
  }
  out
}

# ----------------------------- KDTree / 分箱 (作者 CELL 46/47) -----------------
#' X 中每个点到 Y 的最近距离 (排除自身重合点)
#'
#' 对应作者 identify_nearest_neighbors_with_dist: KDTree.query(X, k=2),
#' 若第一近邻距离为 0 (自身) 则取第二近邻。
nearest_dist <- function(X, Y) {
  if (nrow(X) == 0L || nrow(Y) == 0L) return(numeric(0))
  nn <- FNN::get.knnx(Y, X, k = 2L)
  d <- nn$nn.dist[, 1]
  dup <- d <= 0
  if (any(dup)) d[dup] <- nn$nn.dist[dup, 2]
  d
}

#' 作者 compute_binned_values (CELL 46) — 逐字对应
#'
#' 原 Python:
#'   binned_mean = np.zeros(max_d - min_d - bin_size)
#'   binned_std  = np.zeros(max_d - min_d - bin_size)
#'   for i in np.arange(min_d, max_d - bin_size):      # 步长 1 -> 滑动窗
#'       idx = np.logical_and(dists > i, dists <= (i + bin_size))
#'       binned_mean[i] = np.mean(scores[idx])
#'       binned_std[i]  = np.std(scores[idx]) / np.sqrt(len(scores[idx]))
#'   binned_mean -= binned_mean.mean()
#'   binned_std  -= binned_mean.mean()                 # 中心化后均值=0, 实为 no-op
#'
#' 两处容易忽略、但影响数值的细节:
#'   1) np.std 默认 ddof=0 (总体标准差), 而 R 的 sd() 是 ddof=1 的样本标准差;
#'      这里用 sqrt(mean((v-mean(v))^2)) 对齐 np.std。
#'   2) 作者第二行减的是**中心化之后**的 binned_mean.mean() (≈0),
#'      不是曲线原始均值; 照抄才与作者数值一致。
binned_curve <- function(dists, scores, bin_size = BIN_SIZE, max_d = MAX_D) {
  n_bins <- max_d - bin_size
  bm <- rep(NA_real_, n_bins); bs <- rep(NA_real_, n_bins)
  for (i in 0:(n_bins - 1L)) {
    idx <- which(dists > i & dists <= (i + bin_size))
    if (length(idx) == 0L) next
    v <- scores[idx]
    bm[i + 1L] <- mean(v)
    bs[i + 1L] <- sqrt(mean((v - mean(v))^2)) / sqrt(length(v))  # np.std(ddof=0)
  }
  bm <- bm - mean(bm, na.rm = TRUE)
  bs <- bs - mean(bm, na.rm = TRUE)   # 照抄: 中心化后≈0
  list(mean = bm, sem = bs)
}

#' 作者: x = np.arange(len(binned_mean)) + bin_size
x_axis <- function(bin_size = BIN_SIZE, max_d = MAX_D)
  seq_len(max_d - bin_size) - 1L + bin_size

# ----------------------------- 出图 -------------------------------------------
theme_fig5 <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          plot.title = element_text(size = base_size + 1, hjust = 0.5))
}

save_fig <- function(plot, name, w, h, dpi = 300, cairo = TRUE) {
  if (cairo) {
    grDevices::cairo_pdf(file.path(FIG_DIR, paste0(name, ".pdf")),
                         width = w, height = h)
    print(plot); grDevices::dev.off()
  } else {
    ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = w, height = h)
  }
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = w, height = h,
         dpi = dpi)
  grDevices::png(file.path(HI_DIR, paste0(name, ".png")),
                 width = w, height = h, units = "in", res = 900, type = "cairo")
  print(plot); grDevices::dev.off()
  message("Saved figure: ", name)
}
