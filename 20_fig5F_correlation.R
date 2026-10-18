# ==============================================================================
# 20_fig5F_correlation.R — Figure 5F
#
# 3 个子图: 二维核密度填充图 (对应原文 sns.kdeplot(x=..., y=..., fill=True))
#   左  : (Oligo inflammation score, Micro activation score)
#   中  : (Oligo inflammation score, Astro activation score)
#   右  : (Micro activation score, Astro activation score)
#
# 忠实于原文的定义表达 (作者 int_analysis.ipynb CELL 46 + 472-474):
#
#   # CELL 472-474
#   x, y = compute_celltype_obs_correlation(
#       adata_imputed[adata_imputed.obs.spatial_clust_annots == "CC"],
#       src, tgt, f"activate_{src}", f"activate_{tgt}", radius=50)
#   R = np.corrcoef(x, y)[0, 1]
#   sns.kdeplot(x=x, y=y, fill=True)
#
#   # CELL 46 —— 两个关键细节
#   def identify_nearest_neighbors_with_idx(X, Y, dist_thresh, min_dist_thresh=15):
#       ind, dists = KDTree(Y).query_radius(X, r=dist_thresh, return_distance=True)
#       ...
#       ind = ind[dists > min_dist_thresh]          # 剔除 <= 15 µm 的邻居
#       return ind, ind_X
#   def compute_celltype_obs_correlation(...):
#       neighbors_X, ind_X = identify_nearest_neighbors_with_idx(
#           curr_X, curr_Y, dist_thresh=radius, min_dist_thresh=min_dist_thresh)
#       curr_expr = obs_Y[neighbors_X]              # 逐邻居取值 (非均值)
#       return obs_X.values[ind_X], curr_expr.values
#
# 区域来源 (重要):
#   作者用 obs.spatial_clust_annots == "CC"。CELLxGENE 发布版**没有**这个字段,
#   但发布了等价的 obs$tissue, 取值与原文 CELL 455 的 spatial_order
#   ['Pia','L2/3','L5','L6','CC','LatSept','Striatum','Ventricle'] 一一对应,
#   其中 "corpus callosum" 即 CC (53,460 细胞, 少突胶质占 66.7%, 是典型白质束)。
#   因此这里**直接使用发布版的区域注释**, 不再自行近似。
# ==============================================================================
if (!exists("FIG5_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure5_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
}

# ---- 1. 数据与激活分数 --------------------------------------------------------
obj <- read_merfish()
meta <- obj$meta
keep <- !is.na(meta$ct)
meta <- meta[keep, ]; xy <- obj$xy[keep, ]
X <- obj$counts[keep, , drop = FALSE]
act <- compute_activation_scores(meta, X, colnames(X))
for (nm in names(act)) meta[[nm]] <- act[[nm]]
message("cells: ", nrow(meta))

# ---- 2. 胼胝体 (CC) —— 直接用发布版的区域注释 --------------------------------
in_cc <- meta$tissue == CC_TISSUE
message(sprintf("CC 区域 (obs$tissue == '%s'): %d 细胞 (%.2f%%), 其中 Olig %.1f%%",
                CC_TISSUE, sum(in_cc), 100 * mean(in_cc),
                100 * mean(meta$ct[in_cc] == "Olig")))
for (ctc in c("Micro", "Astro", "Olig")) {
  message(sprintf("   %-6s n=%6d", ctc, sum(in_cc & meta$ct == ctc)))
}

# ---- 3. 配对 (作者 compute_celltype_obs_correlation 的定义) -------------------
#' 逐邻居配对: 区域内每个 src 细胞 vs 其 radius 内 (且距离 > min_dist) 的 tgt 细胞
pair_neighbors <- function(src, tgt, mask_src, radius = CORR_RADIUS,
                           min_dist = MIN_DIST_THRESH) {
  msrc <- which(meta$ct == src & mask_src)
  mtgt <- which(meta$ct == tgt)
  if (length(msrc) == 0L || length(mtgt) == 0L)
    return(list(x = numeric(0), y = numeric(0)))
  Xs <- xy[msrc, , drop = FALSE]
  Ys <- xy[mtgt, , drop = FALSE]
  # 半径内邻居数上限 (防止个别细胞邻居过多拖慢循环)
  kmax <- min(nrow(Ys), 400L)
  nn <- FNN::get.knnx(Ys, Xs, k = kmax)
  sx <- meta[[paste0("activate_", tolower(src))]][msrc]
  ty <- meta[[paste0("activate_", tolower(tgt))]][mtgt]
  # 性能优化: 用逻辑下标一次取出所有合法邻居对, 替代逐行 for 循环。
  # nn$nn.index[sel] 与 nn$nn.dist[sel] 的元素一一对应 (同为列优先展开),
  # 因此 x 与 y 天然配对; 只需再算出每个命中元素所属的 src 行号。
  sel <- nn$nn.dist <= radius & nn$nn.dist > min_dist
  if (!any(sel)) return(list(x = numeric(0), y = numeric(0)))
  pos <- which(sel)                                  # 列优先的线性下标
  src_row <- ((pos - 1L) %% nrow(Xs)) + 1L           # 对应的 src 细胞行号
  list(x = sx[src_row], y = ty[nn$nn.index[sel]])
}

PAIRS <- list(
  list(src = "Olig",  tgt = "Micro", xlab = "Oligo inflammation score",
       ylab = "Micro activation score", R_paper = 0.46,
       xlim = c(-2, 4.5), ylim = c(-1.8, 3.5)),
  list(src = "Olig",  tgt = "Astro", xlab = "Oligo inflammation score",
       ylab = "Astro activation score", R_paper = 0.47,
       xlim = c(-2, 5.5), ylim = c(-1.8, 3.5)),
  list(src = "Micro", tgt = "Astro", xlab = "Micro activation score",
       ylab = "Astro activation score", R_paper = 0.46,
       xlim = c(-2, 3.5), ylim = c(-1.8, 3.5))
)

res <- list(); pdf_dat <- list()
for (i in seq_along(PAIRS)) {
  p <- PAIRS[[i]]
  pr <- pair_neighbors(p$src, p$tgt, in_cc)
  R <- if (length(pr$x) > 30L) cor(pr$x, pr$y) else NA_real_
  message(sprintf("  %-6s -> %-6s  R = %+0.3f  (n=%d, 原文 R=%.2f)",
                  p$src, p$tgt, R, length(pr$x), p$R_paper))
  res[[i]] <- data.frame(source = p$src, target = p$tgt, R = R,
                         n_pairs = length(pr$x), R_paper = p$R_paper,
                         radius_um = CORR_RADIUS, region = CC_TISSUE,
                         n_cells_region = sum(in_cc),
                         min_dist = MIN_DIST_THRESH)
  pdf_dat[[i]] <- data.frame(panel = i, source = p$src, target = p$tgt,
                             x = pr$x, y = pr$y)
}
res_df <- do.call(rbind, res)
print(res_df, digits = 3, row.names = FALSE)
write.csv(res_df, file.path(TAB_DIR, "fig5F_correlation_R.csv"),
          row.names = FALSE)
# 只保存摘要统计, 避免逐配对文件过大
summ <- do.call(rbind, lapply(seq_along(pdf_dat), function(i) {
  d <- pdf_dat[[i]]
  data.frame(panel = i, source = PAIRS[[i]]$src, target = PAIRS[[i]]$tgt,
             n_pairs = nrow(d), R = res_df$R[i],
             x_mean = mean(d$x), x_sd = sd(d$x),
             y_mean = mean(d$y), y_sd = sd(d$y))
}))
write.csv(summ, file.path(TAB_DIR, "fig5F_pairs_summary_R.csv"),
          row.names = FALSE)

# ---- 4. 绘图 (二维核密度填充, 与原文 sns.kdeplot(fill=True) 对应) -------------
pdat <- do.call(rbind, pdf_dat)

# 填充色:
# 原文 Fig 5F 只写了 sns.kdeplot(x=x, y=y, fill=True), 没指定 cmap,
# 所以用的是 seaborn 的默认 2D-KDE 填充梯度 (蓝系, 由浅到深多级)。
# KDE_FILL_COLORS 直接取自该 cmap 的 10 个等距代表色 (见 00_common.R)。
#
# 让外观更贴近原文的四个要点:
#   1) **平滑度**: seaborn 的 scipy gaussian_kde 带宽比 ggplot 默认更大,
#      原文轮廓因此非常圆润; 这里显式放大 h (h = c(bw_x, bw_y)) 对齐。
#   2) **色阶**: breaks 用幂次间隔, 大部分层级落在低密度段, 等高线收得紧;
#   3) **留白**: 最外一层填白色, 相当于把极稀疏外围擦掉;
#   4) **坐标范围**: 沿用各面板的 xlim/ylim (与原文一致)。
kd_pow  <- 5.0
bw_mult <- 1.15     # 相对 ggplot 默认带宽的放大倍数 (对齐 seaborn 的圆润轮廓)

mkF <- function(i) {
  d <- subset(pdat, panel == i)
  n_lv <- length(KDE_FILL_COLORS) + 1L
  brk <- seq(0, 1, length.out = n_lv)^kd_pow
  fill_vals <- c("white", KDE_FILL_COLORS)   # 最外层白 = 留白
  hx <- bw_mult * sd(d$x)
  hy <- bw_mult * sd(d$y)
  ggplot(d, aes(x = x, y = y)) +
    geom_density_2d_filled(
      contour_var = "ndensity",
      breaks = brk,
      h = c(hx, hy),
      show.legend = FALSE
    ) +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    coord_cartesian(xlim = PAIRS[[i]]$xlim, ylim = PAIRS[[i]]$ylim) +
    annotate("text", x = PAIRS[[i]]$xlim[1], y = PAIRS[[i]]$ylim[2],
             label = sprintf("italic(R) == %s",
                             formatC(res_df$R[i], format = "f", digits = 2)),
             hjust = -0.12, vjust = 1.6, parse = TRUE, size = 3.8) +
    labs(x = PAIRS[[i]]$xlab, y = PAIRS[[i]]$ylab) +
    theme_fig5(base_size = 9)
}

library(cowplot)
pF_body <- plot_grid(mkF(1), mkF(2), mkF(3), nrow = 1, align = "h",
                     axis = "tb")
ttl <- ggdraw() +
  draw_label("Figure 5F (reproduced in R): Correlation of activation / inflammation scores in corpus callosum",
             x = 0.008, y = 0.80, hjust = 0, vjust = 0.5,
             fontface = "bold", size = 11, colour = "black") +
  draw_label(sprintf(paste0("Region: obs$tissue == \"corpus callosum\" ",
                            "(%.0f k cells, 66.7%% oligodendrocytes); ",
                            "per-neighbour pairing within %g um ",
                            "(min_dist > %g um)"),
                     sum(in_cc) / 1000, CORR_RADIUS, MIN_DIST_THRESH),
             x = 0.008, y = 0.24, hjust = 0, vjust = 0.5,
             size = 7.5, colour = "grey25")

pF <- plot_grid(ttl, pF_body, ncol = 1, rel_heights = c(0.10, 1))
save_fig(pF, "fig5F_correlation_R", 9.5, 3.7)
message("Fig 5F (R) 完成 ✔")
