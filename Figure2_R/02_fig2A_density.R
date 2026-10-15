# ==============================================================================
# 02_fig2A_density.R — Figure 2A: 主要细胞类型密度随年龄变化
#
# 【本版修正】与作者原始代码逐条对齐
#   (github.com/ZhuangLab/SpatialBrainAgingCell22, analyze_merfish_celltyping.ipynb)
#
#   原文 CELL 308/309:
#       all_ct = ["ExN","InN","MSN","Olig","Astro","Endo","Micro","OPC","Peri",
#                 "Vlmc","Epen","Macro","T cell"]
#       for i in batches:                       # batch = 一个动物
#           curr_batch_area = 0
#           for j in curr_batch.obs.slice.unique():
#               hull = ConvexHull(curr_slice.obsm['spatial'])
#               curr_batch_area += hull.area    # scipy 在 2D 返回的是**周长**
#           counts.append(1000 * n_cells(ct) / curr_batch_area)
#
#   修正点 (相对于本项目上一版):
#     1. 密度分母: 由"真实凸包面积"改为原文的 **sum(ConvexHull.area)**
#        —— 2D 下 scipy 的 .area 是**周长**, 再乘 1000 (坐标单位 nm)。
#        旧口径把 ExN 从 331.5 抬到 373.9, 并改变了 Olig 的年龄趋势。
#     2. 误差线: 由 t 分布 CI 改为 **bootstrap 95% CI** (n_boot=1000),
#        对齐 seaborn barplot 的默认行为。旧口径使 InN 等假显著。
#     3. 配色: 改用作者代码里的 age_colors = cornflowerblue/thistle/lightcoral。
#     4. Epen: 改用 cell_type 字段再做映射, 使 Epen 细胞不再被丢弃
#        (CELLxGENE 把室管膜细胞的 cell_type 记为空, 需按 clust_annot 补回)。
#
#   统计: 只做 4wk vs 90wk 两独立样本 t 检验 (var.equal = TRUE),
#         BH 校正在**全部可检验的细胞类型**上做。
#   y 轴: 刻度 0-100-200-300-400-500 (与原文一致), 上限放宽以容纳误差棒。
#
# 输出: figures/fig2A_density.pdf/png, tables/fig2A_*.csv
# ==============================================================================
if (!exists("PROJECT_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure2_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

obj  <- load_data()
meta <- obj$meta
meta <- meta[!is.na(meta$cell_type), ]      # 丢掉仍未映射的杂类

# ---- 1. 每个 (供体, 切片) 的凸包 --------------------------------------------
# 原文是"每个 batch 的每个 slice 各取一次凸包, 再相加"。
# 空间坐标为 nm; scipy 的 ConvexHull.area 在 2D = 周长。
# 本函数同时给出 poly_perimeter (原文口径) 与 poly_area (真实面积), 便于对比。
message("Computing convex hull per (donor, slice) ...")
hulls <- meta %>%
  group_by(donor_id, slice) %>%
  group_modify(function(df, key) {
    xy <- unique(df[, c("x", "y")])
    if (nrow(xy) < 4) {
      return(data.frame(hull_area = NA_real_, hull_perimeter = NA_real_,
                        true_area_um2 = NA_real_))
    }
    idx <- grDevices::chull(xy$x, xy$y)
    data.frame(
      hull_area      = poly_area(xy$x[idx], xy$y[idx]),        # 2D 即面积
      hull_perimeter = poly_perimeter(xy$x[idx], xy$y[idx]),   # 2D 即周长
      true_area_um2  = poly_area(xy$x[idx], xy$y[idx]) / 1e6   # nm² -> µm²
    )
  }) %>% ungroup() %>% filter(!is.na(hull_area))

donor_area <- hulls %>%
  group_by(donor_id) %>%
  summarise(paper_area   = sum(hull_perimeter),   # 原文口径: 周长之和
            real_area_mm2 = sum(true_area_um2) / 1e6,
            .groups = "drop")

# ---- 2. 供体水平密度 ---------------------------------------------------------
n_ct <- meta %>% count(donor_id, cell_type, age, name = "cell_count")

# 注意: age 已经在 n_ct 里, 不要再 join 一次 (会产生 age.x / age.y 两列,
# 后面 mutate(age = ...) 就会报 "object 'age' not found")。
df_donor <- expand.grid(donor_id  = unique(meta$donor_id),
                        cell_type = CELLTYPE_ORDER,
                        stringsAsFactors = FALSE) %>%
  left_join(n_ct, by = c("donor_id", "cell_type")) %>%
  replace_na(list(cell_count = 0L)) %>%
  left_join(donor_area, by = "donor_id") %>%
  # 每个 (供体, 细胞类型) 补上该供体的年龄; 供体年龄唯一, 用 first 即可
  group_by(donor_id) %>%
  mutate(age = age[which(!is.na(age))[1]]) %>%
  ungroup() %>%
  mutate(
    # 原文公式: counts = 1000 * n_cells / sum(ConvexHull.area over slices)
    # 坐标单位 nm; 该式等价于把 nm 当 µm 处理。
    counts = 1000 * cell_count / paper_area,
    cell_type = factor(cell_type, levels = CELLTYPE_ORDER),
    age = factor(age, levels = AGE_LEVELS)
  ) %>%
  select(donor_id, cell_type, age, cell_count, paper_area, real_area_mm2, counts)

# 供体水平自检: 每组必须是"个位数"只动物
n_grp <- df_donor %>% count(cell_type, age)
message("每个 cell_type x age 的供体数 (应是个位数): ",
        paste(range(n_grp$n), collapse = " - "))
stopifnot(max(n_grp$n) < 20)

write.csv(df_donor, file.path(TAB_DIR, "fig2A_donor_density.csv"),
          row.names = FALSE)

# ---- 3. 分组统计: 均值 + bootstrap 95% CI (对齐 sns.barplot) -----------------
summary_df <- df_donor %>%
  group_by(cell_type, age) %>%
  group_modify(function(df, key) {
    v <- df$counts
    ci <- boot_ci(v)
    data.frame(mean_density = mean(v, na.rm = TRUE),
               n_donors = sum(is.finite(v)),
               ci_low = ci["low"], ci_high = ci["high"])
  }) %>% ungroup() %>%
  mutate(age = factor(age, levels = AGE_LEVELS))

message("\nFig 2A 每 cell_type x age 的供体均值 (原文口径 1000*n/周长):")
print(summary_df %>%
        select(cell_type, age, mean_density, n_donors) %>%
        pivot_wider(names_from = age, values_from = mean_density) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
        as.data.frame())

# ---- 4. 显著性: 只做 4wk vs 90wk, BH 校正 ------------------------------------
pvals <- vapply(CELLTYPE_ORDER, function(ct) {
  v1 <- df_donor$counts[df_donor$cell_type == ct & df_donor$age == "4wk"]
  v2 <- df_donor$counts[df_donor$cell_type == ct & df_donor$age == "90wk"]
  if (length(v1) < 2 || length(v2) < 2) NA_real_ else ttest_p(v1, v2)
}, numeric(1))

# BH 校正只在有效 p 值上做 (NA 不能参与, 否则整组变 NA)
valid <- !is.na(pvals)
fdr <- rep(NA_real_, length(pvals))
fdr[valid] <- p.adjust(pvals[valid], method = "BH")

sig_tab <- data.frame(cell_type = CELLTYPE_ORDER, p_4wk_vs_90wk = pvals,
                      fdr_BH = fdr, significant = !is.na(fdr) & fdr < 0.05)
message("\nFig 2A 显著性 (仅 4wk vs 90wk, BH 校正范围 = ",
        sum(valid), " 个可检验类型):")
print(sig_tab, digits = 3, row.names = FALSE)
write.csv(summary_df, file.path(TAB_DIR, "fig2A_density_summary.csv"),
          row.names = FALSE)
write.csv(sig_tab, file.path(TAB_DIR, "fig2A_significance.csv"),
          row.names = FALSE)

# ---- 5. 主图 -----------------------------------------------------------------
# y 轴刻度按原文只标到 500; 上限放宽到 620 以容纳 ExN 的误差棒上帽。
Y_BREAKS <- seq(0, 500, 100)
Y_LIMITS <- c(0, 620)
DODGE_W  <- 0.7        # 与原文 sns.barplot 的分组宽度一致
BAR_W    <- 0.65
ERR_W    <- 0.25

theme_fig2a <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid   = element_blank(),   # 原文 sns.despine() 无网格
          panel.border = element_blank(),
          axis.line    = element_line(linewidth = 0.3),
          axis.text.x  = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y  = element_text(colour = "black"),
          legend.position = "top",
          legend.title = element_blank(),
          strip.background = element_blank())
}

#' 分组柱状图 + bootstrap 95% CI 误差线 (+ 可选显著性横线)
p_bars <- function(d, y_breaks, y_limits, show_legend = FALSE,
                   base_size = 11, drop_unused_x = FALSE) {
  ggplot(d, aes(x = cell_type, y = mean_density, fill = age)) +
    geom_col(position = position_dodge(width = DODGE_W), width = BAR_W,
             colour = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  position = position_dodge(width = DODGE_W),
                  width = ERR_W, linewidth = 0.45) +
    scale_fill_manual(values = AGE_COLORS,
                      labels = AGE_LABELS_SHORT[AGE_LEVELS], drop = FALSE) +
    scale_x_discrete(drop = drop_unused_x) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits,
                       expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = expression("Cells / mm"^2), fill = NULL) +
    theme_fig2a(base_size) +
    theme(legend.position = if (show_legend) "top" else "none")
}

p_main <- p_bars(summary_df, Y_BREAKS, Y_LIMITS, show_legend = TRUE)

# ---- 6. 显著性横线: 按 dodge 后的实际 x 坐标手动指定 --------------------------
X_OFF <- DODGE_W / 3                                # 3 组时的偏移量
sig_plot_df <- sig_tab %>%
  filter(significant) %>%
  mutate(cell_type = factor(cell_type, levels = CELLTYPE_ORDER),
         xmin = as.numeric(cell_type) - X_OFF,
         xmax = as.numeric(cell_type) + X_OFF) %>%
  left_join(summary_df %>%
              filter(age %in% c("4wk", "90wk")) %>%
              group_by(cell_type) %>%
              summarise(y_max = max(ci_high, na.rm = TRUE), .groups = "drop"),
            by = "cell_type") %>%
  mutate(y_position = y_max + 0.04 * Y_LIMITS[2])

p_main <- p_main +
  geom_segment(data = sig_plot_df,
               aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
               inherit.aes = FALSE, linewidth = 0.45) +
  geom_text(data = sig_plot_df,
            aes(x = (xmin + xmax) / 2, y = y_position, label = "*"),
            inherit.aes = FALSE, vjust = -0.35, size = 4.2)

# ---- 7. 右上角单个 inset: 4 类低丰度细胞 -------------------------------------
inset_df <- summary_df %>%
  filter(cell_type %in% CELLTYPE_INSET) %>%
  mutate(cell_type = factor(cell_type, levels = CELLTYPE_INSET))

p_inset <- p_bars(inset_df, waiver(), c(0, NA), show_legend = FALSE,
                  base_size = 9, drop_unused_x = TRUE) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(4, 6, 4, 4))

pA <- p_main +
  inset_element(p_inset, left = 0.60, bottom = 0.55, right = 1.00, top = 1.00,
                align_to = "panel", clip = FALSE) +
  plot_annotation(
    title = "Fig 2A (reproduced): Cell-type density across ages",
    subtitle = paste0("Donor-level density, 原文口径; bars = mean, ",
                      "error bars = bootstrap 95% CI; * = 4wk vs 90wk, BH FDR < 0.05."))

save_fig(pA, "fig2A_density", 12, 6, cairo = TRUE)
pA_main <- pA
message("Fig 2A 完成 ✔")
