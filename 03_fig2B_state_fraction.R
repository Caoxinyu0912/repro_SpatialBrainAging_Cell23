# ==============================================================================
# 03_fig2B_state_fraction.R — Figure 2B: Olig/Endo/Micro/Astro 各细胞状态
#                             组成的 100% 堆叠柱状图随年龄变化
#
# 版面 (按原文 2 x 2 排列):
#   - 左上 Olig | 右上 Endo | 左下 Micro | 右下 Astro (从左到右、从上到下);
#   - 每个大类一个分面, 分面内 3 根竖柱: Juvenile / Young adult / Old;
#   - 各状态纵向堆叠至 100%, 同大类内深浅渐变配色;
#   - 状态名直接以对应颜色文字标在各分面内堆叠柱左侧 (无单独图例);
#   - 年龄组注释 (Juvenile / Young adult / Old, 竖排文字) 仅出现在第一行
#     两个分面的柱子上方;
#   - 图中不标注显著性星号 (统计表仍输出到 tables/)。
#
# 统计: 状态比例先按供体计算, 跨供体均值 ± 95% CI (输出 tables/);
#       Juvenile vs Old t 检验 + BH FDR (输出 tables/)。
# 输出: figures/fig2B_*.pdf/png, tables/fig2B_*.csv
# ==============================================================================
if (!exists("PROJECT_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure2_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
}
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

obj <- load_data()
meta <- obj$meta

# 2x2 排列顺序: Olig -> Endo -> Micro -> Astro (从左到右、从上到下)
PARENT_2x2 <- c("Olig", "Endo", "Micro", "Astro")

state_map <- bind_rows(lapply(names(STATES_OF), function(p)
  data.frame(parent = p, state = STATES_OF[[p]])))

# ---- 1. 每供体: 各状态细胞数 / 大类总数 -> 状态比例 ---------------------------
meta_states <- meta %>%
  filter(clust_annot %in% ALL_STATES) %>%
  left_join(state_map, by = c("clust_annot" = "state"))

frac <- meta_states %>%
  count(donor_id, age, parent, clust_annot, name = "n_state") %>%
  complete(donor_id, nesting(age, parent, clust_annot), fill = list(n_state = 0L)) %>%
  group_by(donor_id, age, parent) %>%
  mutate(n_parent = sum(n_state)) %>%
  ungroup() %>%
  mutate(fraction = if_else(n_parent > 0, n_state / n_parent, NA_real_),
         parent = factor(parent, levels = PARENT_2x2),
         clust_annot = factor(clust_annot, levels = ALL_STATES),
         age = factor(age, levels = AGE_LEVELS))

# ---- 2. 统计: 均值 ± 95%CI (输出表); Juvenile vs Old t 检验 + BH FDR ----------
stat <- frac %>%
  filter(!is.na(fraction)) %>%
  group_by(parent, clust_annot, age) %>%
  summarise(mean = mean(fraction), sd = sd(fraction), n = n(), .groups = "drop") %>%
  mutate(ci_half = ifelse(n > 1, qt(0.975, pmax(n - 1, 1)) * sd / sqrt(n), 0),
         lower = mean - ci_half, upper = mean + ci_half,
         age = factor(age, levels = AGE_LEVELS))

pvals <- sapply(ALL_STATES, function(st) {
  v1 <- frac$fraction[frac$clust_annot == st & frac$age == "4wk"]
  v2 <- frac$fraction[frac$clust_annot == st & frac$age == "90wk"]
  ttest_p(v1, v2)
})
fdr <- p.adjust(pvals, method = "BH")
sig_tab <- data.frame(state = ALL_STATES, p_4wk_vs_90wk = pvals, fdr_BH = fdr,
                      significant = fdr < 0.05)
message("\nFig 2B 统计检验 (Juvenile vs Old, t-test + BH FDR) — 仅输出表格, 图中不标星号:")
print(sig_tab, digits = 3)

write.csv(stat %>% arrange(parent, clust_annot, age),
          file.path(TAB_DIR, "fig2B_state_fraction_summary.csv"), row.names = FALSE)
write.csv(sig_tab, file.path(TAB_DIR, "fig2B_significance.csv"), row.names = FALSE)

# ---- 3. 100% 堆叠柱状图 (2x2, 无图例, 状态名标注在柱左侧) ----------------------
stack <- stat %>%
  select(parent, clust_annot, age, mean) %>%
  group_by(parent, age) %>%
  arrange(clust_annot, .by_group = TRUE) %>%        # 状态 1,2,3 自下而上堆叠
  mutate(cum_top = cumsum(mean), cum_bot = lag(cum_top, default = 0),
         mid = (cum_top + cum_bot) / 2) %>%
  ungroup()

# 状态名文字: 每个分面内、堆叠柱左侧, 颜色 = 状态色 (与原文一致)
state_lab <- stack %>%
  filter(age == "4wk") %>%
  distinct(parent, clust_annot, mid) %>%
  mutate(x = 0.38)

# 年龄注释文字: 竖排, 位于柱顶上方 (仅第一行分面使用)
age_lab <- data.frame(age = factor(AGE_LEVELS, levels = AGE_LEVELS),
                      x = seq_along(AGE_LEVELS), y = 1.06,
                      label = AGE_LABELS_SHORT[AGE_LEVELS])

# 0% / 100% 标度: 不画在坐标轴上, 而是以文字放在**该行第一列分面的右缘**
# (也就是第一列与第二列之间的空隙左侧)。
#   三根柱中心在 x = 1,2,3, 柱宽 0.58 => 最右柱右缘 3.29; scale 范围 0.30-3.70
#   加 5% expansion 后右边界约 3.54。x 必须 <= limits 上界, 否则该行会被
#   scale 丢弃 (这正是之前标签不显示的原因), 所以放在 3.56、左对齐、贴右缘。
SCALE_LAB_X <- 3.56
scale_lab_df <- function(parent_first)
  data.frame(parent = factor(parent_first, levels = PARENT_2x2),
             x = SCALE_LAB_X, y = c(0, 1), label = c("0%", "100%"))

p_bar2 <- function(parents_show, show_age_lab = TRUE, base_size = 10) {
  d <- stack %>% filter(parent %in% parents_show) %>%
    mutate(age_x = as.numeric(age))
  labs_st <- state_lab %>% filter(parent %in% parents_show)
  # 0%/100% 只画在该行第一个 (左) 分面
  first_parent <- parents_show[1]
  sl <- scale_lab_df(first_parent)
  p <- ggplot(d, aes(x = age_x, y = mean, fill = clust_annot, group = clust_annot)) +
    geom_col(width = 0.58, color = "white", linewidth = 0.22,
             position = position_stack(reverse = TRUE)) +   # state-1 在底部
    geom_text(data = labs_st, aes(x = x, y = mid, label = clust_annot, color = clust_annot),
              inherit.aes = FALSE, size = 2.6, hjust = 0.5) +
    # 0% / 100% 文字 (只出现在该行左分面的右缘, 即两列之间的空隙处)
    geom_text(data = sl, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 2.9, color = "black") +
    scale_color_manual(values = STATE_COLORS, guide = "none") +
    scale_fill_manual(values = STATE_COLORS, guide = "none") +
    facet_wrap(~ parent, ncol = 2) +
    scale_x_continuous(breaks = seq_along(AGE_LEVELS),
                       labels = rep("", 3),          # 年龄只标在柱顶 (竖排), 底部不重复
                       limits = c(0.30, 3.70), expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, labels = NULL,
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),            # 大类名靠左侧彩色状态名标明
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # 去掉所有坐标轴线与刻度 (标度改用上面的文字)
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          # 顶部留白: 竖排的年龄标签会伸出面板上沿, 不留空间会被 plot 标题/说明压住
          plot.margin = margin(14, 4, 2, 2),
          panel.spacing = unit(2.4, "lines"))
  if (show_age_lab) {
    p <- p + geom_text(data = age_lab,
                       aes(x = x, y = y, label = label),
                       inherit.aes = FALSE, angle = 90, hjust = 0, size = 2.6) +
      scale_y_continuous(breaks = NULL, labels = NULL,
                         expand = expansion(mult = c(0, 0.22)))
  }
  p
}

pB_top <- p_bar2(c("Olig", "Endo"), show_age_lab = TRUE)    # 第一行: Olig | Endo
pB_bot <- p_bar2(c("Micro", "Astro"), show_age_lab = FALSE) # 第二行: Micro | Astro

pB <- pB_top / pB_bot +
  # 与 Fig 2A / 2C 统一的标题 + 说明 (三个 panel 版面格式保持一致)
  # 注意: 图中不标显著性星号, 因此说明里也不写 "* = ..." 的图例,
  #       显著性结果只输出到 tables/fig2B_significance.csv。
  plot_annotation(
    title = "Fig 2B (reproduced): Cell-state composition across ages",
    subtitle = paste0("Mean state fraction per donor, stacked to 100%.  ",
                      "State colours shared with Fig 2C."),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 7, color = "grey25",
                                               lineheight = 1.15)))

save_fig(pB, "fig2B_state_fraction", 8.6, 7.4)
message("Fig 2B 完成 ✔")
