# ==============================================================================
# 06_alt_2d_kde_compare.R — 备选方案对比: 二维 KDE 轮廓 vs 一维 KDE 小提琴
#
# 目的 (供报告"失败环节分析 / 改进建议"引用):
#   原文 Figure 2C 的形状来自 seaborn 的二维核密度填充
#   (`sns.kdeplot(x=expr, y=jitter, fill=True)`)。本项目最初也尝试在 R 中
#   忠实复刻这一画法: 对 (表达量, 抖动坐标) 做二维高斯核密度估计, 再取最外层
#   填充等高线。结论是**在本数据上不可行**, 原因如下:
#
#   MERFISH 表达矩阵极度稀疏 —— 例如 Olig-3 中 C4b 只有约 4% 细胞为阳性,
#   其余 96% 的表达值为 0。这 96% 构成一个 "表达量 = 0" 的**点质量**, 它在水平
#   抖动方向上是近似均匀分布的, 因此二维密度场是一个"沿 y 方向近乎平顶的长条"。
#   任何能够包含足够概率质量的等高线都会被拉成横跨整个抖动区间的细长条,
#   无法得到原文那种紧凑的纺锤/菱形 (原文用的是插补后更密集的表达矩阵,
#   见 README 差异分析第 6 条)。若要强行得到紧凑形状, 只能把等高线水平抬到
#   很高, 那样又会把形状压缩成一根看不见的细线。
#
#   因此最终采用一维 KDE (violin_outline(), 见 00_common.R): 形状与原文一致,
#   但完全不受零膨胀点质量影响, 且不需要手动抖动。
#
# 本脚本输出对比图, 作为上述分析的证据留档, 不参与 Figure 2 的最终拼版。
# 输出: figures/fig2C_alt_2d_kde_vs_1d.png
# ==============================================================================
if (!exists("PROJECT_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure2_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
}
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

set.seed(20230105)
obj <- load_data()
meta <- obj$meta; counts <- obj$counts

p <- "Olig"
genes <- intersect(GENES_TO_SHOW[[p]], c("C4b", "Neat1", "Il33"))  # panel 内, 可对比
states <- STATES_OF[[p]]
rows <- which(meta$clust_annot %in% states)

one_cell <- function(g, st, method) {
  sel <- rows[meta$clust_annot[rows] == st]
  v <- sqrt(as.numeric(counts[sel, g]))
  if (method == "2D KDE (jitter + contour)") {
    y <- .norm01(jitter(rep(0, length(v)), amount = 0.35))
    o <- diamond_outline(v, y, lev_q = 0.10)
    if (is.null(o)) return(NULL)
    # 轮廓已被归一化到 [0,1]; 为了可视化相对宽度, 映射回 0..1 原样使用
    return(data.frame(state = st, gene = g, method = method, y = o$y, w = o$w))
  }
  o <- violin_outline(v)
  if (is.null(o)) return(NULL)
  data.frame(state = st, gene = g, method = method, y = o$y, w = o$w)
}

res <- bind_rows(lapply(c("2D KDE (jitter + contour)", "1D KDE (final method)"),
                        function(m)
  bind_rows(lapply(genes, function(g) bind_rows(lapply(states, function(s)
    one_cell(g, s, m)))))))

poly <- res %>%
  group_by(state, gene, method) %>%
  group_modify(function(df, key) {
    sh <- data.frame(w = df$w, y = df$y)
    violin_poly(sh, centre = 0, half_width = 0.42)
  }) %>% ungroup()

p_cmp <- ggplot(poly, aes(x = x, y = y, fill = state, group = interaction(state, gene, method))) +
  geom_polygon(color = "grey30", linewidth = 0.2) +
  facet_grid(gene ~ method, scales = "free_y") +
  scale_fill_manual(values = STATE_COLORS, guide = "none") +
  scale_x_continuous(breaks = NULL) +
  labs(x = NULL, y = expression(sqrt("expression")),
       title = "Fig 2C shape engine comparison (Oligodendrocyte states)",
       subtitle = paste0("Left: faithful 2D KDE contour (as in the original seaborn code) — ",
                         "zero-inflation stretches every shape into a bar spanning the ",
                         "whole jitter range.\nRight: 1D KDE violin used in the final figure.")) +
  theme_paper(base_size = 9) +
  theme(strip.text.y = element_text(angle = 0, hjust = 1, face = "italic", size = 8),
        strip.placement = "outside", strip.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 7.5, color = "grey25"),
        panel.spacing = unit(0.6, "lines"))

ggsave(file.path(FIG_DIR, "fig2C_alt_2d_kde_vs_1d.png"), p_cmp,
       width = 8, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "fig2C_alt_2d_kde_vs_1d.pdf"), p_cmp, width = 8, height = 5)
# 高分辨率版本 (供报告排版)
dir.create(HI_DIR, showWarnings = FALSE, recursive = TRUE)
grDevices::png(file.path(HI_DIR, "fig2C_alt_2d_kde_vs_1d.png"),
               width = 8, height = 5, units = "in", res = 900, type = "cairo")
print(p_cmp)
grDevices::dev.off()
message("对比图已保存: figures/fig2C_alt_2d_kde_vs_1d.png (+ 900 dpi hi-res)")
