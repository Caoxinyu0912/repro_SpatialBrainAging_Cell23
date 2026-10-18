# ==============================================================================
# 10_fig5E_distance_curve.R — Figure 5E
#
# 2 行 x 3 列: 行 = Astro(上) / Micro(下); 列 = Juvenile / Young Adult / Old
# 每个子图 3 条曲线: 到 Endo / VLMC / Olig 的距离 vs 激活分数
#
# 计算部分完全沿用作者 compute_binned_values 的定义 (CELL 46):
#   窗口 i 覆盖 (i, i+bin_size], 步长 1 (滑动窗);
#   binned_mean[i] = mean(scores[idx]);
#   binned_std[i]  = np.std(scores[idx], ddof=0) / sqrt(n);
#   binned_mean -= mean;  binned_std -= mean(binned_mean 中心化后);
#   x = arange(len) + bin_size
#
# 但**坐标轴**按发表的原图 Figure 5E 设置, 而不是按仓库里的 plt.ylim:
#   原图 x 轴刻度只有 20 与 80, 曲线从 ~10 铺到 ~85;
#   y 轴 上排 -0.2~0.3 (间隔 0.1), 下排 -0.05~0.15 (间隔 0.05);
#   曲线在面板边界处被裁掉 (原图也是这么裁的)。
#
# 为此本脚本输出两张图:
#   ① fig5E_distance_curve_R  —— **主轴**, 按发表原图的坐标轴
#      (bin_size=10 -> x=10..89; xlim 6~93; 刻度 20/80)
#   ② fig5E_authorstyle_R     —— 按仓库代码原样
#      (bin_size=50 -> x=50..99; xlim 47.5~101.5; 刻度 60/80/100;
#       Astro ylim -0.2~0.5, Micro ylim -0.05~0.12)
#      这张与 Figure5_authorcode/figures/fig5E_authorcode.png 逐点一致 (max|Δ|=7e-16)
# ==============================================================================
if (!exists("FIG5_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure5_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
}

# 主轴: 原图口径; 对照: 仓库代码口径
BIN_PAPER   <- 20L     # x = bin_size..(max_d-1)
MAX_D_PAPER <- 90L
if (nzchar(Sys.getenv("FIG5_BIN_PAPER")))
  BIN_PAPER <- as.integer(Sys.getenv("FIG5_BIN_PAPER"))
BIN_AUTHOR  <- 50L     # 作者 CELL 464 的值
MAX_D_AUTHOR <- 100L

if (nzchar(Sys.getenv("FIG5_BIN_SIZE")))
  BIN_AUTHOR <- as.integer(Sys.getenv("FIG5_BIN_SIZE"))

# ---- 1. 数据 + 激活分数 (只算一次) -------------------------------------------
allg <- unique(unlist(GENE_SETS))
obj <- read_merfish()                    # 全部 374 基因 (对照基因池需要)
meta <- obj$meta
keep <- !is.na(meta$ct)
meta <- meta[keep, ]; xy <- obj$xy[keep, ]
X <- obj$counts[keep, , drop = FALSE]
message("cells: ", nrow(meta), " | panel genes: ", ncol(X))

message("计算激活分数 (scanpy score_genes 原算法) ...")
act <- compute_activation_scores(meta, X, colnames(X))
for (nm in names(act)) meta[[nm]] <- act[[nm]]

chk <- do.call(rbind, lapply(c("Micro", "Astro", "Olig"), function(ctc) {
  key <- c(Micro = "activate_micro", Astro = "activate_astro",
           Olig = "activate_olig")[[ctc]]
  data.frame(cell_type = ctc, t(sapply(AGE_LEVELS, function(a)
    mean(meta[[key]][meta$ct == ctc & meta$age == a]))),
    check.names = FALSE)
}))
names(chk)[-1] <- AGE_LEVELS
message("\n各年龄组平均激活分数 (对照原文 Fig 5A 的年龄趋势):")
print(chk, digits = 4, row.names = FALSE)
write.csv(chk, file.path(TAB_DIR, "fig5A_activation_check_R.csv"), row.names = FALSE)

# ---- 2. 计算某 bin_size 下的全部曲线 -----------------------------------------
compute_curves <- function(bin_size, max_d) {
  xs <- x_axis(bin_size, max_d)
  dat <- list(); tab <- list()
  for (ctc in c("Astro", "Micro")) {
    key <- if (ctc == "Astro") "activate_astro" else "activate_micro"
    m_ct <- meta$ct == ctc
    for (a in AGE_LEVELS) {
      m <- m_ct & meta$age == a
      for (tg in TARGETS) {
        d <- nearest_dist(xy[m, , drop = FALSE],
                          xy[meta$ct == tg$ct, , drop = FALSE])
        bc <- binned_curve(d, meta[[key]][m], bin_size = bin_size, max_d = max_d)
        dat[[length(dat) + 1L]] <- data.frame(
          cell_type = ctc, age = a, target = tg$label,
          distance_um = xs, mean = bc$mean, sem = bc$sem,
          stringsAsFactors = FALSE)
        tab[[length(tab) + 1L]] <- data.frame(
          cell_type = ctc, age = a, target = tg$label,
          n_cells = sum(m), n_target = sum(meta$ct == tg$ct),
          median_dist = median(d))
      }
    }
  }
  list(df = do.call(rbind, dat), summary = do.call(rbind, tab))
}

# ---- 3. 画图函数: 两种口径共用 ----------------------------------------------
COL <- c(Endo = "#F0E68C", VLMC = "#FFD700", Oligo = "#708090")

plot_style <- function(df, ylim, ybrk, xlim, xbrk, age_lab, xlab_title) {
  df$age <- factor(df$age, levels = AGE_LEVELS)
  df$cell_type <- factor(df$cell_type, levels = c("Astro", "Micro"))
  df$target <- factor(df$target, levels = c("Endo", "VLMC", "Oligo"))
  mk <- function(ctc, bottom = FALSE) {
    d <- subset(df, cell_type == ctc)
    ggplot(d, aes(x = distance_um, y = mean, color = target, fill = target)) +
      geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem),
                  alpha = 0.10, colour = NA) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~ age, nrow = 1, labeller = labeller(age = age_lab)) +
      scale_color_manual(values = COL, name = NULL, drop = TRUE) +
      scale_fill_manual(values = COL, name = NULL, drop = TRUE) +
      scale_x_continuous(breaks = if (bottom) xbrk else NULL) +
      scale_y_continuous(breaks = ybrk[[ctc]]) +
      # 用 coord_cartesian 而不是 scale_*(limits=):
      # scale 的 limits 会把超范围值置 NA 导致曲线断开, coord 只是裁切,
      # 与 matplotlib 的 plt.ylim / 默认 xlim 行为一致。
      coord_cartesian(xlim = xlim, ylim = ylim[[ctc]], expand = FALSE) +
      labs(x = if (bottom) xlab_title else NULL,
           y = paste0(ctc, " activation score")) +
      theme_fig5(base_size = 9) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = 9),
            legend.position = "none",
            panel.spacing = unit(0.9, "lines"))
  }
  p_top <- mk("Astro", bottom = FALSE)
  p_bot <- mk("Micro", bottom = TRUE) +
    theme(plot.margin = margin(2, 5.5, 2, 5.5))
  list(top = p_top, bot = p_bot)
}

save_style <- function(st, title, subtitle, name) {
  library(cowplot)
  body <- plot_grid(st$top, st$bot, ncol = 1, align = "v", axis = "lr",
                    rel_heights = c(1, 1))
  ttl <- ggdraw() +
    draw_label(title, x = 0.012, y = 0.78, hjust = 0, vjust = 0.5,
               fontface = "bold", size = 11, colour = "black") +
    draw_label(subtitle, x = 0.012, y = 0.24, hjust = 0, vjust = 0.5,
               size = 7.5, colour = "grey25")
  p <- plot_grid(ttl, body, ncol = 1, rel_heights = c(0.095, 1))
  save_fig(p, name, 9.5, 6.2)
}

# ---- 4. ① 主轴: 按发表原图的坐标轴 -------------------------------------------
res_paper <- compute_curves(BIN_PAPER, MAX_D_PAPER)
write.csv(res_paper$df, file.path(TAB_DIR, "fig5E_distance_curve_R.csv"),
          row.names = FALSE)
write.csv(res_paper$summary,
          file.path(TAB_DIR, "fig5E_distance_summary_R.csv"), row.names = FALSE)

st1 <- plot_style(
  res_paper$df,
  ylim  = list(Astro = c(-0.2, 0.3), Micro = c(-0.05, 0.15)),
  ybrk  = list(Astro = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3),
               Micro = c(-0.05, 0, 0.05, 0.10, 0.15)),
  xlim  = c(10, 92),                # 使刻度 20/80 落在与原图相同的位置
  xbrk  = c(20, 80),
  age_lab = AGE_LABELS,             # Juvenile / Young Adult / Old
  xlab_title = expression("Distance to specified cells ("*mu*"m)"))
save_style(st1,
           "Figure 5E (reproduced in R): Activation score vs distance to neighbouring cells",
           paste0("Mean (line) \u00b1 SEM (shade); each curve is centred by subtracting ",
                  "its mean across all distances.  bin_size = ", BIN_PAPER,
                  " \u00b5m; axes follow the published figure."),
           "fig5E_distance_curve_R")

# ---- 5. ② 对照: 按仓库代码原样 (bin_size = 50, x = 50..99) -------------------
res_auth <- compute_curves(BIN_AUTHOR, MAX_D_AUTHOR)
write.csv(res_auth$df, file.path(TAB_DIR, "fig5E_authorstyle_R.csv"),
          row.names = FALSE)

st2 <- plot_style(
  res_auth$df,
  ylim  = list(Astro = c(-0.2, 0.5), Micro = c(-0.05, 0.12)),
  ybrk  = list(Astro = seq(-0.2, 0.5, by = 0.1),
               Micro = seq(-0.04, 0.12, by = 0.02)),
  xlim  = c(47.5, 101.5),
  xbrk  = c(60, 80, 100),
  age_lab = c("4wk" = "4wk", "24wk" = "24wk", "90wk" = "90wk"),
  xlab_title = expression("Distance to specified cells ("*mu*"m)"))
save_style(st2,
           "Figure 5E (author-code style): bin_size = 50, x = 50..99",
           paste0("Same computation as the author's CELL 464/467-470 ",
                  "(bin_size = 50, max_d = 100). Numerically identical to the ",
                  "Python author-code run (max |\u0394| = 7e-16)."),
           "fig5E_authorstyle_R")

message("Fig 5E (R) 完成 ✔  主轴 bin=", BIN_PAPER, " | 对照 bin=", BIN_AUTHOR)
