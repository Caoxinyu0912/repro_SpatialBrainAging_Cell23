# ==============================================================================
# 04_fig2C_violin.R — Figure 2C: 各细胞状态示例基因表达小提琴矩阵
#
# 版面 (与原文 Figure 2C 一致):
#   - 四个大类 (Olig / Micro / Endo / Astro) 横向并排, 每个大类一个
#     "基因 (行) x 状态 (列)" 矩阵: 基因为行 (行标签为斜体基因名),
#     细胞状态为列 (列标题只画在最上面一个基因行上方);
#   - 每个 gene x state 一格一个对称小提琴, 填充色 = 该状态配色
#     (与 Fig 2B 的深浅渐变完全一致);
#   - 每类 8 个基因, 共 32 个 (与原文 genes_to_show 一字不差)。
#   - 无 y 轴刻度数字 (每行独立缩放, 数字无意义; 原文亦不显示)。
#     精确数值见 tables/fig2C_expression_summary.csv。
#
# 数据来源:
#   - MERFISH 374 基因 panel 内的 15 个基因: MERFISH 细胞原始计数, sqrt 变换
#     (对应原文 layer='sqrt_norm');
#   - panel 外的 15 个基因: 用**同状态** (联合聚类注释 clust_annot) 的
#     snRNA-seq 细胞实测表达近似 (sqrt(CP10K)), 基因名后带 * 标注。
#
#   关于未使用 BrainAgingSpatialAtlas_Imputed.h5ad 的说明:
#   原文 Fig 2C 用的是 snRNA-seq 逐细胞插补矩阵 (adata_imputed 的
#   layers['sqrt_norm'])。CELLxGENE 发布的该文件 X 已被**按基因中心化**
#   (约 73% 元素为负), 且 layers 为空、var/means 加回后仍有负值, 无法还原
#   原始 sqrt(CP10K); 强行截断负值会在 0 处引入虚假点质量, 使小提琴形状失真。
#   因此改用同状态 snRNA 实测表达 —— 透明、可复现, 且与插补矩阵共享同一套
#   状态注释 (clust_annot), 比较仍然成立。
#
# 形状算法: violin_shape() / violin_poly() (见 00_common.R) —— 对表达量做
#   一维高斯核密度估计, 再上下对称填充。原文用二维核密度填充轮廓, 但 MERFISH
#   表达高度零膨胀, 「表达量 = 0」的点质量会把二维等高线拉成长条, 详见
#   00_common.R 中 violin_outline() 的说明与 README 的差异分析。
#
# 输出: figures/fig2C_violin_states.pdf/png, tables/fig2C_expression_summary.csv
# ==============================================================================
if (!exists("PROJECT_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure2_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
}
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

VIOLIN_BW_MULT <- 1      # 核带宽倍数
HALF_WIDTH     <- 0.44   # 小提琴半宽 (相对状态间距)
set.seed(20230105)       # 固定随机种子, 保证可重复

obj <- load_data()
meta <- obj$meta
counts <- obj$counts

all_target_genes <- unique(unlist(GENES_TO_SHOW, use.names = FALSE))
# 注意: 判断"能否用 MERFISH 实测值"要看 **counts 的列** (只保留了 15 个基因),
# 不能用 obj$genes —— 那是全部 374 个 panel 基因的符号表。
# 用 obj$genes 会误判 30 个基因都在 panel 内, 于是 panel 外的 15 个基因
# 既拿不到 MERFISH 值、也不会走 snRNA 分支, 最终从图里整行消失。
genes_in_panel   <- intersect(all_target_genes, colnames(counts))
genes_from_snrna <- setdiff(all_target_genes, genes_in_panel)
message("Fig 2C: 共 ", length(all_target_genes), " 个基因 | ",
        "MERFISH 实测 ", length(genes_in_panel), " 个: ",
        paste(genes_in_panel, collapse = ", "))
message("        由 snRNA-seq 补充 ", length(genes_from_snrna), " 个: ",
        paste(genes_from_snrna, collapse = ", "))

snrna <- load_snrna_states(genes_needed = genes_from_snrna, states = ALL_STATES)

#' dgCMatrix (行 = 细胞) 每行除以该细胞 library size 再乘 1e4 -> CP10K
cp10k <- function(m, lib) {
  m@x <- m@x / lib[m@i + 1L] * 1e4
  m
}

# ---- 1. 组装每个 gene x state 的表达向量 --------------------------------------
PARENT_ORDER_C <- c("Olig", "Micro", "Endo", "Astro")
df_all <- list()

for (p in PARENT_ORDER_C) {
  genes   <- GENES_TO_SHOW[[p]]
  genes_m <- intersect(genes, colnames(counts))
  genes_s <- intersect(genes, colnames(snrna$counts))

  # MERFISH: 该大类细胞 x 目标基因 (sqrt 原始计数)
  rows <- which(meta$clust_annot %in% STATES_OF[[p]])
  em <- as.matrix(counts[rows, genes_m, drop = FALSE])
  dm <- data.frame(state = meta$clust_annot[rows], em,
                   check.names = FALSE, stringsAsFactors = FALSE)

  # snRNA-seq: 同状态细胞 x panel 外基因 (sqrt(CP10K))
  srows <- which(snrna$meta$clust %in% STATES_OF[[p]])
  if (length(srows) > 0 && length(genes_s) > 0) {
    sm <- as.matrix(cp10k(snrna$counts[srows, genes_s, drop = FALSE],
                          snrna$meta$total_counts[srows]))
    ds <- data.frame(state = snrna$meta$clust[srows], sm,
                     check.names = FALSE, stringsAsFactors = FALSE)
  } else ds <- NULL

  df_long <- bind_rows(
    if (length(genes_m))
      pivot_longer(dm, all_of(genes_m), names_to = "gene", values_to = "value") %>%
        mutate(value = sqrt(value)) else NULL,
    if (!is.null(ds) && length(genes_s))
      pivot_longer(ds, all_of(genes_s), names_to = "gene", values_to = "value") %>%
        mutate(value = sqrt(value)) else NULL) %>%
    mutate(parent = p,
           state = factor(state, levels = STATES_OF[[p]]),
           gene  = factor(gene,  levels = genes),
           label = factor(ifelse(gene %in% genes_m, as.character(gene),
                                 paste0(as.character(gene), "*")),
                          levels = ifelse(genes %in% genes_m, genes,
                                          paste0(genes, "*"))))
  df_all[[p]] <- df_long
  message(sprintf("Fig 2C [%s] MERFISH: %s | snRNA: %s | 细胞数: %d", p,
                  paste(genes_m, collapse = ","), paste(genes_s, collapse = ","),
                  length(rows)))
}
df_all <- bind_rows(df_all)

# ---- 1b. 强制基因行的显示顺序 ----------------------------------------------
# 背景: label 原本是一列**全局** factor, 四个大类共用一套层级。之前这些层级
#   是在 group_by/pivot 之后用 ifelse 现场构造的, 顺序与 GENES_TO_SHOW 不一致;
#   facet_grid 按"层级顺序"排面板, 于是实测 Endo 面板第一行是 B2m、
#   Astro 面板第一行是 C4b —— 与你看到的现象一致。
#
# 关键点: factor 一旦创建, 其层级顺序就固定下来, 事后对子集调用 droplevels()
#   只会**删除未使用层级**, 不会重排。所以必须用正确的顺序**重建整列**。
#
# 修法: 按 PARENT_ORDER_C 逐类把每个基因映射到 GENES_TO_SHOW 的下标, 拼出
#   各面板的 32 个 label; 再把这 32 个标签的**去重顺序**作为整列的层级
#   (Olig 8 → Micro 8 → Endo 8 → Astro 8)。这样渲染每个大类时,
#   facet_grid 的行顺序就是该大类的 genes_to_show 顺序。
label_chr <- rep(NA_character_, nrow(df_all))
for (p in PARENT_ORDER_C) {
  genes_p  <- GENES_TO_SHOW[[p]]
  # 判断"是否在 MERFISH panel 内"必须用 counts 的列名 (obj$genes 不是面板列表)
  genes_m  <- intersect(genes_p, colnames(counts))
  levels_p <- ifelse(genes_p %in% genes_m, genes_p, paste0(genes_p, "*"))

  idx <- which(df_all$parent == p)
  stopifnot(length(idx) > 0)
  g_chr <- as.character(df_all$gene[idx])
  pos <- match(g_chr, genes_p)                # 1..8, 未命中为 NA
  if (anyNA(pos)) stop("[", p, "] 出现表外基因: ",
                       paste(unique(g_chr[is.na(pos)]), collapse = ","))

  label_chr[idx] <- levels_p[pos]
}
stopifnot(!anyNA(label_chr))
df_all$label_chr <- label_chr                 # 先存字符串

# 关键: label 必须**按 parent 分别设定层级**。
# 若做成整列共用一个 30 层级 factor, B2m (Micro 第 8 行, 同时是 Endo 第 3 行)
# 与 C4b (Olig 第 5 行, 同时是 Astro 第 8 行) 会停在"第一次出现"的位置,
# 于是 Endo 第一行变成 B2m、Astro 第一行变成 C4b。
# 因此下面在绘图循环里 (每个大类单独构建数据) 用 factor(label_chr, levels = 该类 8 个)
# 重建层级。

# 自检: 每个大类的行顺序必须与原文 genes_to_show 完全一致
for (p in PARENT_ORDER_C) {
  genes_p  <- GENES_TO_SHOW[[p]]
  genes_m  <- intersect(genes_p, colnames(counts))
  levels_p <- ifelse(genes_p %in% genes_m, genes_p, paste0(genes_p, "*"))
  idx <- which(df_all$parent == p)
  lv <- levels(factor(df_all$label_chr[idx], levels = levels_p))
  if (!identical(lv, levels_p)) {
    message("[DIAG] ", p)
    message("  期望: ", paste(levels_p, collapse = " | "))
    message("  实际: ", paste(lv, collapse = " | "))
    stop("[", p, "] 行顺序与 genes_to_show 不一致")
  }
}
message("基因行顺序自检通过: 每个大类的行顺序 = genes_to_show")

# ---- 2. 每个 gene x state 计算 1D KDE 小提琴多边形 ----------------------------
message("\nComputing KDE violin polygons ...")
# 同一大类同一基因的所有状态共用 y 轴范围, 便于同行互相对比
env_tab <- df_all %>%
  group_by(parent, label_chr) %>%
  summarise(from = min(value, na.rm = TRUE), to = max(value, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(pad = (to - from) * 0.10 + 1e-9,
         from = pmax(0, from - pad),     # 表达量非负, 下界不小于 0
         to   = to + pad)
violin_env <- setNames(lapply(seq_len(nrow(env_tab)), function(i)
  c(env_tab$from[i], env_tab$to[i])),
  paste(env_tab$parent, env_tab$label_chr))

#' 单组表达向量 -> 小提琴轮廓 data.frame(w, y); 全常量组给极小抖动
shape_or_flat <- function(v, from, to) {
  if (length(unique(v)) < 2L) v <- v + abs(stats::rnorm(length(v), 0, 2e-3))
  sh <- violin_shape(v, bw_mult = VIOLIN_BW_MULT, from = from, to = to)
  if (is.null(sh)) return(data.frame(w = numeric(0), y = numeric(0)))
  data.frame(w = sh$w, y = sh$y)
}

shapes <- df_all %>%
  group_by(parent, label_chr, state, gene) %>%
  group_modify(function(df, key) {
    env <- violin_env[[paste(as.character(key$parent), as.character(key$label_chr))]]
    shape_or_flat(df$value, env[1], env[2])
  }) %>% ungroup()

# 面板内相对列位置 (1..k): 每个大类单独编号
state_pos <- do.call(rbind, lapply(PARENT_ORDER_C, function(p)
  data.frame(parent = p, state = STATES_OF[[p]],
             state_idx = seq_along(STATES_OF[[p]]), stringsAsFactors = FALSE)))

poly <- shapes %>%
  ungroup() %>%
  left_join(state_pos, by = c("parent", "state")) %>%
  group_by(parent, label_chr, state, gene, state_idx) %>%
  group_modify(function(df, key)
    violin_poly(data.frame(w = df$w, y = df$y), centre = key$state_idx,
                half_width = HALF_WIDTH)) %>%
  ungroup() %>% filter(is.finite(x), is.finite(y))
message("  polygon vertices: ", nrow(poly))

# ---- 3. 绘制: 基因(行, 斜体标签) x 状态(列, 标题) -----------------------------
plot_list <- list()
for (p in PARENT_ORDER_C) {
  genes_p  <- GENES_TO_SHOW[[p]]
  genes_m  <- intersect(genes_p, colnames(counts))
  levels_p <- ifelse(genes_p %in% genes_m, genes_p, paste0(genes_p, "*"))
  # 每个大类单独构建数据, 并**按该类的 8 个层级**重建 label factor;
  # 这样 facet_grid 的行顺序 = genes_to_show, 且不会受其它大类
  # 共享基因 (B2m / C4b) 的层级位置影响。
  pl <- poly %>% filter(parent == p) %>%
    mutate(label = factor(label_chr, levels = levels_p))
  stopifnot(identical(levels(droplevels(pl$label)), levels_p))
  states_p <- STATES_OF[[p]]
  top_gene <- levels_p[1]
  header <- data.frame(state = factor(states_p, levels = states_p),
                       x = seq_along(states_p),
                       label = factor(top_gene, levels = levels_p),
                       y = max(pl$y[pl$label == top_gene], na.rm = TRUE))
  pp <- ggplot() +
    geom_polygon(data = pl,
                 aes(x = x, y = y, fill = state,
                     group = interaction(label, state)),
                 color = "grey35", linewidth = 0.18) +
    geom_text(data = header, aes(x = x, y = y, label = state),
              inherit.aes = FALSE, vjust = -0.45, size = 2.6, color = "black") +
    facet_grid(label ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = STATE_COLORS, guide = "none") +
    scale_x_continuous(limits = c(0.5, length(states_p) + 0.5),
                       breaks = NULL, expand = c(0.04, 0)) +
    # 顶部留 26% 给"状态列标题", 底部留 4%
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.26))) +
    labs(x = NULL, y = expression(sqrt("expression")), title = p) +
    theme_paper(base_size = 9) +
    # 基因名改为**横向**文字 (angle = 0), 放在每行左侧; hjust = 1 让文字右对齐、
    # 紧贴分面左缘。横向后每个标签要占 ~48pt 宽度, 因此把 plot.margin 左侧加大。
    theme(strip.text.y = element_text(angle = 0, hjust = 1, face = "italic",
                                      size = 6.5, colour = "black"),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1, face = "italic",
                                           size = 6.5, colour = "black"),
          panel.spacing.y = unit(0.05, "lines"),
          # 每行基因量程不同 (free_y), 逐行画刻度会让相邻行互相压叠; 原文亦不显示。
          # 左侧留出足够空白: 竖排的基因名 + y 轴标题 "√expression" 都需要空间,
          # 否则基因名会被 y 轴标题压住。
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # 去掉坐标轴线 (原文该 panel 只有图形, 没有轴线)
          axis.line = element_blank(),
          axis.title.y = element_text(size = 6.5),
          plot.title = element_text(face = "bold.italic", size = 10, hjust = 0.5),
          # 左侧留白要同时容纳横向基因名 (~48pt) 与 y 轴标题
          plot.margin = margin(16, 4, 2, 46))
  plot_list[[p]] <- pp
}

pC <- wrap_plots(plot_list, nrow = 1) +
  plot_annotation(
    title = "Fig 2C (reproduced): Expression of example genes across cell states",
    subtitle = paste0("One violin per gene x state; fill = state colour (as in Fig 2B). ",
                      "Genes marked * are absent from the MERFISH panel and are shown as ",
                      "sqrt(CP10K) of snRNA-seq nuclei with the same joint-cluster state."),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 7, color = "grey25")))
save_fig(pC, "fig2C_violin_states", 13, 6.0)

# ---- 4. 表达汇总表 (报告可直接引用) -------------------------------------------
expr_tab <- df_all %>%
  group_by(parent, gene, label_chr, state) %>%
  summarise(n_cells = n(),
            pct_pos = mean(value > 0),
            mean_sqrt = mean(value),
            median_sqrt = median(value),
            .groups = "drop") %>%
  rename(label = label_chr) %>%
  arrange(parent, gene, state)
write.csv(expr_tab, file.path(TAB_DIR, "fig2C_expression_summary.csv"), row.names = FALSE)
message("  saved table: fig2C_expression_summary.csv (", nrow(expr_tab), " rows)")

# 供 05_combine_figure.R 使用
pC_main <- pC
message("\nFig 2C 完成 ✔")
