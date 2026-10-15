# ==============================================================================
# 05_combine_figure.R — 将 Fig 2A / 2B / 2C 组合成整幅 Figure 2
#
# 版面 (参照原文 Figure 2 的排布思路):
#   - 第一行: 2A (左) | 2B (右);
#   - 第二行: 2C (整行, 四个基因 x 状态矩阵横向并排)。
# 输出: figures/Figure2_reproduced.png / .pdf
#
# 说明:
#   1. 依次 source 三个绘图脚本 (它们分别生成全局变量 pA / pB / pC, 并已单独保存
#      各自的图)。
#   2. pA 本身是 patchwork 对象 (Fig 2A 主版面 + 低丰度 inset + T cell inset),
#      直接嵌进外层 patchwork 会触发 "Need 3 panels, but together nrow and ncol
#      only provide 2"; 因此这里用 patchworkGrob() 把它固化成单个 grob 再手工排版,
#      既能保留 Fig 2A 内部的 inset 布局, 又能放进整幅图。
#   3. 画布长宽比必须接近各子图的自然比例, 否则会被拉伸到标题被裁切 / 右侧溢出。
#      经渲染验证: 16 x 11 in 且 2C 独占下半行时最稳。
# ==============================================================================
if (!exists("PROJECT_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure2_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
}
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid)
})

for (f in c("02_fig2A_density.R", "03_fig2B_state_fraction.R", "04_fig2C_violin.R")) {
  message("\n>>> sourcing ", f)
  source(file.path(PROJECT_DIR, f))
}

if (!exists("pA")) stop("pA 未生成, 请检查 02_fig2A_density.R")
if (!exists("pB")) stop("pB 未生成, 请检查 03_fig2B_state_fraction.R")
if (!exists("pC")) stop("pC 未生成, 请检查 04_fig2C_violin.R")

# ---- 把三块子图固化成 grob ----------------------------------------------------
gA <- patchworkGrob(pA)   # Fig 2A: 主版面 + 两个 inset (patchwork 对象)
gB <- patchworkGrob(pB)   # Fig 2B: 四个分面 (patchwork 对象)
gC <- patchworkGrob(pC)   # Fig 2C: 四个大类横向并排 (patchwork 对象)

W_IN <- 16
H_IN <- 11

draw_final <- function() {
  grid.newpage()
  # 4 行: 标题 / 上半 2A|2B / 下半 2C / 脚注留白
  lyt <- grid.layout(4, 4, heights = unit(c(0.075, 0.45, 0.43, 0.045), "npc"))
  pushViewport(viewport(layout = lyt))

  # --- 标题行 (跨越 4 列) ---
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  grid.text("Figure 2 (reproduced in R): Cell-type and cell-state changes across ages",
            gp = gpar(fontsize = 14, fontface = "bold"), y = 0.70)
  grid.text("Allen et al., Cell 186(1):194-208, 2023  |  data: CELLxGENE (MERFISH + snRNA-seq)",
            gp = gpar(fontsize = 8.5, col = "grey30"), y = 0.16)
  popViewport()

  # --- 上半: 2A (左两列) | 2B (右两列) ---
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  grid.draw(gA); popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3:4))
  grid.draw(gB); popViewport()

  # --- 下半: 2C (整行) ---
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:4))
  grid.draw(gC); popViewport()

  popViewport()
}

# PDF: 用 cairo_pdf —— Fig 2A 使用 Arial 无衬线字体, 而默认的 pdf() 设备
# 不支持这种字体 (会报 "invalid font type"); cairo_pdf 能正确嵌入。
grDevices::cairo_pdf(file.path(FIG_DIR, "Figure2_reproduced.pdf"),
                     width = W_IN, height = H_IN)
draw_final()
invisible(dev.off())

png(file.path(FIG_DIR, "Figure2_reproduced.png"), width = W_IN, height = H_IN,
    units = "in", res = 300)
draw_final()
invisible(dev.off())

# 高分辨率版本, 供 DOCX / PPTX 排版 (300 dpi 放大查看会发虚)
dir.create(HI_DIR, showWarnings = FALSE, recursive = TRUE)
grDevices::png(file.path(HI_DIR, "Figure2_reproduced.png"),
               width = W_IN, height = H_IN, units = "in", res = 900,
               type = "cairo")
draw_final()
invisible(dev.off())

message("\n整幅 Figure 2 组合完成 ✔  -> figures/Figure2_reproduced.png / .pdf",
        "\n                          + figures_hires/Figure2_reproduced.png (900 dpi)")
