# ==============================================================================
# run_all.R — 一键复现 Figure 5E / 5F (R 版)
#
# 用法 (在 RStudio 中 Source 本文件, 或命令行):
#   Rscript run_all.R
# 输出: Figure5_R/figures/ (含 PDF + 300dpi PNG) 与 figures_hires/ (900dpi)
#       Figure5_R/tables/  (逐点曲线与相关系数 CSV)
# ==============================================================================
script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                       error = function(e) getwd())
if (!file.exists(file.path(script_dir, "00_common.R")))
  script_dir <- "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure5_R"
setwd(script_dir)

t0 <- Sys.time()
for (f in c("10_fig5E_distance_curve.R", "20_fig5F_correlation.R")) {
  cat("\n", strrep("=", 72), "\n>>> ", f, "\n", strrep("=", 72), "\n", sep = "")
  source(file.path(script_dir, f), echo = FALSE)
}
cat(sprintf("\nFigure 5E/5F (R) 复现完成, 耗时 %.1f 分钟\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

graphics.off()
rp <- file.path(script_dir, "Rplots.pdf")
if (file.exists(rp)) unlink(rp)
