# ==============================================================================
# run_all.R — 一键运行完整复现流程 (数据准备 -> Fig 2A -> 2B -> 2C -> 拼版)
# 用法: 在 RStudio 中打开本文件, 点击 Source; 或命令行:
#   Rscript run_all.R
# ==============================================================================
script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
if (!file.exists(file.path(script_dir, "00_common.R"))) script_dir <- getwd()
setwd(script_dir)
cat("Working directory:", getwd(), "\n")

t0 <- Sys.time()
for (f in c("01_prepare_data.R", "02_fig2A_density.R",
            "03_fig2B_state_fraction.R", "04_fig2C_violin.R",
            "05_combine_figure.R", "06_alt_2d_kde_compare.R")) {
  cat("\n", strrep("=", 70), "\n", f, "\n", strrep("=", 70), "\n", sep = "")
  source(f, echo = TRUE, max.deparse.length = 200)
}
cat("\n全流程完成! 总耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "分钟\n")
cat("图片输出目录:", file.path(script_dir, "figures"), "\n")
cat("统计表输出目录:", file.path(script_dir, "tables"), "\n")

# 命令行下 Rscript 若在无图形设备时被触发绘图, 会生成一个多余的 Rplots.pdf;
# 这里统一关闭设备并清掉该文件, 保持输出目录整洁。
graphics.off()
rp <- file.path(script_dir, "Rplots.pdf")
if (file.exists(rp)) unlink(rp)