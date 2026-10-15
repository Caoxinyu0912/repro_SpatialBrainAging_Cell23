# ==============================================================================
# 01_prepare_data.R — 数据准备与缓存
#   1) MERFISH h5ad: 细胞元信息 + panel 内目标基因原始计数 -> cache/merfish_cache.rds
#   2) snRNA-seq h5ad (若存在): 11 个状态细胞的 panel 外基因计数
#      -> cache/snrna_states_cache.rds (供 Fig 2C 使用)
# ==============================================================================

# 定位 00_common.R (兼容 RStudio 交互运行 / Rscript / run_all.R 三种方式)
if (!exists("PROJECT_DIR")) {
  cand <- c(getwd(),
            "C:/Users/17721/Desktop/R-for-bioinformatics/Homework_AI_powered_R/repro/Figure2_R")
  for (d in cand) if (file.exists(file.path(d, "00_common.R"))) {
    source(file.path(d, "00_common.R")); break
  }
  if (!exists("PROJECT_DIR")) stop("找不到 00_common.R, 请先设置工作目录到 Figure2_R/")
}

obj <- load_data()

message("\n========== MERFISH 数据概览 ==========")
message("细胞数: ", nrow(obj$meta))
message("panel 基因数: ", length(obj$genes),
        " | 缓存保留: ", length(obj$genes_kept))
message("\n各年龄细胞数:")
print(table(obj$meta$age))
message("\n各供体 × 年龄:")
print(table(obj$meta$donor_id, obj$meta$age))
message("\n主要细胞类型 (cell_type):")
print(table(obj$meta$cell_type))

message("\n========== snRNA-seq 数据准备 (Fig 2C 补充基因) ==========")
if (file.exists(SNRNA_FILE) || file.exists(SNRNA_CACHE_FILE)) {
  all_target <- unlist(GENES_TO_SHOW, use.names = FALSE)
  genes_from_snrna <- setdiff(all_target, obj$genes)
  sn <- load_snrna_states(genes_needed = genes_from_snrna, states = ALL_STATES)
  message("snRNA 目标细胞数: ", nrow(sn$meta))
  print(table(sn$meta$clust, sn$meta$age))
  message("snRNA 保留基因: ", paste(sn$genes, collapse = ", "))
} else {
  message("未找到 snRNA-seq 文件, 跳过 (Fig 2C 将只使用 MERFISH panel 内基因)。")
  message("如需完整 32 基因, 请从 CELLxGENE 下载 752ed29d-bbb0-45bd-9c57-19edf68ef73e.h5ad")
}

message("\n数据准备完成 ✔  缓存目录: ", CACHE_DIR)
