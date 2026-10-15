# Figure 2 复现项目（R 语言）

复现论文 **Allen W.E., Blosser T.R., Sullivan Z.A., Dulac C., Zhuang X. (2023).
"Molecular and spatial signatures of mouse brain aging at single-cell resolution."
Cell 186(1):194–208.e18. DOI: [10.1016/j.cell.2022.12.010](https://doi.org/10.1016/j.cell.2022.12.010)**
中的 **Figure 2**（细胞类型与细胞状态组成随年龄的变化）。

- 原文代码（Python）：`SpatialBrainAgingCell22-main/`（Zhuang Lab）
- 数据（均为 CELLxGENE 同一 collection 31937775-0602-4e52-a799-b6acdd2bac2e）：

  | 文件 | 内容 | 用途 |
  |---|---|---|
  | `BrainAgingSpatialAtlas_MERFISH.h5ad` | 378,918 细胞 × 374 基因 panel | Fig 2A/2B/2C 主数据 |
  | `BrainAgingSpatialAtlas_snRNAseq.h5ad` | 79,667 核 × 20,929 基因（全转录组） | Fig 2C 补充 panel 外的 15 个基因 |
  | `BrainAgingSpatialAtlas_MERFISH_LPS.h5ad` | LPS 处理数据（同 374 panel） | 本项目未使用（对应原文 Fig 5–6） |
  | `BrainAgingSpatialAtlas_Imputed.h5ad` | 17 GB snRNA 插补数据 | 本项目未使用（内存不足，见「与原文的差异」） |

  数据来源与可访问性详见上一级目录的 `数据来源与可访问性说明.md`。

## 环境要求

- R ≥ 4.2（本机已装 R 4.6.1 + RStudio）
- R 包：`hdf5r, Matrix, dplyr, tidyr, ggplot2, patchwork, RColorBrewer, scales`（本机均已安装）
- 无需 Python / reticulate（直接用 hdf5r 解析 h5ad）

## 运行方法

在 RStudio 中打开 `run_all.R`，点击 **Source**；或在终端执行：

```
Rscript run_all.R
```

首次运行会解析 h5ad 并生成缓存（MERFISH 约 1 分钟；snRNA-seq 2.5 GB 分块读取约 2–5 分钟）；
之后重复运行直接读缓存，**约 1.3 分钟**完成全部 6 个脚本。

## 各面板说明（与原文格式一致）

### Fig 2A — 细胞类型密度 (cells/mm²) 随年龄变化
- 主版面：全部 13 个细胞类型；**一级 inset** 嵌在主版面右上空白区
  （放大低丰度类型 VLMC / Epen / Macro / T cell）；**二级 inset** 嵌在一级 inset
  右上空白区（仅放大 T cell）——与原文的嵌套放大版式一致（因版面问题并未实现）；
- 显著差异 (Juvenile vs Old, 独立样本bootstrap t 检验, BH FDR < 0.05) 的类型画**方括号**
  （只跨 4wk 与 90wk 两根柱）并在顶端标 `*`；不显著的类型不画方括号、不标 `ns`；
- **柱宽与间距规则**（`DODGE_W` / `BAR_W` / `ERR_CAP`，脚本内有 `stopifnot` 自检）：
  - 同一细胞类型的 3 根柱为"一组"，组内相邻柱中心距 = `DODGE_W/3`，
    组间空隙 = `1 - DODGE_W`，要求 **组间空隙 > 组内中心距**（即 `DODGE_W < 0.75`），
    这样一眼能分辨哪些柱属于同一细胞类型；当前 `DODGE_W = 0.72`（组内 0.24 / 组间 0.28）；
  - 柱宽 `BAR_W = DODGE_W/3*0.75 = 0.18`，柱间仍留 25% 空隙；
  - 误差棒帽宽 `ERR_CAP = BAR_W*0.36`，明显窄于柱宽，避免被误看成第二根柱子；
- **y 轴上限放到 600**（刻度仍只标到 500）：若沿用 520，ExN 的 95% CI 上帽
  （24wk 546、4wk 527）会被裁掉，表现为"误差棒标注不全"；
- 年龄配色直接采样自原文图例：**#7296D8 (Juvenile) / #D2C0D2 (Young adult) /
  #DE8A8A (Old)**。

### Fig 2B — 细胞状态组成的 100% 堆叠柱状图
- 每个大类 (Olig/Micro/Endo/Astro) 一个分面，分面内 3 根竖柱
  (Juvenile / Young adult / Old)，状态纵向堆叠至 100%；
- 同大类内状态**深浅渐变**配色（逐块采样自原文，见 `tools/sample_paper_colors.py`）：
  Oligo `#211F1F → #919396 → #CCCCD2`；Micro `#96185A → #C91778 → #FC1191`；
  Endo `#8A8454 → #BAB46C → #EAE48A`；Astro `#246642 → #2A8A54`；
- 显著变化的状态在 Juvenile 与 Old 柱的对应段内标 `*`（段高 < 5% 时不标）。

### Fig 2C — 基因 × 状态表达小提琴矩阵（每类 8 个基因，共 32 个）
- 四个大类横向并排；每个大类内**基因为行**（斜体基因名作行标签）、
  **状态为列**（状态名作列标题）；
- 小提琴填充色与 Fig 2B 的状态渐变一致；
- **数据来源**（图中带 `*` 的基因为 snRNA-seq 来源）：
  - MERFISH panel 内的 15 个基因：MERFISH 细胞原始计数，sqrt 变换（对应原文 `layer='sqrt_norm'`）；
  - panel 外的 15 个基因（Ank2, Trim2, Robo1, Spock1, Dgki, Dpp10, Meg3, Ly6a,
    Rbm20, Tmem209, Lsamp, Slc1a2, Luzp2, Trpm3, Brinp3）：用**同状态**
    （联合聚类注释 `clust_annot`）的 snRNA-seq 细胞实测表达近似，sqrt(CP10K) 归一化。
    原文做法是基于 snRNA-seq 对每个 MERFISH 细胞做插补（结果在 17 GB 的
    `BrainAgingSpatialAtlas_Imputed` 数据集中，本机内存无法使用）。

## 脚本说明

| 脚本 | 作用 |
|---|---|
| `00_common.R` | 公共配置：路径、配色、32 基因列表、h5ad 读取函数（分块、内存友好）、核密度小提琴引擎 |
| `01_prepare_data.R` | 读取 MERFISH + snRNA-seq → 缓存到 `cache/` |
| `02_fig2A_density.R` | Fig 2A：密度柱状图 + 嵌套 inset + bracket 显著性 |
| `03_fig2B_state_fraction.R` | Fig 2B：100% 堆叠柱状图 |
| `04_fig2C_violin.R` | Fig 2C：基因 × 状态小提琴矩阵（32 基因） |
| `05_combine_figure.R` | 用 patchwork 将 A/B/C 拼成整幅 Figure 2 |
| `06_alt_2d_kde_compare.R` | 备选方案对照：二维 KDE 轮廓 vs 一维 KDE 小提琴（失败环节留证） |
| `run_all.R` | 一键运行全部流程 |

`tools/` 目录存放报告与幻灯片的生成脚本（Python，非复现流程必需）：
`make_report.py`（复现报告 DOCX）、`make_slides_v3.py`（汇报 PPTX）、
`make_comparison.py`（原文 vs 复现并排对比图）、
`fix_pdf_images.py`（修复导出 PDF 的图片分辨率，见下）、
`make_report_pdf.ps1`（一键 DOCX→PDF 并自动修复图片）、
`paper_fig2a_recompute.py`（按原文代码复算 2A 密度，用于核对）、
`verify_boot_ci.py`（核实 bootstrap 次数 / 逐组执行 / seed 影响）、
`diag_2c_zero.py`（诊断 2C 底部宽横条的来源）、
`audit_report.py`（报告数字与 CSV 逐项核对）、
`robustness_check.py`（**稳健性检验**，见下）。

### 稳健性检验（重要，报告 5.6 节）

`tools/robustness_check.py` 在供体级明细表上同时做四种检验并统一 BH 校正：

```bash
python tools/robustness_check.py     # -> tables/fig2A_robustness.csv
```

| 方法 | BH 显著类型 |
|---|---|
| 独立样本 t 检验（等方差）← 原文与本项目采用 | Olig、OPC、VLMC |
| Welch t（不等方差） | **无**（Olig 0.052、VLMC 0.055、OPC 0.084） |
| Mann-Whitney U | **无**（Olig/OPC/VLMC/Epen 均 0.052） |
| 置换检验（20,000 次） | Olig、OPC、VLMC |

**即“恰好 3 个显著类型”只在原文所用的等方差 t 检验下成立，且处于临界位置。**
直接原因是两组方差极不齐（Olig 的 4wk 组 sd = 7.6 vs 90wk 组 sd = 63.1，相差 8.3 倍），
而每组只有 3~5 只动物。报告已在 5.6 节与差异表 D10 中如实披露。

> 注意：生成报告的 python 环境没有 scipy，因此该检验结果先落盘为
> `tables/fig2A_robustness.csv`，再由 `make_report.py` 读取成表。

### 报告三张数据表均由 CSV 动态生成

`make_report.py` 里 `fig2a_rows()` / `fig2b_rows()` / `fig2c_rows()` /
`robustness_rows()` 分别读取
`fig2A_*.csv`、`fig2B_state_fraction_summary.csv` + `fig2B_significance.csv`、
`fig2C_expression_summary.csv`、`fig2A_robustness.csv`。
报告中的数字不再手抄，口径一变图文不会脱节（此前 2A 曾因手抄旧口径值而不一致）。

### 报告 PDF 的图片分辨率（重要）

LibreOffice 把 DOCX 导出为 PDF 时会**对内嵌位图做降采样**——实测一张
10800×5400 的图被压到 737×634，有效分辨率只剩 **144 dpi**，放大查看明显发虚。

因此报告 PDF 必须走两步：

```powershell
# 在 repro 目录下执行
pwsh -File Figure2_R\tools\make_report_pdf.ps1
```

1. `libreoffice-kit` 把 DOCX 转成 PDF；
2. `tools/fix_pdf_images.py` 逐页找出被降采样的小图，按同样的显示位置
   插入 `figures_hires/` 的高分辨率版本（缩放到 **600 dpi**）并删除小图。

修复后：5 张报告插图全部为 600 dpi，PDF 约 3.5 MB。
（不要直接把 900 dpi 源图原样插入——PDF 会涨到 700 MB 以上。）

另外，`make_report.py` 的正文渲染已支持 `**粗体**` 标记：会解析成真正的
Word 加粗 run，而不是把 `**` 当普通字符输出（此前正文里出现了字面量的 `**`）。

## 输出文件

- `figures/`（300 dpi，日常查看 / 幻灯片）：
  - `fig2A_density`、`fig2B_state_fraction`、`fig2C_violin_states`、`fig2C_violin_by_age`、
    `Figure2_reproduced`（各含 .pdf + .png）
  - `Figure2_side_by_side.png`（原文 Figure 2 与复现并排对比）
  - `paper_Figure2_original.png`（从论文 PDF 裁出的原文 Figure 2，900 dpi 裁剪）
  - `fig2C_alt_2d_kde_vs_1d.png`（二维 KDE 与一维 KDE 的形状对比，失败环节证据）
- `figures_hires/`（**900 dpi**，供 DOCX / PPTX 排版；由 `save_fig()` 自动额外导出）
- `tables/`：各面板统计量（均值、95% bootstrap CI、t 检验 p 值、BH 校正 FDR），可直接放入报告
- `cache/`：数据缓存（删除对应 `.rds` 后重跑可重建）

## 统计方法（与原文一致）

- 密度/比例先在**每只供体**水平计算，再按年龄分组：柱/误差棒 = 跨供体均值 ± 95% CI（t 分布）；
- Juvenile (4wk) vs Old (90wk) 独立样本 t 检验（等方差，同 scipy `ttest_ind` 默认），BH 校正；
- 样本量：MERFISH 4wk n=4 / 24wk n=3 / 90wk n=5 供体；snRNA-seq 仅 4wk 与 90wk（各 2 供体）。

## 复现结果概览（与原文一致的主要发现）

- 抑制性神经元（InN）密度随衰老**显著下降**；Olig 显著**增加**、OPC 显著**下降**；
  Endo、Micro、Peri、Epen 显著变化；ExN / MSN 无显著变化；T 细胞上升但未达显著
  （p = 0.22，与原文「substantial increase … did not reach statistical significance」一致）；
- 状态组成随年龄显著改变：Olig-1 → Olig-2/Olig-3 分化（少突胶质成熟与激活），
  Endo 各状态比例变化显著；激活状态 Olig-3、Micro-3、Astro-2 有升高趋势；
- 状态特异基因表达：Olig-3 高表达 C4b（阳性率 3.5% → 85.9%）/Il33；Micro-3 高表达 Trem2/B2m；
  Astro-2 高表达 Gfap/C4b；Endo-3 高表达 Xdh/Fmo2 —— 与原文「衰老激活态」结论一致。

## 与原文的差异（报告「对比与差异分析」部分可直接引用）

1. **数据版本**：使用 CELLxGENE 处理版导出数据，供体-切片数与原文略有出入
   （本版 `slice` 字段仅 3 个取值，凸包面积只能按供体合并计算）。
2. **密度归一化**：原文 Python 代码使用 `scipy.ConvexHull.area`（二维时实为**周长**）×1000；
   本脚本默认用真正的凸包面积（µm²→mm²），物理意义正确。如需对齐原文口径，
   将 `02_fig2A_density.R` 中 `DENSITY_METHOD <- "area"` 改为 `"paper"`。
3. **Fig 2C 的 15 个 panel 外基因**：原文为 snRNA-seq 插补（逐细胞）；
   本项目用同状态 snRNA 细胞实测分布近似（图中以 `*` 标注），分布形态可能略有差异；
   snRNA-seq 仅有 Juvenile 与 Old 两个年龄。
4. **表达归一化**：原文 MERFISH 用 sqrt(计数/细胞体积)；本导出数据无体积信息，用 sqrt(计数) 近似；
   snRNA 基因用 sqrt(CP10K)。
5. **Fig 2C 的形状算法**：原文用二维核密度填充轮廓（`sns.kdeplot(x=expr, y=jitter, fill=True)`）；
   本项目在 R 中**完整实现了这条路线**（`00_common.R` 的 `kde2d_binned` / `.ms_lines` / `.ms_ring`），
   但 MERFISH 表达高度零膨胀（Chrm3 仅 4.8% 阳性、Bgn 14.0%、Fmo2 18.2%），
   「表达量 = 0」的点质量沿抖动方向近似均匀，使二维密度场沿 y 方向近乎平顶，
   任何能包含足够概率质量的等高线都会被拉成横跨整个抖动区间的细长条。
   因此最终改用**一维核密度小提琴**（`violin_outline()`，并把求值下界限制在 0 以避免虚假平底），
   形状构造原理与原文一致，且不受零膨胀影响。对照证据见
   `figures/fig2C_alt_2d_kde_vs_1d.png`（由 `06_alt_2d_kde_compare.R` 生成）。
6. **细胞类型命名**：CELLxGENE 中血管软脑膜细胞标签为 `Vlmc`，原文代码写作 `VLMC`；
   本脚本在读取阶段统一归一化为 `VLMC`，否则 Fig 2A 中该类型密度会全为 0。
7. snRNA-seq 各状态细胞数较少（如 Endo-2 仅 39 个核），其小提琴分辨率有限。
8. **未使用 17 GB 插补矩阵的原因（已实测验证）**：原文 Fig 2C 用的是 snRNA-seq
   逐细胞插补矩阵（`adata_imputed` 的 `layers['sqrt_norm']`）。本项目检查了
   CELLxGENE 发布的 `BrainAgingSpatialAtlas_Imputed.h5ad`（378,918 细胞 × 7,276 基因），
   发现：
   - 该文件的 `layers` 为空，无法按原代码取 `sqrt_norm` 层；`X` 已被**按基因中心化**
     （抽样检查约 73% 元素为负，`X` 值域如 C4b ∈ [-0.61, 9.20]）；
   - 把 `var/means` 加回去仍存在负值（C4b 最小 −0.46、Il33 最小 −0.64），
     无法还原原始 `sqrt(CP10K)`；
   - 若强行把负值截断为 0，会在 0 处注入一个巨大的**虚假点质量**，使低表达基因的
     小提琴出现"底部横杠 + 细尖刺"的失真形状。

   因此最终改用**同状态 snRNA-seq 实测表达**近似。这个替代方案是透明的、可复现的，
   并且与插补矩阵共享同一套状态注释（`clust_annot`），所以状态间的比较仍然成立；
   代价是带 `*` 的基因分布比原文更离散（见第 3 条）。

## 修改指引

- 换数据路径 / 输出路径：改 `00_common.R` 顶部路径变量
- 换密度计算口径：`02_fig2A_density.R` 的 `DENSITY_METHOD`
- 增减 Fig 2C 基因：`00_common.R` 的 `GENES_TO_SHOW` / `EXTRA_GENES`
  （新增 MERFISH panel 内基因后需删除 `cache/merfish_cache.rds` 重跑）
- 调小提琴平滑度 / 宽度：`04_fig2C_violin.R` 的 `VIOLIN_BW_MULT` / `HALF_WIDTH`
- 重新生成报告 / 幻灯片：`python tools/make_report.py`、`python tools/make_slides_v3.py`

## 引用

- 论文：Allen et al., Cell 2023, DOI 10.1016/j.cell.2022.12.010（PMC10024607）
- 原代码：https://github.com/ZhuangLab/SpatialBrainAgingCell22
- 数据：https://cellxgene.cziscience.com/collections/31937775-0602-4e52-a799-b6acdd2bac2e
  （MERFISH: c93d78c2…；snRNA-seq: 752ed29d…；LPS: 9a610e48…）

---

## 版面修订记录

三个 panel 的最终版式如下（脚本 `02` / `03` / `04` / `05`）。

### Fig 2A（统计口径，已按作者原始代码逐条对齐）

作者代码 `analyze_merfish_celltyping.ipynb` CELL 308/309：

```python
all_ct = ["ExN","InN","MSN","Olig","Astro","Endo","Micro","OPC","Peri",
          "Vlmc","Epen","Macro","T cell"]
for i in np.arange(1, adata_combined.obs.batch.max()):
    curr_batch = adata_combined[adata_combined.obs.batch == i]
    curr_batch_area = 0
    for j in curr_batch.obs.slice.unique():
        hull = ConvexHull(curr_slice.obsm['spatial'])
        curr_batch_area += hull.area        # 2D 下 scipy 返回的是**周长**
    for ct in all_ct:
        counts.append(1000 * curr_batch[curr_batch.obs.cell_type == ct].shape[0]
                      / curr_batch_area)
```

逐条对应关系：

| 原文口径 | 本版 R 实现 |
|---|---|
| `batch` = 一个动物 | `donor_id`（CELLxGENE 里没有 `batch` 字段） |
| 面积 = 该 batch 各 `slice` 凸包之和 | 同左；`spatial_coords` 单位 nm，`1000*n/perimeter` 即把 nm 当 µm |
| 统计单位 = batch（动物） | 每动物一个观测值，n = 4 / 3 / 5；`stopifnot(max(n) < 20)` 兜底 |
| 检验 | 只做 4wk vs 90wk，`t.test(var.equal = TRUE)`，`p.adjust(method="BH")` |
| **误差线** | **bootstrap 95% CI**（`boot_ci()`，`n_boot = 1000`，`set.seed(20230105)`），对齐 `sns.barplot` 默认 |
| 配色 | `age_colors = ['cornflowerblue','thistle','lightcoral']` → `#6495ED / #D8BFD8 / #F08080` |
| inset | 单张，`order = ['Vlmc','Epen','Macro','T cell']` |

**本版修正的三个问题**（上一版与原文不一致的根因）：

1. **密度分母口径错了** —— 上一版用真实凸包**面积**；原文用的是
   `ConvexHull.area`，而 scipy 在二维返回的是**周长**。这一条把 ExN 从
   373.9 拉到 331.5、Olig 从 133 拉到 112，并改变了 Olig 的年龄趋势
   （旧口径 4wk→24wk 几乎持平，原文口径单调上升）。
2. **误差线类型错了** —— 上一版用 t 分布 CI；`sns.barplot` 默认是
   **bootstrap** 95% CI。n = 4/3/5 时两者差别明显，直接把 InN
   （p = 0.033 → FDR 0.107）等从"显著"变成"不显著"。
3. **Epen 被整类丢掉** —— 上一版读的是 `cell_type_annot`，CELLxGENE 在这个
   字段里把室管膜细胞记为 `cell`，映射后为 NA 被过滤，于是 Epen 全为 0。
   改为读 `cell_type` 再做本体名 → 简写映射，并按 `clust_annot == 'Epen'`
   补回（原文也是这样从 Astro 里拆出来的），Epen 恢复为
   4wk 1.78 / 24wk 1.81 / 90wk 7.21。

**修正后的结果与原文一致**（`tables/fig2A_significance.csv`）：

| 细胞类型 | 4wk | 24wk | 90wk | p | BH FDR |
|---|---|---|---|---|---|
| ExN | 331.5 | 293.9 | 295.3 | 0.650 | 0.650 |
| InN | 91.3 | 80.8 | 64.8 | 0.033 | 0.107 |
| MSN | 205.4 | 160.0 | 146.5 | 0.126 | 0.244 |
| **Olig** | 112.2 | 173.5 | 276.4 | 0.0014 | **0.0172** |
| Astro | 114.3 | 88.3 | 106.9 | 0.589 | 0.638 |
| Endo | 102.1 | 60.4 | 85.9 | 0.290 | 0.343 |
| Micro | 57.3 | 36.6 | 44.0 | 0.156 | 0.244 |
| **OPC** | 37.5 | 20.6 | 22.6 | 0.0040 | **0.0172** |
| Peri | 30.2 | 24.2 | 25.8 | 0.169 | 0.244 |
| **VLMC** | 9.6 | 8.0 | 5.7 | 0.0027 | **0.0172** |
| Epen | 1.8 | 1.8 | 7.2 | 0.049 | 0.128 |
| Macro | 3.2 | 4.2 | 1.7 | 0.145 | 0.244 |
| T cell | 0.4 | 0.4 | 0.9 | 0.245 | 0.319 |

**显著类型 = Olig / OPC / VLMC 三类，与原文图上的 3 个 `*` 完全一致。**
（上一版标出 7 个 `*`，是口径错误造成的假阳性。）

- 显著性横线按 dodge 后的实际 x 坐标手动指定：
  `xmin = as.numeric(cell_type) - 0.2333`、`xmax = + 0.2333`；
- y 轴刻度 **0-100-200-300-400-500**（与原文一致），上限 620 以容纳
  ExN 的 bootstrap CI 上帽；
- 导出：`save_fig(..., cairo = TRUE)` 走 `cairo_pdf`（默认 `pdf()` 设备
  不支持 Arial，会报 `invalid font type`），尺寸 12 × 6 in。

与 Python 版（`../Figure2_py/`，逐行照搬原文代码）的交叉验证：
两版 13 类的均值、p 值、BH FDR 与显著类型集合**完全相同**。

### Fig 2B
- **2×2 排列**（左上 Olig | 右上 Endo | 左下 Micro | 右下 Astro），与原文一致；
- 状态名以对应颜色文字直接标在各分面堆叠柱左侧，不另设图例；
- 年龄组注释（竖排）出现在每一行的柱子上方；
- **0% / 100% 标度用文字放在该行第一列分面的右缘**（即两列之间的空隙处），
  不画坐标轴线；`x = 3.56` 必须落在 `scale_x` limits 之内，
  否则该行会被 scale 丢弃、标签不显示；
- **图中不标显著性星号**，说明文字里也不提 `*`；
  显著性结果只输出到 `tables/fig2B_significance.csv`；
- 与 Fig 2A / 2C 统一带有图题与说明（`plot_annotation`）。

### Fig 2C（已修正：此前有 15 个基因整行消失）

- "基因为行、状态为列"的矩阵版式；**基因行顺序与原文 `genes_to_show` 完全一致**
  （Olig: Ank2, Trim2, Robo1, Spock1, C4b, Neat1, Dgki, Il33；
  Micro: Dpp10, Meg3, Ctss, Chrm3, Trem2, Trim30a, Itgb2, B2m；
  Endo: Ly6a, Apoe, B2m, Rbm20, Tmem209, Xdh, Fmo2, Bgn；
  Astro: Lsamp, Gpc5, Slc1a2, Luzp2, Trpm3, Brinp3, Gfap, C4b），共 30 个唯一基因；
- **本版修正的关键 bug**：判断"某基因能否用 MERFISH 实测值"时必须看
  **`colnames(counts)`**（缓存里只保留了 15 个基因），
  上一版写的是 `obj$genes` —— 那是全部 **374** 个 panel 基因的符号表。
  于是 30 个基因都被误判为"在 panel 内"，panel 外的 15 个基因
  既取不到 MERFISH 值、也不会走 snRNA 分支，**整行从图里消失**。
  修正后：30 个基因全部出现，MERFISH 实测 15 个、snRNA 补充 15 个；
- **基因名用横向文字**（`angle = 0`，右对齐贴分面左缘），便于阅读；
  左侧留白 46 pt 以容纳横向基因名与 y 轴标题；
- 状态列标题只画在最上一行基因上方，并为标题预留顶部空间，避免与基因行标签压叠；
- 因每个基因行量程不同（`free_y`），不逐行画 y 刻度（原文亦不显示），
  也不画坐标轴线；精确数值见 `tables/fig2C_expression_summary.csv`；
- 配色与 Fig 2B 共用同一组 `STATE_COLORS`（**逐块采样自原文**，
  见 `tools/sample_paper_colors.py`）；
- **数据来源（实测判定）**：
  * 在 MERFISH 374 基因 panel 内的 15 个基因 → **MERFISH raw 分子计数的 sqrt**。
    这条路径能完整重现原文 Fig 2C 的全部方向性模式（C4b↑Olig-3、Il33↑Olig-3、
    Trem2↑Micro-3、B2m↑Micro-3、Xdh/Fmo2↑Endo-3、Gfap↑Astro-2、C4b↑Astro-2）；
  * panel 外的 15 个基因 → 同状态（`clust_annot`）snRNA-seq 细胞实测
    `sqrt(CP10K)`，基因名后带 `*`。
- `BrainAgingSpatialAtlas_Imputed.h5ad` 经实测**不能**用于作图：虽然它含全部
  7,276 个基因（30 个目标基因都在），但 C4b / Il33 / Trem2 / Gfap 等
  在同状态细胞上的均值**几乎全是 0**（该矩阵是 Z-score 后的 kNN 平均，
  非正值被截断），小提琴会退化成一条线。详见报告"与原文对比"一节。

### 整图组合 (`05_combine_figure.R`)
- 版面：**上半 2A（左）| 2B（右），下半 2C 占满整行**；
- 组合方式用 `patchworkGrob()` + grid 视口手动拼版
  （直接嵌套 patchwork 会触发 `wrap_dims` 报错）。
