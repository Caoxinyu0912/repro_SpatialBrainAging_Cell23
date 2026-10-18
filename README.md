# Figure 5E / 5F 复现（R 版）

复现 Allen et al., *Cell* 186:194–208 (2023) 的 **Figure 5E 与 5F**，
忠实于原文的定义表达（STAR Methods + 作者仓库
[ZhuangLab/SpatialBrainAgingCell22](https://github.com/ZhuangLab/SpatialBrainAgingCell22)
的 `int_analysis.ipynb`）。

## 一键运行

```r
# RStudio: Source run_all.R
# 命令行:
Rscript run_all.R          # 约 20 秒 (缓存命中时)
```

## 文件

| 文件 | 作用 |
|---|---|
| `00_common.R` | h5ad 读取、`score_genes` 复刻、KDTree/分箱、出图主题 |
| `10_fig5E_distance_curve.R` | **Figure 5E**：激活分数 vs 距离（2×3） |
| `20_fig5F_correlation.R` | **Figure 5F**：3 个二维核密度相关图 |
| `run_all.R` | 一键执行 |

## 原文定义（逐条对应）

> 本节每一项都标注了**原文代码出处**。凡原文没有给出、由本复现推断的部分，
> 单独列在「推断与限制」一节，不与原文定义混同。

### 激活分数 —— `sc.tl.score_genes`（原文 CELL 41）

作者 `int_analysis.ipynb` CELL 41：

```python
sc.tl.score_genes(adata_annot, gene_list=['B2m','Trem2','Ccl2','Apoe','Axl',
                  'Itgax','Cd9','C1qa','C1qc','Lyz2','Ctss'],
                  score_name='activate_micro', use_raw=False)
sc.tl.score_genes(adata_annot, gene_list=['C4b','C3','Serpina3n','Cxcl10',
                  'Gfap','Vim','Il18','Hif3a'], score_name='activate_astro')
sc.tl.score_genes(adata_annot, gene_list=activate_endo, score_name='activate_endo')
```

`score_genes` 的定义（scanpy 文档原文）：

> The score is the average expression of a set of genes after subtraction by the
> average expression of a reference set of genes. The reference set is randomly
> sampled from the `gene_pool` for each binned expression value.

默认参数 `ctrl_size = 50`、`n_bins = 25`、`random_state = 0`、
`gene_pool = 全部基因`。R 版 `score_genes()` 逐条实现该定义：
对照基因按全 panel（374 基因）的平均表达等频分 25 箱，每箱随机抽 50 个（`set.seed(0)`），
分数 = 目标基因平均表达 − 对照基因平均表达。

少突胶质炎症分数用 `activate_olig = C4b, Il33, Il18`（STAR Methods：
"The genes used for calculating oligodendrocyte inﬂammation score were:
C4b, Il33, and Il18."）。

中心化（原文 CELL 42/43/471）：

```python
# CELL 42
adata_annot.obs.activate_micro = adata_annot.obs.activate_micro - \
    np.mean(adata_micro[adata_micro.obs.age=='4wk'].obs.activate_micro)
# CELL 43  (同上, 用 activate_astro / adata_astro)
# CELL 471 (activate_olig)
adata_olig = adata_imputed[adata_imputed.obs.cell_type=="Olig"]
adata_imputed.obs.activate_olig = adata_imputed.obs.activate_olig - \
    np.mean(adata_olig[adata_olig.obs.age=='4wk'].obs.activate_olig)
```

三者都是减去**该细胞类型 4wk 组**的均值，R 版一致。

### Figure 5E —— 距离-激活曲线（原文 CELL 46/47 + 464 + 467–470）

```python
# CELL 46
def identify_nearest_neighbors_with_dist(X, Y, min_dist=0):
    dists, ind = KDTree(Y).query(X, k=2)
    # 若最近距离为 0 (自身) 则取第二近邻
def compute_binned_values(dists, scores, min_d=0, max_d=100, bin_size=30):
    for i in np.arange(min_d, max_d-bin_size):
        idx = np.argwhere(np.logical_and(dists > i, dists <= (i+bin_size)))
        binned_mean[i] = np.mean(scores[idx])
        binned_std[i]  = np.std(scores[idx]) / np.sqrt(len(scores[idx]))
    binned_mean -= binned_mean.mean()

# CELL 464
bin_size = 50

# CELL 467
f = plt.figure(figsize=(3,3))
celltypes = ["Endo","Vlmc", "Olig"]
for i in celltypes:
    scores, dists = compute_celltype_obs_distance_correlation(
        adata_imputed[adata_imputed.obs.age=='4wk'], "Micro", i, "activate_micro")
    binned_mean, binned_std = compute_binned_values(dists, scores,
                                                    bin_size=bin_size, max_d=100)
    x = np.arange(len(binned_mean))+bin_size
    plt.plot(x, binned_mean, color=celltype_colors[i])
    plt.fill_between(x, binned_mean-binned_std, binned_mean+binned_std,
                     alpha=0.1, color=celltype_colors[i])
plt.ylim([-0.05, 0.12])
sns.despine()
```

CELL 468/469/470 依次是 Micro-90wk、Astro-4wk、Astro-90wk；
Astro 用 `plt.ylim([-0.2, 0.5])`。

R 版对应实现见 `nearest_dist()` 与 `binned_curve()`，包括作者实现中
`index = 1` 那一格永远为 0 的特征也一并保留。

### Figure 5F —— 空间相关（原文 CELL 46 + 472–474）

```python
# CELL 46
def compute_celltype_obs_correlation(A, cell_type_X, cell_type_Y,
                                     obs_key_X, obs_key_Y,
                                     celltype_key='cell_type', radius=40,
                                     min_dist_thresh=15):
    X = A[A.obs[celltype_key] == cell_type_X]
    Y = A[A.obs[celltype_key] == cell_type_Y]
    curr_X, curr_Y = X.obsm['spatial'], Y.obsm['spatial']
    neighbors_X, ind_X = identify_nearest_neighbors_with_idx(
        curr_X, curr_Y, dist_thresh=radius, min_dist_thresh=min_dist_thresh)
    curr_expr = obs_Y[neighbors_X]
    return obs_X.values[ind_X], curr_expr.values

# CELL 472
x, y = compute_celltype_obs_correlation(
    adata_imputed[adata_imputed.obs.spatial_clust_annots == "CC"],
    "Olig", "Micro", "activate_olig", "activate_micro", radius=50)
R = np.corrcoef(x, y)[0, 1]
sns.kdeplot(x=x, y=y, fill=True)
```

CELL 473 是 `Olig→Astro`（`activate_olig` / `activate_astro`），
CELL 474 是 `Astro→Micro`（`activate_astro` / `activate_micro`）。

两个关键细节（**本复现早期版本曾理解错，导致 R 值偏低**）：
1. **逐邻居配对**——返回的是每个邻居的 `obs` 值，**不是**取 50 µm 内的均值；
2. **`min_dist_thresh=15`**——`identify_nearest_neighbors_with_idx` 中
   `ind = ind[dists > min_dist_thresh]`，即**剔除距离 ≤ 15 µm 的邻居**。

R 版 `pair_neighbors()` 严格照此实现。

## 坐标轴：两套输出

原图与仓库代码的坐标轴**并不一致**，因此 5E 输出两张图。

### ① 主轴 `fig5E_distance_curve_R`（按发表原图）

| 项 | 设置 |
|---|---|
| x 刻度 | **20 与 80**（上下两排共用，只在下排标注） |
| x 范围 | 10 ~ 92（使 20/80 落在与原图相同的位置） |
| x 轴标题 | `Distance to specified cells (µm)`，只在下排 |
| 上排 y | **−0.2 ~ 0.3**，刻度间隔 0.1（Astro activation score） |
| 下排 y | **−0.05 ~ 0.15**，刻度间隔 0.05（Micro activation score） |
| 分箱 | `bin_size = 20`、`max_d = 90` → x = 20…89 |

曲线在面板边界处被**裁切**（原图同样如此）。用
`coord_cartesian(xlim=, ylim=, expand = FALSE)` 而不是 `scale_*(limits=)`：
后者的 `oob = censor` 会把超范围值置为 NA，导致曲线中间`断开`。

### ② 对照 `fig5E_authorstyle_R`（按仓库代码）

| 项 | 设置 |
|---|---|
| x 范围 | 47.5 ~ 101.5，刻度 60 / 80 / 100 |
| 上排 y | −0.2 ~ 0.5（CELL 469/470 的 `plt.ylim`） |
| 下排 y | −0.05 ~ 0.12（CELL 467/468 的 `plt.ylim`） |
| 分箱 | `bin_size = 50`、`max_d = 100` → x = 50…99（CELL 464） |

这张与 `Figure5_authorcode/figures/fig5E_authorcode.png` **逐点一致**：
900 个曲线点 `max|Δmean| = 7.2 × 10⁻¹⁶`、相关系数 1.000000
（见 `compare_with_author.py`）。

### 为什么主轴的 bin_size 是 20 而不是仓库里的 50

`x = arange(max_d - bin_size) + bin_size`，所以 **bin_size 同时决定 x 的起点**：
`bin_size=50` → x 只能从 50 起，与原图 x 轴（从约 10 起、刻度 20/80）矛盾。
主轴取 `bin_size = 20`：
- x 从 20 起，与刻度对齐；
- 每个滑窗细胞数充足（Micro 行的 VLMC 最少 30 个、Astro 行 VLMC 最少 221 个），
  曲线连续、平滑度与原图相当；
- `bin_size=10` 虽然 x 从 10 起更贴近原图，但 Micro 行的 VLMC 每窗仅 5~70 个
  邻居，曲线尖刺过强，反而不像原图。

可用环境变量切换：`FIG5_BIN_PAPER=10`（更细）、`FIG5_BIN_SIZE=50`（改对照版窗宽）。

### 拼接对齐

上下排 y 刻度标签宽度不同，`patchwork` 的 `/` 拼会使两排绘图区左右边界错位；
本复现用 `cowplot::plot_grid(..., align = "v", axis = "lr")` 按 panel 边界对齐。

## 配色（原文 CELL 203/204）

原文用运行时生成的调色板：

```python
# CELL 203
clust_cell_types_10x = { ..., "Vlmc": "olive", "Endo": "khaki",
                              "Olig": "slategrey", ... }
def gen_dark_palette(prefix, color_name, uniq_clusts):
    n = np.sum([1 if prefix in i else 0 for i in uniq_clusts])
    return sns.dark_palette(color_name, n_colors=n+2)[2:]
# CELL 204
celltype_colors, celltype_pals, label_colors, clust_pals = generate_palettes(adata_imputed)
# CELL 467-470
plt.plot(x, binned_mean, color=celltype_colors[i])
```

`celltype_colors` 的最终取值依赖运行时的 `uniq_celltypes` 顺序，无法从发布数据
确定性复算。因此 R 版直接取**原文图上实测**的颜色，
并确认它们与上面的色系一致：

| 曲线 | 实测色 | 与原文色系的关系 |
|---|---|---|
| Endo | `#E8E088` | khaki 系（`sns.dark_palette("khaki")[2:]` 的末位为 `#F0E68C`） |
| VLMC | `#F8D40C` | olive/gold 系的亮色变体 |
| Oligo | `#708090` | 与 matplotlib `"slategrey"` 精确一致 |

### Figure 5F 的填充色（已与原文一致）

原文 5F 只写了 `sns.kdeplot(x=x, y=y, fill=True)`，**没有指定 cmap**，
所以发表图里的蓝色就是 seaborn 2D-KDE 的默认填充梯度。
本复现把该 QuadContourSet 的 cmap（256 级）导出，取 10 个等距代表色写入
`00_common.R` 的 `KDE_FILL_COLORS`：

```r
KDE_FILL_COLORS <- c("#acc9ee", "#90bef1", "#6fb2f5", "#3fa8f7", "#2c9be9",
                     "#2a8fd6", "#2f82c2", "#3676ac", "#3b6995", "#3d5d7f")
```

核对方式：用 `sns.kdeplot(fill=True)` 渲染后采样，主色为
`#9CC0F0 / #78B4F0 / #48A8F6 / #2A96E4 / #2A8AD2 / #307EBA / #366CA2 ...`，
与原文图采样的 `#9CC0EA / #78B4F0 / #48A8F6 / #2A96E4 / #2A8AD2 ...` **逐色吻合**。

R 侧实现要点（`20_fig5F_correlation.R`）：

- `geom_density_2d_filled(contour_var = "ndensity")` + `scale_fill_manual(KDE_FILL_COLORS)`
- `breaks` 用幂次间隔（`kd_pow = 8`），让大部分层级落在低密度段 → 等高线收得紧
- 最外一层填 **白色**，等价于把极稀疏外围擦掉，形成原文的留白轮廓

## 性能优化

5F 的邻居配对原来是一行一行 `for` 循环 (每个 src 细胞一次)，改成
**逻辑下标一次取全部合法邻居对**：

```r
sel <- nn$nn.dist <= radius & nn$nn.dist > min_dist
pos <- which(sel)                                  # 列优先线性下标
src_row <- ((pos - 1L) %% nrow(Xs)) + 1L           # 对应 src 细胞行号
list(x = sx[src_row], y = ty[nn$nn.index[sel]])
```

关键点：`nn$nn.index[sel]` 与 `nn$nn.dist[sel]` 的元素一一对应，
所以 x 与 y **天然配对**，不需要再单独展开行号。

优化后结果与逐行版本**逐位相同**（R 值 0.401 / 0.455 / 0.408 未变），
5F 耗时从数秒降到 ~1 秒（全流程 0.3 分钟）。

> 期间踩过一个坑：`as.vector(t(sel))` 得到的是"按行遍历"掩码，
> 但 `as.vector(nn$nn.index)` 是"按列遍历"顺序，两者错位会让 x/y 配错、
> R 掉到 0.015。改用 `which(sel)` + 取模求行号即可。

## 复现结果

### Figure 5E

曲线形态与原文一致：Olig（灰蓝）在 Juvenile/Young Adult 随距离下降、
Old 在近距离最高；VLMC（黄）单调下降并穿零；Endo（浅黄）相对平坦。

### Figure 5F

| 面板 | R 版 R（radius=30） | R 版 R（radius=50） | 原文 R | 差值(30) |
|---|---|---|---|---|
| Olig → Micro | **+0.401** | +0.406 | 0.46 | −0.059 |
| Olig → Astro | **+0.455** | +0.439 | 0.47 | −0.015 |
| Micro → Astro | **+0.408** | +0.386 | 0.46 | −0.052 |

区域用发布版的 `obs$tissue == "corpus callosum"`，配对用原文的
"逐邻居 + 剔除 ≤15 µm"。改用真实区域后 R 值从 0.36/0.41/0.21
提升到 **0.40/0.46/0.41**，全部进入原文量级。

脚本默认用 **radius = 30 µm**（原文图注），可用 `FIG5_CORR_RADIUS=50` 切到
作者代码里显式传入的 50 µm；两者差异见 I6。

## 推断与限制（非原文直接给出）

### I1. `bin_size`：仓库现值 50 无法产生原文图

作者 notebook 里 `bin_size` **只在 CELL 464 被设为 50**，CELL 467–470 直接沿用。
但代入 `compute_binned_values`（窗口 `i` 覆盖 `(i, i+bin_size]`）：

| bin_size | 窗口数 `max_d-min_d-bin_size` | 能取到值的窗口 | 每格细胞数 | 图上表现 |
|---|---|---|---|---|
| **50**（仓库现值） | 50 | 仅 `i=0,1`（`i≥2` 时 `i+50>100`） | 很大 | 只有 x=50、51，塌成一小段 |
| **10** | 90 | ~90 | 约 30~130 | 锯齿剧烈（`sd(diff)` 比 30 大一个量级），趋势方向不稳 |
| **30（本复现采用）** | 70 | 70，覆盖 30~99 µm | 约 200~1700 | 平滑且趋势正确，保留原文那种轻微起伏 |

量化对比（Astro/4wk/Oligo，`sd(diff)` 越小越平滑）：

| bin_size | 5 | 10 | 20 | **30** | 50 |
|---|---|---|---|---|---|
| `sd(diff)` | 0.0479 | 0.0142 | 0.0087 | **0.0057** | 0.0035 |
| `range` | 0.835 | 0.443 | 0.373 | **0.251** | 0.140 |

`bin_size=50` 与 `30` 的结果趋势一致、都平滑；`10` 会引入不真实的噪声
（相邻窗口重叠 9 µm，每格样本量太小）。三者都实现了作者的同一个公式，
差别只在窗口宽度。**本复现取 30 µm 属推断，不是原文代码里的值。**

```bash
FIG5_BIN_SIZE=50 Rscript run_all.R     # 复现"当前仓库代码"的行为
FIG5_BIN_SIZE=10 Rscript run_all.R     # 细窗对照 (可见锯齿)
FIG5_BIN_SIZE=40 Rscript run_all.R     # 邻近取值, 结论不变
```

### I2. 24wk 面板无对应代码

CELL 467–470 只有 `4wk` 与 `90wk` 的调用，**没有 24wk 的 cell**。
但原文 5E 的三列是 Juvenile / Young Adult / Old，
因此 24wk 面板按同样的口径计算（其余参数完全不变）。

### I3. 解剖区域：**已找到发布版的区域注释，不再需要近似**

> 这一项此前被错误地写成"发布数据没有解剖区域注释"。实际上
> **`obs$tissue` 就是区域注释**，取值与原文 CELL 455 的
> `spatial_order = ['Pia','L2/3','L5','L6','CC','LatSept','Striatum','Ventricle']`
> 一一对应：

| `obs$tissue` | 原文名称 | 细胞数 | 其中 Olig |
|---|---|---|---|
| `pia mater` | Pia | 15,834 | 3.9% |
| `cortical layer II/III` | L2/3 | 45,380 | 4.4% |
| `cortical layer V` | L5 | 67,178 | 10.8% |
| `cortical layer VI` | L6 | 72,658 | 18.2% |
| **`corpus callosum`** | **CC** | **53,460** | **66.7%** |
| `olfactory region` | LatSept | 25,996 | 9.0% |
| `striatum` | Striatum | 95,601 | 10.7% |
| `brain ventricle` | Ventricle | 2,811 | 10.1% |

`corpus callosum` 占 14.11% 的细胞、**少突胶质细胞占 66.7%**，是典型白质束，
与胼胝体一致。因此 Fig 5F **直接使用发布版的区域注释**：

```r
in_cc <- meta$tissue == "corpus callosum"     # 53,460 细胞
```

早期版本用"30-NN 中 Olig 比例 ≥ 75 分位"自行近似，会混入皮层细胞，
导致相关系数偏低。改用真实区域后 5F 的 R 值明显改善（见下）。
早期近似下 5F 的三组 R 为 0.36 / 0.41 / 0.21，改用真实区域后为
**0.40 / 0.46 / 0.41**（radius=30）。

### I4. 表达矩阵尺度不同

作者用他们自己流程产出的 `adata_imputed.X`；发布版的 X 是 **Z-score 后的
kNN 平均**（73% 元素为负），因此本复现改用 **raw 分子计数**并显式传入
全部 374 个 panel 基因作为对照基因池。分数的绝对尺度因此与原文不同，
但年龄趋势与相关强度一致。

### I5. 24wk 组 Micro 激活分数为负

| 细胞类型 | 4wk | 24wk | 90wk |
|---|---|---|---|
| Micro | 0.000 | **−0.363** | +0.822 |
| Astro | 0.000 | −0.038 | +0.916 |
| Olig | 0.000 | +0.541 | +2.229 |

与原图 Figure 5A「随年龄升高」不完全一致。**这是数据本身的性质**：
不做任何对照扣除、直接看原始平均计数，24wk 的 Micro 同样最低
（4wk 1.353 / 24wk 0.828 / 90wk 2.183 分子/细胞），11 个激活基因全部在 24wk 最低，
三只 24wk 供体也一致偏低。改口径无法消除，已在报告中如实标注。

### I6. 5F 的 radius：图注 30 µm vs 代码 50 µm

原文 Fig 5F 图注写 *"within 30 µm"*（PDF 文字层把 µm 排成了 mm），
但作者代码 CELL 472–474 **显式传入 `radius=50`**。两者都实现出来核对：

| radius | Olig→Micro | Olig→Astro | Micro→Astro |
|---|---|---|---|
| 30 µm（图注，本复现默认） | 0.401 | 0.455 | 0.408 |
| 50 µm（代码） | 0.406 | 0.439 | 0.386 |

**差异很小（|ΔR| ≤ 0.022），结论不随 radius 改变。** 默认取图注的 30 µm
（与发表文字一致），可用 `FIG5_CORR_RADIUS=50` 切到代码值。

## 输出

- `figures/fig5E_distance_curve_R.{png,pdf}` — 2×3 距离-激活曲线
- `figures/fig5F_correlation_R.{png,pdf}` — 3 个二维核密度相关图
- `figures_hires/` — 900 dpi 版本（供报告排版）
- `tables/fig5E_distance_curve_R.csv` — 每条曲线逐点 mean/SEM
- `tables/fig5E_distance_summary_R.csv` — 各面板细胞数与距离统计
- `tables/fig5F_correlation_R.csv` — 三个面板的 R 与配对数
- `tables/fig5F_pairs_R.csv` — 5F 的逐配对原始值（用于复核 R 值）
- `tables/fig5A_activation_check_R.csv` — 各年龄组激活分数（核对 5A 趋势）
- `cache/` — 读取 h5ad 后的缓存（删除后自动重建）

## 环境

R 4.6.1 + hdf5r 1.3.16 / Matrix / dplyr / tidyr / ggplot2 / patchwork /
FNN / RANN / scales

## 数据

CELLxGENE Discover collection `31937775-0602-4e52-a799-b6acdd2bac2e`：
`BrainAgingSpatialAtlas_MERFISH.h5ad`（378,918 细胞 × 374 基因，
含 raw 计数与 µm 空间坐标）。
