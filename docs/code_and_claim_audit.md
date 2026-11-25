# 原文与作者代码核对记录

## 原文事实

- 论文标题：*Spatiotemporal co-transcriptomics uncovers coordinated rhizobia-legume
  stress programs supporting nitrogen fixation in Medicago*
- 期刊：Cell, 2026
- DOI：10.1016/j.cell.2026.07.033
- 平台：Stereo-seq V1.1 mixed-probe panel
- 原始研究：P1–P5 五个阶段、20 个切片、7 个芯片、20,841 个空间单元/细胞
- 空间 bin：60 × 60 DNB spots，约 30 μm
- 聚类：Seurat 3.2.1、Harmony、2,000 HVG、top 20 PCs、resolution 0.5

## 作者代码审计

作者仓库 `XiangSH111/Legume_Spatial_Atlas` 在 `main` 分支提供 12 个脚本：

- `Giotto.R`
- `Harmony.R`
- `Cell_cluster.R`
- `Plot_UMAP.R`
- `Dotplot.R`
- `AverageExpression.R`
- `RCTD.R`
- `Monocle3.R`
- `RNA_Velocity.py`
- `Mfuzz.R`
- `WGCNA.R`
- `Sankey.R`

核对结论：

1. 仓库没有附带 `data/` 目录，也没有锁定环境的容器或 conda 文件。
2. `Cell_cluster.R` 硬编码了 `BIN60FP200000343BR_C4_5_P4.csv`，不是公开 `.h5ad`。
3. `Plot_UMAP.R`、`Dotplot.R`、`AverageExpression.R` 等需要
   `Medicago_nodule.RDS`，仓库没有提供该对象。
4. `RCTD.R`、`Monocle3.R`、`RNA_Velocity.py` 等引用多个私有输入路径。
5. 因此，作者脚本适合作为参数、可视化方式和分析思路参考，但不能作为本次公开数据的
   一键复现入口。

## 本仓库采用的替代方案

- 使用 Python 3.12、Scanpy 1.12.6 和 Notebook，避免 R 旧版本依赖。
- 直接读取作者公开的 `.h5ad`，保留其 `counts` 层。
- 固定分析参数和随机种子；每个输出都能从原始下载文件重新生成。
- 使用论文 Table S1-5 与 Medicago Mt5.0 官方注释交叉解析标记基因。
- 对作者标签只做结果对照，不作为独立聚类的输入。
- 在报告中明确所有与原文不完全一致的结果。

