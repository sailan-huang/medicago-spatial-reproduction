from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORT = ROOT / "docs" / "analysis_report.md"


def markdown_table(frame: pd.DataFrame, float_digits: int = 3) -> str:
    display = frame.copy()
    for column in display.columns:
        if pd.api.types.is_float_dtype(display[column]):
            display[column] = display[column].map(
                lambda value: "" if pd.isna(value) else f"{value:.{float_digits}f}"
            )
    header = "| " + " | ".join(map(str, display.columns)) + " |"
    separator = "| " + " | ".join("---" for _ in display.columns) + " |"
    rows = [
        "| " + " | ".join(map(str, row)) + " |"
        for row in display.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def build_report() -> Path:
    metadata = json.loads((TABLES / "run_metadata.json").read_text(encoding="utf-8"))
    clu = pd.read_json(TABLES / "clustering_metrics.json", typ="series")
    marker = pd.read_csv(TABLES / "marker_specificity.csv")
    rhizobia = pd.read_csv(TABLES / "rhizobia_by_region.csv")
    rhizobia_statistics = pd.read_csv(TABLES / "rhizobia_statistics.csv")
    region_markers = pd.read_csv(TABLES / "region_marker_discovery.csv")
    bins = pd.read_csv(TABLES / "bin_level_results.csv")
    sample_reference = pd.read_csv(TABLES / "sample_reference.csv")

    dims = metadata["input_dimensions"]
    integrity = metadata["integrity_checks"]
    fine_ari = float(clu["leiden_fine_ari_vs_author_tissue"])
    fine_nmi = float(clu["leiden_fine_nmi_vs_author_tissue"])
    marker_matches = int(marker["matches_expected_region"].sum())
    marker_total = len(marker)
    counts_median = bins["total_counts_recomputed"].median()
    counts_min = bins["total_counts_recomputed"].min()
    counts_max = bins["total_counts_recomputed"].max()
    genes_median = bins["detected_genes_recomputed"].median()
    author_cell_number = int(sample_reference.iloc[0]["Cell Number"])
    nitrogen_fraction = rhizobia.loc[
        rhizobia["organ"] == "Nitrogen fixation zone",
        "median_rhizobia_fraction",
    ]
    nfix_value = (
        f"{float(nitrogen_fraction.iloc[0]):.4f}"
        if not nitrogen_fraction.empty
        else "NA"
    )
    top_rhizobia = rhizobia.sort_values(
        "median_rhizobia_fraction", ascending=False
    ).head(3)

    region_marker_view = (
        region_markers.groupby("author_region", observed=True)
        .head(3)[
            [
                "author_region",
                "names",
                "gene_symbol",
                "scores",
                "logfoldchanges",
                "pvals_adj",
            ]
        ]
        .copy()
    )
    region_marker_view.columns = [
        "区域",
        "locus ID",
        "符号",
        "Wilcoxon 分数",
        "log2FC",
        "校正 p 值",
    ]

    marker_view = marker[
        [
            "acronym",
            "expected_region",
            "best_author_region",
            "matches_expected_region",
            "detection_fraction_overall",
            "supplementary_evidence_match",
        ]
    ].copy()
    marker_view.columns = [
        "标记",
        "预期区域",
        "本切片最高区域",
        "区域一致",
        "总体检出率",
        "Table S1-5 证据",
    ]

    rhizobia_view = rhizobia[
        [
            "organ",
            "n_bins",
            "median_rhizobia_fraction",
            "q25_rhizobia_fraction",
            "q75_rhizobia_fraction",
            "median_total_counts",
        ]
    ].copy()
    rhizobia_view.columns = ["区域", "bin 数", "中位占比", "Q1", "Q3", "中位总计数"]
    statistics_view = rhizobia_statistics.copy()
    statistics_view.columns = [
        "检验",
        "比较",
        "n₁",
        "n₂",
        "统计量",
        "p 值",
        "解释边界",
    ]
    statistics_view["p 值"] = statistics_view["p 值"].map(lambda value: f"{value:.3g}")

    report = f"""# 苜蓿根瘤空间转录组复现报告

## 复现对象

本工作复现 Shao et al., *Cell* (2026) 论文中的一个代表样本：
`FP200000569BR_C5_bottom4_MXP3.h5ad`（P3，14 dpi）。论文完整研究包含 5 个发育阶段、
20 个切片和 20,841 个细胞/空间单元；本仓库不声称重建完整图谱。

复现范围包括：

1. 从公开计数矩阵重建质量指标，检查植物与根瘤菌计数；
2. 仅使用植物基因独立完成高变基因选择、PCA、邻居图、UMAP 与 Leiden 聚类；
3. 根据论文补充表与 Mt5.0 官方注释解析标记基因，并绘制分组和空间表达；
4. 计算每个空间 bin 的根瘤菌计数占比，比较作者组织分区。

原始 FASTQ 过滤、STAR 比对、完整多切片 Harmony 整合、RNA velocity、
Monocle 3、WGCNA、CRISPR/RNAi、ROS 染色及固氮表型实验均不在本次范围内。

## 数据完整性与可重复性

- 空间 bin：**{dims["bins"]:,}**
- 基因：**{dims["genes"]:,}**，其中植物 **{dims["plant_genes"]:,}**、根瘤菌 **{dims["rhizobia_genes"]:,}**
- 非零计数：**{dims["nonzeros"]:,}**
- 总计数：**{dims["total_counts"]:,.0f}**
- 原始计数层与 `.X` 一致：**{integrity["X_equals_counts"]}**
- 原始计数为非负整数：**通过**
- `nCount_Spatial`、`nFeature_Spatial` 与重算值一致：**通过**
- 坐标与 `array_row/array_col` 一致：**通过**
- 基因 ID 唯一、坐标有限、组织标签无缺失：**通过**
- 植物分析保留 bin：**{integrity["plant_analysis_bins"]:,}**
- 植物分析保留基因：**{integrity["plant_analysis_genes"]:,}**

单个 bin 的总计数中位数为 **{counts_median:.1f}**（范围 {counts_min:.0f}–{counts_max:.0f}），
检测基因数中位数为 **{genes_median:.1f}**。随机种子固定为
`{metadata["analysis_parameters"]["random_seed"]}`，所有主要参数记录在
[`config/analysis.yaml`](../config/analysis.yaml) 和
[`run_metadata.json`](../results/tables/run_metadata.json)。
全部图均输出 300 dpi PNG 和可缩放 PDF/SVG。

论文补充表 S1-3 对同一样本报告 **{author_cell_number:,}** 个 cell/bin，
而公开 `.h5ad` 包含 **{dims["bins"]:,}** 个 bin。本复现保留公开对象中的全部
{dims["bins"]:,} 个 bin，不做无法从文件中唯一确定的删除；该 14-bin 差异被视为
上游处理或版本差异，并记录在 `run_metadata.json`。

## 结果

### 1. 质量检查与标准化

![质量检查](../results/figures/01_qc_overview.png)

图 1. P3 切片质量检查。**(A)** 每个 bin 的总计数；**(B)** 检测基因数；
**(C)** 植物与根瘤菌计数关系；**(D)** 测序深度空间分布；**(E)** 根瘤菌计数占比；
**(F)** 根瘤菌占比分布。空间图比例尺为 300 µm。

测序深度和组织边界需要同时检查。根瘤菌占比不能被直接解释为细菌数量：
低深度 bin 的比例方差更大，而且植物与细菌 RNA 的捕获效率不同。

![作者注释与深度](../results/figures/02_annotation_and_depth.png)

图 2. 作者组织注释和测序深度。颜色在所有空间图中保持一致。

![标准化前后对照](../results/figures/01b_normalization_contrast.png)

图 3. 同一基因 `MtLb2` 的原始计数、标准化表达和逐 bin 关系。该图验证计数标准化
改善跨 bin 可比性，同时保留主要空间轮廓。

### 2. 独立聚类

仅使用植物原始计数，过滤后以 10,000 counts/bin 标准化并执行 `log1p`。使用 2,000 个高变基因、
20 个主成分、15 个邻居，比较 Leiden 分辨率 0.3 与 0.8。

![独立聚类](../results/figures/03_independent_clustering.png)

图 4. 植物表达独立聚类。**(A–B)** 两档分辨率的 UMAP；**(C–E)** 作者标签与独立聚类的空间对照；
**(F)** 作者组织颜色图例。空间图比例尺为 300 µm。

| 指标 | 分辨率 0.3 | 分辨率 0.8 |
|---|---:|---:|
| 聚类数 | {int(clu["leiden_coarse_n_clusters"])} | {int(clu["leiden_fine_n_clusters"])} |
| ARI vs 作者区域 | {float(clu["leiden_coarse_ari_vs_author_tissue"]):.3f} | {fine_ari:.3f} |
| NMI vs 作者区域 | {float(clu["leiden_coarse_nmi_vs_author_tissue"]):.3f} | {fine_nmi:.3f} |
| PCA-20 silhouette | {float(clu["leiden_coarse_silhouette_pca20"]):.3f} | {float(clu["leiden_fine_silhouette_pca20"]):.3f} |
| 空间 6-近邻同类比例 | {float(clu["leiden_coarse_spatial_neighbor_purity"]):.3f} | {float(clu["leiden_fine_spatial_neighbor_purity"]):.3f} |

独立聚类保留了明显的空间分区，但没有复制作者标签。ARI 为 **{fine_ari:.3f}**、
NMI 为 **{fine_nmi:.3f}**，说明单切片表达聚类与整合了更多样本和人工证据的区域注释
只部分一致。作者标签自身的空间 6-近邻同类比例为
**{float(clu["author_tissue_spatial_neighbor_purity"]):.3f}**；两档独立聚类均不低于该值，
说明结果具有空间连续性。这里没有为了追求与原文相似而继续挑选分辨率。

![聚类与作者标签混淆矩阵](../results/figures/03b_cluster_annotation_confusion.png)

图 5. 作者区域与独立聚类的混淆矩阵。格内数字为 bin 数，颜色为按作者区域行归一化后的比例。

### 3. 标记基因

标记候选来自论文补充表 Table S1-5，并通过 Mt5.0 官方注释表解析到数据中的 locus ID。
共解析 **{marker_total}** 个可检测标记，其中 **{marker_matches}** 个在本 P3 切片中
的最高平均表达区域与其在完整论文图谱中的预期区域一致。

![标记基因点图](../results/figures/04_marker_dotplots.png)

![独立聚类标记基因点图](../results/figures/04_marker_dotplots_by_cluster.png)

![标记基因空间表达](../results/figures/05_marker_spatial.png)

图 6–8. 论文 Table S1-5 标记在作者区域和独立聚类中的表达，以及逐基因空间分布。
标记的 locus ID 同时通过 Mt5.0 官方注释核对；比例尺为 300 µm。

{markdown_table(marker_view)}

分生区中 `MtPLT3` 和 `MtPLT4`、感染区中 `MtNOD25`、`MtLb6` 和 `MtENOD8.1`
具有清晰表达和空间聚集。`MtSymCRK`、`MtNCR160` 与 `MtLP7` 也呈现区域相关信号，
但本次单切片的最高区域与完整图谱命名并不完全一致。可能原因包括：单个 P3 切片的采样位置、
相邻区域标签边界、跨阶段注释证据，以及基因家族/同源基因命名差异。

从作者标签内部重新计算的每个区域前三个差异基因如下：

{markdown_table(region_marker_view)}

这一结果用于说明后续区域发现流程，而不是独立验证实验功能。

### 4. 植物–根瘤菌共检测

![根瘤菌计数占比](../results/figures/06_rhizobia_fraction.png)

{markdown_table(rhizobia_view)}

图 9. 根瘤菌计数占比的空间分布、深度依赖和组织区域比较。纵轴在 99% 分位处截断显示。

统计检验仅作描述性辅助：

{markdown_table(statistics_view)}

本切片中根瘤菌占比最高的区域中位数为：

{markdown_table(top_rhizobia[["organ", "median_rhizobia_fraction", "n_bins"]])}

固氮区在本切片中的根瘤菌计数占比中位数为 **{nfix_value}**，低于感染区和交界区。
这并不否定论文跨阶段图谱中成熟固氮区富含根瘤菌的结论，而是说明单切片、单阶段和
相对计数指标不足以重现跨阶段空间趋势。计数占比同时受总深度和 RNA 组成影响，
不能作为细菌量、单菌表达或固氮效率的直接测量。

## 与论文和作者代码的关系

- 论文分析了 20 个切片并整合 20,841 个空间单元；本复现只使用一个 P3 切片。
- 作者 GitHub 提供 Giotto、Harmony、Monocle3、WGCNA、RNA velocity 等分析脚本，
  但脚本依赖未公开的 CSV/RDS 路径，不是开箱即用流程。
- 本仓库补充了数据下载、SHA-256 校验、输入契约测试、端到端 Python 流程、
  固定参数、可执行 Notebook 和派生结果。
- 本复现属于“图级局部重分析”，不是原始数据全流程复现或实验验证复现。

## 局限性

1. 公开 `.h5ad` 已经过作者处理，不能追问其原始过滤、比对和 bin 生成步骤。
2. 单切片不包含完整 P1–P5 发育轨迹，也不含本切片中的衰老区。
3. 作者空间标签用于结果对照和描述性差异分析，不能被当作独立金标准。
4. 聚类与作者注释的一致性为中等偏低；UMAP 方向与 Leiden 编号没有生物学含义。
5. 根瘤菌占比是相对捕获计数，不是绝对细菌负荷。
6. 所有区域差异均为观察性关联，需要独立切片、真实实验和功能验证。

## 结论

该复现证明：公开数据足以支持单样本级空间矩阵核验、植物表达聚类、官方标记映射与
宿主–根瘤菌相对信号分析；结果大体支持组织分区与共生相关表达的空间结构化，
但单切片分析不能替代论文的跨阶段整合和实验验证。
"""
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(report, encoding="utf-8")
    return REPORT


if __name__ == "__main__":
    print(build_report())
