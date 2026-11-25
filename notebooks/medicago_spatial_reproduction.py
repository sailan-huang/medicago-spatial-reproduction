# %% [markdown]
# # 苜蓿根瘤空间转录组：单样本局部复现
#
# 本 Notebook 复现 Shao et al., *Cell* (2026) 苜蓿根瘤共转录组论文中的四项核心分析：
#
# 1. 质量检查与标准化；
# 2. 植物基因的独立聚类、UMAP 和空间定位；
# 3. Mt5.0 标记基因映射、点图和空间表达；
# 4. 植物与根瘤菌计数拆分及根瘤菌计数的空间占比。
#
# 分析样本是公开的 P3 切片 `FP200000569BR_C5_bottom4_MXP3.h5ad`。这是论文 20 个切片中的
# 一个代表切片，不是完整图谱重算。原始 FASTQ、比对、文章中的遗传学和湿实验不在范围内。

# %%
from pathlib import Path
import os
import sys

import pandas as pd
from IPython.display import Image, Markdown, display

project_root = Path.cwd()
if project_root.name == "notebooks":
    project_root = project_root.parent
os.chdir(project_root)
sys.path.insert(0, str(project_root))

from analysis.reproduction import (
    Paths,
    add_qc_metrics,
    configure_logging,
    configure_publication_style,
    load_and_validate,
    load_config,
    load_marker_mapping,
    load_sample_reference_metadata,
    load_supplementary_marker_evidence,
    plot_annotation_and_depth,
    plot_cluster_confusion,
    plot_clustering,
    plot_marker_genes,
    plot_normalization_contrast,
    plot_qc_overview,
    plot_rhizobia_analysis,
    prepare_plant_analysis,
    rank_region_markers,
    resolve_inputs,
    resolve_markers,
    resolve_supplementary,
    transfer_analysis_results,
    write_run_metadata,
)

configure_logging()
configure_publication_style()
paths = Paths()
config = load_config()
sample_path, annotation_path = resolve_inputs(config)
supplementary_path = resolve_supplementary(config)
sample_reference = load_sample_reference_metadata(
    supplementary_path,
    config["data"]["sample_filename"].removesuffix(".h5ad"),
)
sample_reference.to_csv(
    paths.table_dir / "sample_reference.csv",
    index=False,
    encoding="utf-8-sig",
)

display(
    Markdown(
        f"""
**输入文件**

- 空间矩阵：`{sample_path.name}`
- Mt5.0 基因注释：`{annotation_path.name}`
- 论文补充表：`{supplementary_path.name}`

三个文件均按 SHA-256 校验，运行前会拒绝被替换或不完整的数据。
"""
    )
)
display(Markdown("作者补充表 S1-3 报告的同一样本信息："))
display(sample_reference)

# %% [markdown]
# ## 1. 输入契约与质量检查
#
# 原始计数层被保留。植物和根瘤菌基因按官方 ID 前缀区分：
# `MtrunA17...` 为植物，`SM-...` 为 *Sinorhizobium meliloti*。

# %%
adata = load_and_validate(sample_path)
adata = add_qc_metrics(adata)
display(adata)
display(adata.obs["organ"].value_counts().rename("bins").to_frame())

integrity = {
    "空间 bin": adata.n_obs,
    "基因数": adata.n_vars,
    "植物基因": int((adata.var["species"] == "Medicago").sum()),
    "根瘤菌基因": int((adata.var["species"] == "Sinorhizobium").sum()),
    "非零计数": int(adata.X.nnz),
    "总计数": int(adata.X.sum()),
}
display(Markdown("\n".join(f"- **{key}**：{value:,}" for key, value in integrity.items())))

plot_qc_overview(adata, paths.figure_dir / "01_qc_overview.png")
display(Image(filename=str(paths.figure_dir / "01_qc_overview.png")))

# %%
plot_annotation_and_depth(
    adata,
    paths.figure_dir / "02_annotation_and_depth.png",
)
display(Image(filename=str(paths.figure_dir / "02_annotation_and_depth.png")))

# %% [markdown]
# ## 2. 独立聚类与空间定位
#
# 聚类只使用植物基因，避免根瘤菌比例直接驱动宿主组织分组。步骤为：
# 过滤低表达基因、10,000 counts/ bin 标准化、`log1p`、2,000 个高变基因、20 个主成分、
# 15 个近邻、Leiden 聚类。随机种子固定为 42。

# %%
plant, plant_mask = prepare_plant_analysis(adata, config)
adata = transfer_analysis_results(adata, plant)

normalization_gene = config["analysis"]["normalization_example"]
plot_normalization_contrast(
    plant,
    normalization_gene["gene_id"],
    normalization_gene["symbol"],
    paths.figure_dir / "01b_normalization_contrast.png",
)
display(Image(filename=str(paths.figure_dir / "01b_normalization_contrast.png")))

metrics = plot_clustering(
    plant,
    config,
    paths.figure_dir / "03_independent_clustering.png",
    paths.table_dir / "clustering_metrics.json",
)
display(Image(filename=str(paths.figure_dir / "03_independent_clustering.png")))
display(Markdown("```json\n" + str(metrics) + "\n```"))

plot_cluster_confusion(
    plant,
    config,
    paths.figure_dir / "03b_cluster_annotation_confusion.png",
)
display(
    Image(
        filename=str(
            paths.figure_dir / "03b_cluster_annotation_confusion.png"
        )
    )
)
display(
    Markdown(
        "除 ARI/NMI 外，本流程还计算 PCA 空间中的 silhouette 分数和空间 6-近邻"
        "同类比例，用于同时评估表达可分性和空间连续性。"
    )
)

discovery = rank_region_markers(
    plant,
    load_marker_mapping(annotation_path),
    paths.table_dir / "region_marker_discovery.csv",
)
display(discovery.groupby("author_region", observed=True).head(5))

# %% [markdown]
# ## 3. 标记基因与区域解释
#
# 标记基因来自论文补充表 Table S1-5，并使用 Mt5.0 官方注释表把基因名解析为数据中的 locus ID。
# 对每个区域保留一个或两个在本次 P3 样本中可检测的标记。

# %%
mapping = load_marker_mapping(annotation_path)
supplementary_evidence = load_supplementary_marker_evidence(supplementary_path)
marker_table = resolve_markers(
    config["marker_genes"],
    config["marker_preferred_loci"],
    mapping,
    set(plant.var_names),
    supplementary_evidence,
)
marker_table.to_csv(
    paths.table_dir / "marker_resolution.csv",
    index=False,
    encoding="utf-8-sig",
)
display(marker_table)

marker_summary = plot_marker_genes(
    plant,
    marker_table,
    config,
    paths.figure_dir / "04_marker_dotplots.png",
    paths.figure_dir / "05_marker_spatial.png",
    paths.table_dir / "marker_specificity.csv",
)
display(marker_summary)
display(Image(filename=str(paths.figure_dir / "04_marker_dotplots.png")))
display(
    Image(
        filename=str(
            paths.figure_dir / "04_marker_dotplots_by_cluster.png"
        )
    )
)
display(Image(filename=str(paths.figure_dir / "05_marker_spatial.png")))

# %% [markdown]
# ## 4. 植物–根瘤菌信号
#
# 每个 bin 的根瘤菌占比定义为：
#
# $$
# \mathrm{rhizobia\ fraction}_i =
# \frac{\mathrm{rhizobia\ counts}_i}
# {\mathrm{plant\ counts}_i + \mathrm{rhizobia\ counts}_i}
# $$
#
# 该指标受捕获效率、RNA 组成和组织密度影响，不能解释为细菌数量或固氮效率。

# %%
rhizobia_summary = plot_rhizobia_analysis(
    adata,
    paths.figure_dir / "06_rhizobia_fraction.png",
    paths.table_dir / "rhizobia_by_region.csv",
    paths.table_dir / "rhizobia_statistics.csv",
)
display(rhizobia_summary)
display(Image(filename=str(paths.figure_dir / "06_rhizobia_fraction.png")))
display(pd.read_csv(paths.table_dir / "rhizobia_statistics.csv"))

# %%
adata.obs.to_csv(
    paths.table_dir / "bin_level_results.csv",
    encoding="utf-8-sig",
)
write_run_metadata(
    adata,
    plant,
    metrics,
    config,
    sample_reference,
    paths.table_dir / "run_metadata.json",
)

# %% [markdown]
# ## 结论边界
#
# - 本复现证明公开矩阵可以被独立读取、重建 QC、重新聚类、映射标记并拆分物种计数。
# - 独立聚类不是论文注释的复制品。分辨率 0.8 的空间结构清晰，但 ARI 只有约 0.20，
#   说明作者注释整合了本次单切片之外的证据；不能通过挑选参数追求图形完全一致。
# - 根瘤菌占比在感染区/交界区较高，在分生区较低。该结果支持“宿主组织分区与细菌信号相关”，
#   但不能推断因果关系。
# - 每个 bin 不是独立生物学重复，因此区域统计只作描述性辅助。
# - 补充表 S1-3 报告 1,270 个 cell/bin，公开矩阵为 1,284 个；本分析保留公开对象
#   的全部 bin，并在元数据中记录该上游差异。
# - 论文完整的 20 个切片、Harmony 整合、RNA velocity、Monocle 3、WGCNA 和实验验证不在此范围内。

