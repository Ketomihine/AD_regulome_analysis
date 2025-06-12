# AD Regulome 分析

该仓库包含用于处理和分析单细胞 ATAC-seq 与 RNA-seq 数据的脚本和工作流程，旨在研究阿尔茨海默病（AD）的调控组。代码按照分析步骤划分在不同目录中，下面依次介绍各目录及其中的文件。

## 根目录文件

- `README.org` – 关于仓库目的的简短 org-mode 说明。
- `.gitmodules` – 定义 `epiclust` 子模块（此处为空）。
- `.gitignore` – 针对 Python 与 R 项目的常规忽略规则。

## `integration_linking_modules`

此目录包含用于峰位点与基因关联、模块推断以及多组学整合的 Python、R 和 shell 脚本。

<details>
<summary>文件列表</summary>

- **aggregate_modules.py** – 聚合各簇的模块计数。
- **auprc.py / auprc_full.py** – 计算预测链接的精确率-召回曲线。
- **coaccessibility.py / coaccessibility.sh** – 合并基因与峰矩阵计算共可及性。
- **concat_bedgraph.sh** – 合并 bedGraph 文件的辅助脚本。
- **count_peaks_in_fragments.sh** – 统计片段文件中的峰覆盖。
- **count_peaks_to_h5ad.py** – 将峰计数转换为 AnnData 对象。
- **deg.R** – 计算模块的差异表达/差异峰。
- **diff_mod.R** – 模块差异分析。
- **epiclust_generate_interactions.py** – 构建 `epiclust` 关联模型的候选峰-基因对。
- **epiclust_linking.py / epiclust_linking.sh** – 在伪汇总数据上运行 `epiclust` 关联模型。
- **epiclust_norm_links.py** – 归一化关联得分并评估性能。
- **epiclust_predict_links.py / epiclust_predict_links.sh** – 使用逻辑回归从训练数据预测关联。
- **filter_tilematrix_to_peaks.py** – 在 TileMatrix 中标注峰。
- **fit_adgwas.R** – 针对关联峰拟合 AD GWAS 富集模型。
- **gene_distance.py** – 计算峰与基因之间距离的工具函数。
- **gene_estimation.py** – 结合可选的关联信息从 ATAC 数据估计基因活性。
- **generate_ct_bedgraph.sh** – 生成各细胞类型的 bigWig/bedGraph 轨迹。
- **go_tf.R** – 对模块进行 GO/TF 富集分析。
- **gtf.py** – 读取 GTF 注释并计算基因长度评分。
- **impute.py** – TF‑IDF 归一化及 NMF 插补工具。
- **integrate.py / integrate.sh** – 使用 Harmony/BBKNN 进行 ATAC 与 RNA 整合并转移标签。
- **linking.py** – 用于距离加权关联的训练与预测框架。
- **load_fragments_from_annot.py / load_fragments_from_annot.sh** – 根据样本注释从片段文件构建 SnapATAC 对象。
- **module_pseudobulk.py** – 创建模块级别的伪汇总计数。
- **module_stats.R / module_stats.py** – 生成模块的统计与汇总。
- **pseudobulk.py / pseudobulk.sh** – 根据整合结果生成伪汇总 AnnData 对象。
- **run_modules.py** – 运行 `epiclust` 模块发现流程。
- **run_seurat_signac_integration.sh** – 预留的 Seurat/Signac 整合脚本。
- **seurat_signac_integration.R** – 实例流程：在 Seurat/Signac 中整合 scRNA 与 scATAC。
- **seurat_signac_multiome.R** – 适用于 10x multiome 数据的 Seurat 脚本。
- **seurat_signac_peaks.R** – 使用 Signac/Cicero 进行峰调用。
- **subtype_enrichment.R** – 比较不同亚型中模块的富集情况。
</details>

## `snATAC.processing`

此目录存放单核 ATAC‑seq 数据的预处理及后续分析脚本。

<details>
<summary>文件列表</summary>

- **1.scATAC.processing.R** – ArchR 主流程，创建项目、执行质控、LSI 降维及初步聚类。
- **2.RemoveDoublets.iter1.R** – 首次去重并重新聚类，包含大量 QC 图。
- **3.RemoveDoublets.iter2.R** – 第二轮去重与聚类。
- **4.callpeak.and.GREAT_annotation.R** – 使用 MACS2 调峰，构建峰矩阵并通过 rGREAT 注释。
- **5.co_accessibility.R** – 利用 ArchR 计算共可及性环。
- **6.TF.candidate.R** – 基于基序富集与基因得分识别候选转录因子。
- **6_2.TF_footprinting.R** – 对候选转录因子进行足迹分析。
- **7.aQTL.calling/** – 等位基因特异 QTL 分析脚本：
  - `1.normalized.R` – 归一化伪汇总峰计数。
  - `2.hg38Tohg19.sh` – 将坐标从 hg38 转换到 hg19。
  - `3.attach.pos.sh` – 将基因组位置附加到归一化矩阵。
  - `4.bgzip_index.sh` – 排序后 bgzip 并索引矩阵文件。
  - `5.PEER.sh` – 提交 PEER 因子估计任务。
  - `6.generate.Fastqtl.running.script.sh` – 生成 FastQTL 任务脚本。
  - `attach.pos.pl` – `3.attach.pos.sh` 使用的 Perl 辅助脚本。
  - `get_Fastqtl.script.pl` – 创建 FastQTL 命令的 Perl 脚本。
  - `run_PEER.R` – PEER 因子估计工具。
- **8.run_differential_02c_NB_wRUV.non_vs_early.R** – 采用 Nebula 比较非 AD 与早期 AD 的差异峰。
- **8_2.run_differential_02c_NB_wRUV.all.R** – 在所有病理分组间的 Nebula 差异分析。
- **9.Cell.composition.R** – 利用 `speckle` 与 `limma` 评估细胞类型组成差异。
- **get_markerMat.PEC_markers_2.byCluster.re.pl** – 按簇汇总标记矩阵的 Perl 脚本。
</details>

## `Process.Morabito.etal.data`

处理 Morabito 等人的外部数据集。

- **1.setup.TSS1.R** – 从原始片段创建 ArchR arrows。
- **2.ArchR.clustering.R** – 进行迭代 LSI 与聚类。
- **3.assign.celltype.R** – 根据标记基因为簇指定细胞类型。
- **4_1.epigenome.chromHMM.R** – 以 ChromHMM 状态注释峰。
- **4_2.erosion.score.R** – 计算各状态的染色质侵蚀评分。

## `epiclust`

该文件夹在仓库中为空，但被定义为指向 [kellislab/epiclust](https://github.com/kellislab/epiclust) 的 git 子模块，内含峰-基因关联及模块算法的实现。克隆子模块即可获得代码。

---

本仓库同时包含 R、Python 与 shell 脚本。许多脚本假设特定的数据路径，运行前可能需要根据实际环境调整。请参阅各脚本内的注释以了解使用细节。
