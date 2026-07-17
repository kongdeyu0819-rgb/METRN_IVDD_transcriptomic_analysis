# Supplementary files for BMC Musculoskeletal Disorders manuscript
# METRN–FTX–PIEZO1 transcriptomic analysis in IVDD (v11.2)

本文件夹包含稿件 Availability 段引用、需随 GitHub 仓库一并发布（并触发 Zenodo 新版本归档）的补充材料。
所有文件均由本会话（2026-07-16/17）生成，对应稿件 §3.3、Supplementary Material 段及 Table S1 / S2。

## 文件清单与稿件对应关系

| 文件名 | 稿件引用 | 说明 |
|--------|----------|------|
| `Supplementary_Fig_S1.png` | Supplementary Fig. S1 | 单细胞 METRN/FTX/PIEZO1 三联图（A 各细胞型阳性率；B/C METRN–FTX、METRN–PIEZO1 散点），300 DPI |
| `Supplementary_Table_S1.csv` | Supplementary Table S1 | GSE199866 各细胞型 METRN/FTX/PIEZO1 阳性率与均值（= helix_inputs/Supp_scRNA_METRN_FTX_PIEZO1.csv） |
| `Supplementary_Table_S2_limma_GSE124272.csv` | Supplementary S2 | GSE124272 limma 全量差异表达结果（FDR） |
| `Supplementary_Table_S2_limma_GSE150408.csv` | Supplementary S2 | GSE150408 limma 全量差异表达结果（FDR） |
| `Supplementary_Table_S2_limma_GSE23130.csv` | Supplementary S2 | GSE23130 limma 全量差异表达结果（FDR） |
| `Supp_scRNA_X3_raw.npy` | （可复现原始数据） | GSE199866 全 14,001 细胞 × 三基因表达矩阵（numpy） |
| `Supp_scRNA_cells.json` | （可复现原始数据） | 上述细胞对应的细胞型注释标签 |

## 推送步骤（在仓库本地克隆中执行）

1. 将本文件夹内全部文件复制到仓库克隆的 `supplementary/` 目录（如仓库无此目录则新建）。
2. 提交并推送：
   ```
   git add supplementary/
   git commit -m "Add Supplementary Fig. S1, Table S1 and S2 limma results for v11.2 manuscript"
   git push origin main
   ```
3. 在 GitHub 仓库页面 **Draft a new release**，版本号 `v2.0.3`（或下一个递增版本），
   发布后 Zenodo GitHub 集成会自动归档新版本；概念 DOI `10.5281/zenodo.21261995` 将指向含本补充文件的最新版本。

> 注意：概念 DOI 永远指向最新版本；打新 release 后稿件中的 DOI 链接无需再改。
> 若仅推到 main 而不打 release，GitHub 链接可见但 Zenodo 归档（v2.0.2）仍不含这些文件，Availability 仍不完整——务必打 release。
