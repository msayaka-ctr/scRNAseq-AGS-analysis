# scRNA-seq Analysis of AGS (Aicardi-Goutières Syndrome)

自己炎症性疾患 **Aicardi-Goutières Syndrome (AGS)** に関するscRNAseqの解析再現。

---

## 概要

AGS患者（P1, P2）と健常コントロール（C1, C2）の末梢血単核細胞（PBMC）を対象に、scRNA-seqデータの統合・クラスタリング・差次発現解析（DEG）・I型インターフェロン応答スコアリングを実施。

**主な発見：**
- NK/MAIT/Tγδ細胞でAGS特異的なI型インターフェロン関連遺伝子（*MX1*, *OAS1*, *IFI44L*）の発現上昇
- AGS患者群ではほぼ全ての免疫細胞クラスターでI型インターフェロンスコアが有意に高値

---

##データ
-GEO:GSE220764
-論文:Batignes, M., Luka, M., Jagtap, S. et al. Pharmacological stabilization of hypoxia-inducible factor 1-α dampens the interferon response and promotes glycolysis in Aicardi-Goutières syndrome. Nat Commun 17, 3379 (2026). https://doi.org/10.1038/s41467-026-69979-9

## 解析フロー

```
統合済みSeuratオブジェクト読み込み
        ↓
サンプルサブセット（C1, C2, P1, P2）
        ↓
UMAP可視化・MacroCluster確認
        ↓
NK/MAIT/Tγδ細胞の抽出
        ↓
DEG解析（AGS vs Ctrl)
        ↓
I型インターフェロン遺伝子セットスコアリング（AddModuleScore）
        ↓
クラスター別t検定・有意差の可視化
```

---

## 主な解析内容

### 1. データの前処理・サブセット化

```r
# AGSサンプルとコントロールの抽出
seurat_obj2 <- subset(seurat_obj, SampleID_paper %in% c("C1", "C2", "P1", "P2_before_treatment"))

# ラベルの整理
levels(seurat_obj2$SampleID_paper)[levels(seurat_obj2$SampleID_paper) == "P2_before_treatment"] <- "P2"
```

### 2. UMAP可視化

- サンプル別・MacroCluster別のUMAP描画
- 細胞数の10,000細胞スケーリング比較

### 3. DEG解析（NK/MAIT/Tγδ細胞）

```r
DEG <- FindMarkers(NK_MAIT_Tgd,
                   ident.1 = c("P1", "P2"),
                   ident.2 = c("C1", "C2"),
                   group.by = "SampleID_paper")

# フィルタリング条件: p_val_adj < 0.05 & |avg_log2FC| >= 1
```

### 4. I型インターフェロンスコアリング

MSigDB の `GOBP_RESPONSE_TO_TYPE_I_INTERFERON` 遺伝子セットを使用し、`AddModuleScore` により各細胞にスコアを付与。クラスター別にCtrl vs AGSをt検定で比較し、有意差をVlnPlotに注釈。

---

## 環境・依存パッケージ

| パッケージ | 用途 |
|---|---|
| [Seurat](https://satijalab.org/seurat/) v5 | scRNA-seq 解析全般 |
| ggplot2 | Volcano Plot 等の可視化 |
| dplyr | データ操作 |

```r
library(Seurat)
library(ggplot2)
library(dplyr)
```

### 動作確認環境

- R >= 4.3
- Seurat >= 5.0

---

## ファイル構成

```
.
├── ADS_scRNAseq.Rmd                              # メイン解析ノートブック
├── GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2026.1.Hs.grp  # MSigDB 遺伝子セット
└── README.md
```

> **Note:** 元データ `Integration#1.rds` は元データ(GEO:GSE220764)よりダウンロード。

---

## 結果サマリー

| 解析 | 結果 |
|---|---|
| DEG (NK/MAIT/Tγδ, AGS vs Ctrl) | 有意な上昇遺伝子に *MX1*, *OAS1*, *IFI44L*, *CD52* など |
| IFNI スコア | AGS群で全MacroClusterにわたって有意に上昇（** p<0.01） |

---

## 作者

**Sayaka Matsushima**
