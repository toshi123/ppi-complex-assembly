# PPI Complex Assembly

ヒトのタンパク質相互作用ネットワーク上で、
既知複合体構造と相同性情報を利用して **巨大なタンパク質複合体を組み上げるプロジェクト** です。

- IntAct に登録された PPI をグラフとして扱い、
- PDB に存在する複合体（および相同な複合体）をテンプレートとして、
- 共通サブユニットの構造重ね合わせにより、より大きな複合体を構築します。

最終的には、構築した複合体のインターフェースに疾患関連変異（GWAS / ClinVar など）をマッピングし、
**インターフェースのタイプと病原性との関係**を解析することを目標としています。

このアイデアは、以下の論文で行った解析を拡張・再構築したものです：

- Tsuji et al., *Scientific Reports* (2015) – multi-interface hub / singlish-interface hub 解析 など

---

## Project Goals

- ヒト PPI（IntAct）から得られる **相互作用ネットワーク** を構築する
- PDB + UniProt + 相同性検索（mmseqs2 など）により、
  - **既知複合体テンプレート** を網羅的にリストアップする
- 共通サブユニットを介した構造重ね合わせにより、
  - **PPI グラフ上で複合体をどんどん拡張**する
- インターフェースを抽出し、
  - 多パートナー界面 / 単一パートナー界面 などに分類
  - 疾患変異のマッピング・物性比較・エンリッチメント解析を行う

---

## Repository Structure

このリポジトリは、だいたい次のような構成を想定しています（進行に応じて変える可能性あり）。

```text
.
├── README.md
├── LICENSE
├── pyproject.toml / requirements.txt / environment.yml
├── .gitignore
├── data/
│   ├── README.md         # 外部DBのダウンロード方法・バージョン情報
│   ├── raw/              # 生データ (IntAct, UniProt, PDBなど) → Git管理外
│   ├── interim/          # 中間生成物 (mmseqs結果など) → Git管理外
│   └── processed/        # 軽量な最終テーブルなど (必要に応じてGit管理)
├── src/
│   └── ppi_complex/
│       ├── __init__.py
│       ├── config.py
│       ├── download/     # IntAct, UniProt, PDB などのDL & 整形
│       ├── homology/     # mmseqs2, BLAST 等による相同性検索
│       ├── assembly/     # gemmi + TMalign/Mican による複合体組み上げ
│       ├── interface/    # インターフェース抽出, ΔSASA, 物性
│       ├── analysis/     # 変異・GWAS のマッピングと統計解析
│       └── utils/        # ログ, パス解決などの共通処理
├── scripts/
│   ├── download_intact.py
│   ├── download_uniprot_human.py
│   ├── download_pdb_seqres.py
│   ├── build_pdb_chain_db.py
│   ├── run_mmseqs_search.sh
│   ├── map_intact_pairs_to_pdb.py
│   └── build_complexes_from_graph.py
├── notebooks/
│   ├── 00_project_overview.ipynb
│   ├── 10_mmseqs_qc.ipynb
│   ├── 20_interface_properties.ipynb
│   └── 30_variant_enrichment.ipynb
├── results/
│   ├── figures/          # 図 (png, pdfなど)
│   ├── tables/           # 軽いTSV/CSV (論文用の最終テーブル)
│   └── logs/             # 実行ログなど
└── config/
    ├── paths.yaml        # データディレクトリ等 (環境依存はここに集約)
    └── project.yaml      # 閾値, バージョン情報など
