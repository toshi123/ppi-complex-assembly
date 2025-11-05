# PPI Complex Assembly プロジェクト共有メモ

※ このドキュメントは、Cursor の AI へのコンテキスト共有用のまとめです。  
リポジトリ構成、環境、スクリプトのたたき台をここに集約しています。

---

## プロジェクト概要

- ヒト PPI（IntAct）をグラフとして扱い、PDB の複合体構造および相同性情報を用いて **大きなタンパク質複合体を組み上げる** プロジェクト。

- 基本アイデア：
  - IntAct の PPI をノード（タンパク質）・エッジ（相互作用）とするグラフにする
  - PDB に存在する複合体（＋30%以上の相同性を持つ相同複合体）をテンプレートとして利用
  - 共通サブユニット（例：タンパク質 A）で複合体を重ね合わせ、ABC…と拡張していく

- 最終目標：
  - 構築した複合体のインターフェースに GWAS / ClinVar などの疾患変異をマッピング
  - **インターフェースのタイプ（多パートナー vs 単一パートナー）と病原性・物性との関係** を解析

---

## リポジトリの想定構成

ルート例：`ppi-complex-assembly/`

```text
ppi-complex-assembly/
├── README.md
├── LICENSE
├── environment.yml          # conda 環境定義
├── .gitignore
├── data/
│   ├── README.md            # DBのダウンロード方法やバージョン管理
│   ├── raw/                 # 生データ (IntAct, UniProt, PDBなど) → Git管理外
│   ├── interim/             # 中間生成物 → Git管理外
│   └── processed/           # 軽量な最終テーブルなど (必要に応じてGit管理)
├── src/
│   └── ppi_complex/
│       ├── __init__.py
│       ├── config.py        # 設定（必要に応じて）
│       ├── download/        # IntAct, UniProt, PDB などのDL & 整形
│       ├── homology/        # mmseqs2, BLAST などによる相同性検索
│       ├── assembly/        # gemmi + TMalign/Mican による複合体組み上げ
│       ├── interface/       # インターフェース抽出, ΔSASA, 物性解析
│       ├── analysis/        # 変異・GWAS のマッピングと統計解析
│       └── utils/           # ログ, パス解決など共通処理
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
│   ├── figures/             # 図 (png, pdf)
│   ├── tables/              # 軽いTSV/CSV (最終結果等)
│   └── logs/                # 実行ログなど
└── config/
    ├── paths.yaml           # データディレクトリなど（環境依存情報はここに）
    └── project.yaml         # 閾値, 外部DBバージョンなど
```

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

## このプロジェクトの背景

このアイデアは、以下の論文で行った解析を拡張・再構築したものです：

- Tsuji et al., *Scientific Reports* (2015) – multi-interface hub / singlish-interface hub 解析 など
