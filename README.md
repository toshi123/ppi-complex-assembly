# PPI Complex Assembly

ヒトのタンパク質相互作用ネットワーク上で、  
既知複合体構造と相同性情報を利用して **大きなタンパク質複合体を組み上げる** プロジェクトです。

---

## プロジェクト概要

- ヒト PPI（IntAct）をグラフとして扱い、PDB の複合体構造および相同性情報を用いて **巨大な複合体モデル** を構築する。
- IntAct の PPI をノード（タンパク質）・エッジ（相互作用）とするグラフにする。
- PDB に存在する複合体（＋30%以上の相同性を持つ相同複合体）を **テンプレート** として利用する。
- 共通サブユニット（例：タンパク質 A）で複合体を重ね合わせ、ABC… と拡張していく。

最終目標は：

- 構築した複合体のインターフェースに GWAS / ClinVar などの疾患変異をマッピングし、
- **インターフェースのタイプ（多パートナー vs 単一パートナー）と病原性・物性との関係** を解析することです。

---

## このプロジェクトの背景

このアイデアは、以下の論文で行った解析を拡張・再構築したものです：

- Tsuji et al., *Scientific Reports* (2015) – multi-interface hub / single-interface hub 解析 など

---

## リポジトリ構成（想定）

ルート例：`ppi-complex-assembly/`

```text
ppi-complex-assembly/
├── README.md
├── LICENSE
├── environment.yml          # conda 環境定義
├── .gitignore
├── data/
│   ├── README.md            # DBのダウンロード方法やバージョン情報
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
│   ├── download_sifts_mapping.py
│   ├── build_pdb_chain_db.py
│   ├── build_human_db.py
│   ├── filter_intact_human.py
│   ├── run_mmseqs_search.sh          # UniProt vs PDBチェーンの検索（予定）
│   ├── map_intact_pairs_to_pdb.py    # IntActペア→PDBテンプレート（予定）
│   └── build_complexes_from_graph.py # 複合体組み上げ（予定）
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
    ├── paths.yaml           # データディレクトリなど（環境依存情報）
    └── project.yaml         # 閾値, 外部DBバージョンなど
````

---

## Project Goals

* ヒト PPI（IntAct）から得られる **相互作用ネットワーク** を構築する。
* PDB + UniProt + 相同性検索（mmseqs2 など）により、

  * **既知複合体テンプレート** を網羅的にリストアップする。
* 共通サブユニットを介した構造重ね合わせにより、

  * **PPI グラフ上で複合体をどんどん拡張** する。
* インターフェースを抽出し、

  * 多パートナー界面 / 単一パートナー界面 などに分類。
  * 疾患変異のマッピング・物性比較・エンリッチメント解析を行う。

---

## インストール（環境構築）

conda を前提に環境を構築します。

```bash
git clone https://github.com/USER/ppi-complex-assembly.git
cd ppi-complex-assembly

# conda 環境の作成
conda env create -f environment.yml
conda activate ppi-complex-assembly
```

`environment.yml` の想定（概要）：

* Python 3.11
* numpy, pandas, scipy, matplotlib, scikit-learn
* gemmi, biopython, freesasa, networkx
* jupyterlab, ipykernel
* mmseqs2 など

---

## Quickstart: データ取得と前処理（現段階）

現時点で実装済みのステップは：

1. 外部DBからの生データダウンロード
2. IntAct からヒト PPI ペア（UniProt ペア）の抽出
3. PDB チェーン FASTA の作成
4. UniProt / PDB チェーンの mmseqs2 DB 作成

です。ここまでが「土台作り」です。

### 0. 前提

```bash
conda activate ppi-complex-assembly
```

`data/raw/` に書き込み可能であること。

### 1. 生データのダウンロード

```bash
# IntAct (全 PSIMITAB)
python scripts/download_intact.py

# UniProt (Human; all entries)
python scripts/download_uniprot_human.py --reviewed false

# PDB SEQRES
python scripts/download_pdb_seqres.py

# SIFTS PDB–UniProt mapping
python scripts/download_sifts_mapping.py
```

ダウンロード後、`data/raw/` に例えば次のようなファイルができている想定です：

* `intact_20251106.txt`
* `uniprot_human_all_20251110.fasta`
* `pdb_seqres_20251111.txt`
* `pdb_chain_uniprot_20251111.csv.gz`

### 2. 現行版ファイルへのシンボリックリンク（任意だが推奨）

解析コード側からは固定名を参照できるように、`*_current.*` というリンクを張ります：

```bash
cd data/raw

ln -s intact_20251106.txt intact_current.txt
ln -s uniprot_human_all_20251110.fasta uniprot_human_current.fasta
ln -s pdb_seqres_20251111.txt pdb_seqres_current.txt
ln -s pdb_chain_uniprot_20251111.csv.gz pdb_chain_uniprot_current.csv.gz
```

以降のスクリプトは、これらの `*_current` を前提に動きます。

### 3. IntAct からヒト PPI ペアを抽出

IntAct PSIMITAB から「ヒト同士の UniProt ペア」を取り出します。

```bash
python scripts/filter_intact_human.py
```

* 入力: `data/raw/intact_current.txt`
* 出力: `data/interim/intact_human_pairs.tsv`

出力ファイルの中身のイメージ：

```text
uniprot_a    uniprot_b
P12345       Q67890
...
```

ペアは（デフォルトでは）無向グラフとして扱うため、`(A,B)` と `(B,A)` は統合されます。

### 4. PDB チェーン FASTA と mmseqs2 DB の作成

```bash
# pdb_seqres_current.txt から PDBチェーンごとの FASTA を作成
# （同時に mmseqs2 DB も作成）
python scripts/build_pdb_chain_db.py --build-mmseqs
```

* 入力:

  * `data/raw/pdb_seqres_current.txt`
* 出力:

  * `data/interim/pdb_chain_current.fasta`
  * `data/interim/mmseqs/pdb_chain_db*` （mmseqs2 DB 一式）

### 5. ヒト UniProt の mmseqs2 DB の作成

```bash
python scripts/build_human_db.py
```

* 入力:

  * `data/raw/uniprot_human_current.fasta`
* 出力:

  * `data/interim/mmseqs/human_db*` （mmseqs2 DB 一式）

ここまでで：

* Query DB: `human_db`
* Target DB: `pdb_chain_db`

が揃った状態になります。

> 次のステップ（予定）
> `scripts/run_mmseqs_search.sh` を用意し、
> `mmseqs search human_db pdb_chain_db ...` による
> UniProt→PDB チェーンの相同性テーブル (`human_vs_pdb.m8`) を生成する。

---

## 今後のパイプライン（予定）

ここから先はまだ実装途中ですが、目標としているフローは以下の通りです。

1. **PPI グラフの構築**

   * `data/interim/intact_human_pairs.tsv` から NetworkX などでグラフ構築。

2. **UniProt→PDB チェーン相同性検索（mmseqs2）**

   * `run_mmseqs_search.sh` で `human_db` vs `pdb_chain_db` を検索。
   * 結果を `human_vs_pdb.m8` として出力。

3. **PDB 複合体テンプレートリストの作成**

   * `pdb_chain_uniprot_current.csv.gz`（SIFTS）と mmseqs 結果を組み合わせ、
   * IntActペア (A,B) に対応する複合体テンプレート (PDBID, chainA, chainB, pident, coverage, ...) のテーブルを作る。

4. **複合体組み上げ**

   * テンプレート付きエッジだけからなるサブグラフを抽出。
   * 共通サブユニットを介して gemmi + TMalign/Mican で複合体を拡張。
   * ステリッククラッシュなどをチェックし、不合理なモデルを除外。

5. **インターフェース解析**

   * FreeSASA + gemmi で ΔSASA を計算。
   * 接触残基・インターフェースパッチの定義。
   * 多パートナー界面 / 単一パートナー界面などに分類。

6. **疾患変異・物性解析**

   * ClinVar / GWAS の変異をインターフェース残基にマッピング。
   * 疾患変異のエンリッチメント解析。
   * インターフェースの物性（疎水性、電荷、パッチサイズなど）との関係を解析。

---

## Status

現段階での実装状況メモです。

* [x] データダウンロードスクリプトの実装

  * [x] IntAct (`download_intact.py`)
  * [x] UniProt human (`download_uniprot_human.py`)
  * [x] PDB SEQRES (`download_pdb_seqres.py`)
  * [x] SIFTS PDB–UniProt mapping (`download_sifts_mapping.py`)
* [x] IntAct ヒト PPI ペア抽出

  * [x] `filter_intact_human.py` による human–human UniProt ペア作成
* [x] PDB チェーン FASTA & mmseqs2 DB

  * [x] `build_pdb_chain_db.py`
* [x] UniProt(Human) mmseqs2 DB

  * [x] `build_human_db.py`
* [ ] UniProt ↔ PDB チェーン相同性検索（mmseqs2）

  * [ ] `run_mmseqs_search.sh` による `human_vs_pdb.m8` の生成
* [ ] PDB 複合体テンプレートリストの作成
* [ ] 複合体組み上げコード（gemmi + TMalign/Mican）
* [ ] インターフェース抽出 & ΔSASA 計算
* [ ] 疾患変異のマッピングと解析

---

## ライセンス / データ利用上の注意

* このリポジトリはコード・設定・ドキュメントのみを含み、
  IntAct / UniProt / PDB / SIFTS / ClinVar / GWAS などの **生データは再配布しません**。
* 各データベースの利用規約・ライセンスに従って、各自の環境でダウンロードしてください。
* ソースコード自体のライセンス（MIT など）は `LICENSE` を参照してください。

```
::contentReference[oaicite:0]{index=0}
```
