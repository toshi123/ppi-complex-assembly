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
│   ├── download_clinvar.py
│   ├── download_gwas_catalog.py
│   ├── build_pdb_chain_db.py
│   ├── build_human_db.py
│   ├── filter_intact_human.py
│   ├── run_mmseqs_search.sh          # UniProt vs PDBチェーン検索
│   ├── convert_human_vs_pdb.sh       # mmseqs convertalis ラッパー
│   ├── check_mapping_mmseqs_sifts_intact.py # IntAct/mmseqs/SIFTS の QC
│   ├── map_intact_pairs_to_pdb.py    # IntActペア→PDBテンプレート一覧
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

# ClinVar VCF (GRCh38; optional, 変異解析用)
python scripts/download_clinvar.py

# GWAS Catalog associations (optional)
python scripts/download_gwas_catalog.py
```

ダウンロード後、`data/raw/` に例えば次のようなファイルができている想定です：

* `intact_20251106.txt`
* `uniprot_human_all_20251110.fasta`
* `pdb_seqres_20251111.txt`
* `pdb_chain_uniprot_20251111.csv.gz`
* `clinvar_grch38_YYYYMMDD.vcf.gz`（任意）
* `gwas_catalog_associations_YYYYMMDD.tsv.gz`（任意）

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

### 6. UniProt ↔ PDB チェーンの相同性検索（mmseqs2）

`run_mmseqs_search.sh` が `mmseqs search` を実行し、`human_db`（Query）と `pdb_chain_db`（Target）から `data/interim/mmseqs/human_vs_pdb` を生成します。

```bash
bash scripts/run_mmseqs_search.sh
```

スクリプトは `--min-seq-id 0.3` と `-e 1e-3` を既定値に固定しており、テンポラリディレクトリ `data/interim/mmseqs/tmp` を併用します。

```1:7:scripts/run_mmseqs_search.sh
mmseqs search \
  data/interim/mmseqs/human_db \
  data/interim/mmseqs/pdb_chain_db \
  data/interim/mmseqs/human_vs_pdb \
  data/interim/mmseqs/tmp \
  --min-seq-id 0.3 \
  -e 1e-3
```

再実行する場合は `human_vs_pdb*` を削除するか別プレフィックスを指定してください。

### 7. `convertalis` による m8 生成 & マッピング QC（任意）

`mmseqs search` の結果 DB から TSV (`human_vs_pdb.m8`) を作るには以下を実行します。

```bash
bash scripts/convert_human_vs_pdb.sh
```

内部では `mmseqs convertalis` を呼び出し、`query,target,pident,...,bits` までの 10 列を出力しています。

```47:54:scripts/convert_human_vs_pdb.sh
"${MMSEQS_BIN}" convertalis \
  "${HUMAN_DB}" \
  "${PDB_CHAIN_DB}" \
  "${RESULT_DB}" \
  "${RESULT_TSV}" \
  --format-output "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits"
```

IntAct / mmseqs / SIFTS の ID が一貫しているか手早く確認したいときは `scripts/check_mapping_mmseqs_sifts_intact.py` を使います。

```bash
python scripts/check_mapping_mmseqs_sifts_intact.py
```

この QC スクリプトは mmseqs の target を SIFTS チェーン (`chain_label`) と突き合わせ、IntAct に出現する UniProt との共通部分をサマリ表示します。

```270:298:scripts/check_mapping_mmseqs_sifts_intact.py
m8_with_sifts = m8_df.merge(
    sifts_chain_uniprot,
    left_on="target",
    right_on="chain_label",
    how="left",
    indicator=True,
)
...
connected = m8_with_sifts[mask_query_in_intact & mask_has_sifts]
print(f"[INFO] # hits with query∈IntAct & SIFTS UniProt : {n_connected_hits}")
```

### 8. IntAct ペアと PDB テンプレートの統合

`map_intact_pairs_to_pdb.py` は IntAct ペア、`human_vs_pdb.m8`、SIFTS を突き合わせ、`data/interim/intact_pairs_with_pdb_templates.tsv` を生成します。相同性の閾値（例：`--min-pident 0.30`、`--min-alnlen 30`）や 1 UniProt あたりのヒット上限（`--max-hits-per-uniprot`）を CLI 引数で調整できます。

```bash
python scripts/map_intact_pairs_to_pdb.py \
  --min-pident 0.30 \
  --min-alnlen 30 \
  --max-evalue 1e-5
```

スクリプト本体は IntAct ペア → SIFTS チェーン → PDB ID を順につなぎ、`uniprot_a/uniprot_b` ごとに同一 PDB に含まれる鎖ペアを列挙します。

```573:614:scripts/map_intact_pairs_to_pdb.py
pair_templates = build_intact_pair_templates(
    intact_df=intact_df,
    uniprot_chain_table=uniprot_chain_table,
    allow_self_chain=args.allow_self_chain,
)
pair_templates.to_csv(output_path, sep="\t", index=False)
print(f"[INFO] Written: {output_path}")
```

出力カラムは `pdb_id`, `chain_id_a/b`, `pident_a/b`, `evalue_a/b` などで構成され、後続の複合体組み上げステップのテンプレート一覧として利用します。

---

## 今後のパイプライン（進捗メモ）

1. **PPI グラフの構築**

   * `data/interim/intact_human_pairs.tsv` の作成までは完了済み。NetworkX でのグラフ構築・解析ユーティリティを今後追加予定。

2. **UniProt→PDB チェーン相同性検索（mmseqs2）**

   * `run_mmseqs_search.sh` と `convert_human_vs_pdb.sh` で `human_vs_pdb.m8` を自動生成できる状態。
   * `check_mapping_mmseqs_sifts_intact.py` により ID マッピングの QC も可能。

3. **PDB 複合体テンプレートリストの作成**

   * `map_intact_pairs_to_pdb.py` が `intact_pairs_with_pdb_templates.tsv` を出力するところまで整備済み。

4. **複合体組み上げ**

   * gemmi + TMalign/Mican による構造拡張ロジックを `build_complexes_from_graph.py` 以下で実装予定。

5. **インターフェース解析**

   * FreeSASA + gemmi で ΔSASA や接触残基を定義するスクリプト／モジュールを追加予定。

6. **疾患変異・物性解析**

   * `download_clinvar.py` / `download_gwas_catalog.py` でデータ取得まで対応済み。  
     変異マッピングと統計解析は `src/ppi_complex/interface` / `analysis` 以下に今後実装予定。

---

## Status

現段階での実装状況メモです。

* [x] データダウンロードスクリプトの実装

  * [x] IntAct (`download_intact.py`)
  * [x] UniProt human (`download_uniprot_human.py`)
  * [x] PDB SEQRES (`download_pdb_seqres.py`)
  * [x] SIFTS PDB–UniProt mapping (`download_sifts_mapping.py`)
  * [x] ClinVar VCF (`download_clinvar.py`)
  * [x] GWAS Catalog associations (`download_gwas_catalog.py`)
* [x] IntAct ヒト PPI ペア抽出

  * [x] `filter_intact_human.py` による human–human UniProt ペア作成
* [x] PDB チェーン FASTA & mmseqs2 DB

  * [x] `build_pdb_chain_db.py`
* [x] UniProt(Human) mmseqs2 DB

  * [x] `build_human_db.py`
* [x] UniProt ↔ PDB チェーン相同性検索（mmseqs2）

  * [x] `run_mmseqs_search.sh` ＋ `convert_human_vs_pdb.sh`
  * [x] `check_mapping_mmseqs_sifts_intact.py` による QC
* [x] PDB 複合体テンプレートリストの作成

  * [x] `map_intact_pairs_to_pdb.py`
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
