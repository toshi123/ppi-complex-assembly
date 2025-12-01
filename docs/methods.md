# METHODS: Pipeline, Thresholds, and Provenance

各ステップを **目的 / 用いるスクリプト / 入力ファイル / 閾値** で明確化し、再現と変更追跡を容易にする。

> 環境: `environment.yml`（mmseqs2, gemmi, freesasa, pandas 等）

---

## Quick Start（最短導線）

```bash
# 0) 環境
mamba env create -f environment.yml
mamba activate ppi-complex-assembly
```
## 0. Raw Downloads

### 目的
必要な外部DBを data/raw/ に取得し、日付付きで保存する。

### 用いるスクリプト
- `scripts/download_intact.py`
- `scripts/download_uniprot_human.py`
- `scripts/download_pdb_seqres.py`
- `scripts/download_sifts_mapping.py`

### 入力ファイル
- 外部DL、ローカルに存在しない前提

### 閾値
- なし（取得のみ）。*_current.* のシンボリックリンクで最新版を参照可。

### 主な生成物
- `data/raw/intact_YYYYMMDD.txt`
- `data/raw/uniprot_human_all_YYYYMMDD.fasta`
- `data/raw/pdb_seqres_YYYYMMDD.txt`
- `data/raw/pdb_chain_uniprot_YYYYMMDD.csv.gz（/current.csv.gz）`

```bash
# 0) 原データの取得（ファイルは *_YYYYMMDD.* で保存）
python scripts/download_intact.py
python scripts/download_uniprot_human.py
python scripts/download_pdb_seqres.py
python scripts/download_sifts_mapping.py
```

## 1. IntAct → ヒト PPI ペア抽出

### 目的
- IntAct（PSIMITAB）からヒト同士のUniProtペアを抽出・正規化。

### 用いるスクリプト
- `scripts/filter_intact_human.py`

### 入力ファイル
- `data/raw/intact_YYYYMMDD.txt`

### 閾値
- 種フィルタ: Human × Human
- 無向化: (A,B) と (B,A) の重複排除

### 出力
- `data/interim/intact_human_pairs.tsv`

```bash
# 1) 前処理 → 検索 → マッピング
python scripts/filter_intact_human.py --in data/raw/intact_YYYYMMDD.txt \
  --out data/interim/intact_human_pairs.tsv
```

## 2. mmseqs2 DB 構築（PDBチェーン / ヒトUniProt）

### 目的
- 相同性検索クエリ/ターゲットDBの準備。

### 用いるスクリプト
- scripts/build_pdb_chain_db.py
- scripts/build_human_db.py

### 入力ファイル
- PDB チェーン: pdb_seqres_*.txt, pdb_chain_uniprot_*.csv.gz
- ヒト UniProt: uniprot_human_all_*.fasta

### 閾値
- なし（DB生成のみ）

### 出力
- `data/interim/mmseqs/pdb_chain_db*`
- `data/interim/mmseqs/human_db*`

```bash
# 2) mmseqs2 DB 構築
python scripts/build_pdb_chain_db.py \
  --seqres data/raw/pdb_seqres_YYYYMMDD.txt \
  --sifts  data/raw/pdb_chain_uniprot_YYYYMMDD.csv.gz \
  --out-prefix data/interim/mmseqs/pdb_chain_db

python scripts/build_human_db.py \
  --fasta data/raw/uniprot_human_all_YYYYMMDD.fasta \
  --out-prefix data/interim/mmseqs/human_db
```

## 3. UniProt ↔ PDB チェーン相同性検索（mmseqs2）

### 目的
- ヒトUniProt（query）とPDBチェーン（target）の対応候補を得る。

### 用いるスクリプト
- `scripts/run_mmseqs_search.sh（mmseqs search）`
  - （必要に応じて）`scripts/convert_human_vs_pdb.sh`（.m8 生成）

### 入力ファイル
- `data/interim/mmseqs/human_db*`, `data/interim/mmseqs/pdb_chain_db*`

### 閾値（既定例）
- `--min-seq-id 0.3`, `-e 1e-3`（スクリプト内既定）
  - `.m8` カラム: query,target,pident,alnlen,...,evalue,bits

### 出力
- `data/interim/mmseqs/human_vs_pdb.m8`

```bash
# 3) UniProt ↔ PDB チェーン相同性検索 
bash scripts/run_mmseqs_search.sh  # -> data/interim/mmseqs/human_vs_pdb.m8
```

## 4. SIFTS 統合 → PDB テンプレート候補

### 目的
- IntAct ペア × mmseqs ヒット × SIFTS（PDB chain ↔ UniProt）を統合し、同一PDB内で接触し得る鎖ペア候補を列挙。

### 用いるスクリプト
- `scripts/check_mapping_mmseqs_sifts_intact.py（QC）`
- `scripts/map_intact_pairs_to_pdb.py`

### 入力ファイル
- `data/interim/intact_human_pairs.tsv`
- `data/interim/mmseqs/human_vs_pdb.m8`
- `data/raw/pdb_chain_uniprot_current.csv.gz`

### 閾値（既定）
- ヒット採用: `pident ≥ 0.25`, `alnlen ≥ 20`, `evalue ≤ 1e-3`
- 自己同一鎖除去: `chain_id_a != chain_id_b`

### 出力
- `data/interim/intact_pairs_with_pdb_templates.tsv`

```bash
# 4) SIFTS 統合 → PDB テンプレート候補
python scripts/map_intact_pairs_to_pdb.py \
  --intact data/interim/intact_human_pairs.tsv \
  --mmseqs data/interim/mmseqs/human_vs_pdb.m8 \
  --sifts  data/raw/pdb_chain_uniprot_current.csv.gz \
  --out    data/interim/intact_pairs_with_pdb_templates.tsv
```

## 5. 物理接触フィルタ（座標参照・並列）

### 目的
- 構造ファイル（生体会合体優先）を取得し、heavy atom 距離にもとづく実接触のみ残す。

### 用いるスクリプト
- `scripts/filter_pdb_contacts_parallel.py`

### 入力ファイル
- `data/interim/intact_pairs_with_pdb_templates.tsv`
- PDB構造キャッシュ: `data/raw/pdb_structures/*.cif.gz`

### 閾値（既定）
- 近接距離: `distance-cutoff = 5.0 Å`
- 最小接触数: `min-contacts = 5`（原子ペア）
- アセンブリ優先: `--use-assembly`
- サイズ制限: `--max-file-size 50 (MB)`
- 並列: `--n-proc 6`（M1 Max 例）

### 出力
- `data/interim/intact_pairs_with_pdb_contacts.tsv`
    - 付帯列例: `min_atom_distance`, `n_atom_contacts`, `iface_res_* / iface_sasa_*`（有効時）

```bash
# 5) 座標参照の接触フィルタ → “tight”
python scripts/filter_pdb_contacts_parallel.py \
  --input data/interim/intact_pairs_with_pdb_templates.tsv \
  --output data/interim/intact_pairs_with_pdb_contacts.tsv \
  --pdb-dir data/raw/pdb_structures \
  --distance-cutoff 5.0 --min-contacts 5 --use-assembly --n-proc 6
```

## 6. “tight” フィルタ（過剰/異常界面の除去）

### 目的
- 過度なoverpackやアーチファクトを除去し、堅い候補集合へ。

### 用いるスクリプト
- `scripts/make_tight_filter.py`

### 入力ファイル
- `data/interim/intact_pairs_with_pdb_contacts.tsv`

### 閾値（既定）
- `max_atom_contacts ≤ 1500`
- `min_res_a ≥ 10`, `min_res_b ≥ 10`
- `1.6 ≤ min_atom_distance ≤ 4.5` (Å)
- `contacts_per_residue ≤ 40`

### 出力
- `data/interim/intact_pairs_with_pdb_contacts.tight.tsv`

```bash
# 6) “tight” フィルタ（過剰/異常界面の除去）
python scripts/make_tight_filter.py \
  --input  data/interim/intact_pairs_with_pdb_contacts.tsv \
  --output data/interim/intact_pairs_with_pdb_contacts.tight.tsv \
  --max-atom-contacts 1500 --min-res-a 10 --min-res-b 10 \
  --min-dist 1.6 --max-dist 4.5 --max-contacts-per-res 40
```

## 7. BSA（ΔSASA）付与 → band 抽出

### 目的
- FreeSASA で BSA（total ΔSASA）を算出し、目的のBSA帯で候補を絞る。

### 用いるスクリプト
- `scripts/add_bsa_with_freesasa.py`（sane/fixed 内包）

### 入力ファイル
- `data/interim/intact_pairs_with_pdb_contacts.tight.tsv`
- `data/raw/pdb_structures/*.cif.gz`

### 閾値（既定）
- BSA帯: `800–3200 Å²`
- `--assembly`, `--n-proc 6`

### 出力
- `data/interim/intact_pairs_with_pdb_contacts.tight_bsa.tsv`
- `data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv`（帯域抽出）

```bash
# 7) ΔSASA/BSA 付与 → band（800–3200 Å²）
python scripts/add_bsa_with_freesasa.py \
  --input  data/interim/intact_pairs_with_pdb_contacts.tight.tsv \
  --output data/interim/intact_pairs_with_pdb_contacts.tight_bsa.tsv \
  --pdb-dir data/raw/pdb_structures --assembly --n-proc 6 \
  --bsa-min 800 --bsa-max 3200

python - <<'PY'
import pandas as pd
i="data/interim/intact_pairs_with_pdb_contacts.tight_bsa.tsv"
o="data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv"
df=pd.read_csv(i,sep="\t")
df=df[df["bsa_total"].notna() & df["bsa_total"].between(800,3200)]
df.to_csv(o,sep="\t",index=False); print("[INFO] kept:",len(df),"->",o)
PY
```

## 8. Top-K テンプレート選抜（MMR: 質×多様性）

## 目的
- 距離・接触・BSAの合成スコアと界面多様性（Jaccard差・同一PDB/Asm/Ligペナルティ）により、各ペア K=3 の代表を選ぶ。

### 用いるスクリプト
- `scripts/select_template_topk_diverse.py`

### 入力ファイル
- `data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv`

### 閾値・重み（既定）
- `K: 3` / `λ: 0.5` / `--min-iface-diff 0.25`
- `--use-sasa-iface`（なければ距離由来にフォールバック）
- `--w-dist 1.0 --w-contacts 1.0 --w-bsa 0.5`
- `--w-iface 1.0 --w-lig 0.3 --w-same-pdb 0.2 --w-same-asm 0.2`

### 出力
- `data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv`

```bash
# 8) Top-K（MMR: 質×多様性）→ QC
python scripts/select_template_topk_diverse.py \
  --input  data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv \
  --output data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
  --k 3 --lambda 0.5 --min-iface-diff 0.25 --use-sasa-iface \
  --w-dist 1.0 --w-contacts 1.0 --w-bsa 0.5 \
  --w-iface 1.0 --w-lig 0.3 --w-same-pdb 0.2 --w-same-asm 0.2
```

## 9. QC（代表の健全性・多様性・カバレッジ）

### 目的
- Top-K 出力の妥当性を自動チェック（不足列は上流から補完）。

### 用いるスクリプト
- `scripts/qc_representatives.py`

### 入力ファイル
- `data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv`
- `--upstream` に `band / bsa / tight` などを列挙

### 閾値（既定）
- suspicious: `min_atom_distance > 4.5 Å` or `n_atom_contacts < 20`

### 出力
- `data/interim/qc_topk.summary.json`
- `data/interim/qc_topk.pairs_counts.tsv`
- `data/interim/qc_topk.suspicious.tsv`

```bash
# 9)  QC（代表の健全性・多様性・カバレッジ）
python scripts/qc_representatives.py \
  --topk data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
  --k 3 --max-dist 4.5 --min-contacts 20 \
  --upstream data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv \
            data/interim/intact_pairs_with_pdb_contacts.tight_bsa.tsv \
            data/interim/intact_pairs_with_pdb_contacts.tight.tsv \
  --out-prefix data/interim/qc_topk
```




