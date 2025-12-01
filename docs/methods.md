# METHODS: Pipeline, Thresholds, and Provenance

各ステップを **目的 / 用いるスクリプト / 入力ファイル / 閾値** で明確化し、再現と変更追跡を容易にする。

> 環境: `environment.yml`（mmseqs2, gemmi, freesasa, pandas 等）

---

## Quick Start（最短導線）

```bash
# 0) 環境
conda env create -f environment.yml
conda activate ppi-complex-assembly
```
## 0. Raw Downloads

### 目的
必要な外部DBを data/raw/ に取得し、日付付きで保存する。
以下の外部データベースから必要な情報を取得

| ファイル | 内容 |
|----------|------|
| IntAct | 実験的に検証されたPPIのカタログ |
| UniProt FASTA | ヒトタンパク質の配列 |
| PDB seqres | PDBエントリに含まれる配列 |
| SIFTS | PDBチェーン ↔ UniProtの対応表（EBIが提供） |

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
---

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

---

## 2. mmseqs2 DB 構築（PDBチェーン / ヒトUniProt）

### 目的
- 相同性検索のために2つのデータベースを作成：

1. **human_db**: ヒトUniProtタンパク質の配列DB（クエリ側）
2. **pdb_chain_db**: PDBチェーンの配列DB（ターゲット側）

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

---

## 3. UniProt ↔ PDB チェーン相同性検索（mmseqs2）

### 目的
- ヒトUniProt（query）とPDBチェーン（target）の対応候補を得る。
- これにより「このヒトタンパク質の構造はPDBのどこに存在するか」を見つける
- IntActにはPPIの情報（タンパク質Aとタンパク質Bが相互作用する）はあるが、**その複合体の立体構造がどこにあるか**は記載されていない。
- この検索により、各UniProtタンパク質に対応するPDBチェーンを見つける。

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

---

## 4. SIFTS 統合 → PDB テンプレート候補

### 目的
- IntAct ペア × mmseqs ヒット × SIFTS（PDB chain ↔ UniProt）を統合し、同一PDB内で接触し得る鎖ペア候補を列挙。

### 用いるスクリプト
- `scripts/check_mapping_mmseqs_sifts_intact.py（QC）`
- `scripts/map_intact_pairs_to_pdb.py`

### 処理の流れ (`map_intact_pairs_to_pdb.py`):
1. IntActペア `(UniProt_A, UniProt_B)` を取得
2. mmseqsから `UniProt_A → PDB_chain_X` のマッピングを取得
3. 同様に `UniProt_B → PDB_chain_Y` のマッピングを取得
4. **同じPDB ID内に `chain_X` と `chain_Y` が両方存在する**場合のみを「複合体テンプレート候補」として採用

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
- Step 4で得られた候補は「同じPDBに両チェーンが存在する」だけで、**実際に物理的に接触しているかは不明**。
- PDB構造を取得し、原子座標から**実際に接触している**候補のみを残す。

### 用いるスクリプト
- `scripts/filter_pdb_contacts_parallel.py`

### 処理 (`filter_pdb_contacts_parallel.py`):

1. PDB構造（BinaryCIF優先）をダウンロード/キャッシュ
2. **生物学的会合体（Biological Assembly）** を展開（`--use-assembly`）
   - 結晶学的な非対称単位ではなく、生理的に意味のある複合体状態を使用
3. チェーンA/B間の重原子（non-hydrogen）距離を計算
4. **閾値フィルタ**:
   - 最短原子間距離 ≤ 5.0 Å
   - 接触原子ペア数 ≥ 5

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

---

## 6. “tight” フィルタ（過剰/異常界面の除去）

### 目的
- Step 5の出力から、**結晶学的アーティファクトや異常な界面**を除去して信頼性の高い候補を得る。
- `min_atom_distance < 1.6 Å`: ファンデルワールス半径以下→構造エラー
- 過大な接触数: 非生理的な界面や複数界面の混同
- 点接触: 偶然の近接であり機能的界面ではない可能性

### 用いるスクリプト
- `scripts/make_tight_filter.py`

### 入力ファイル
- `data/interim/intact_pairs_with_pdb_contacts.tsv`

### 閾値（既定）
| 条件 | 閾値 | 意図 |
|------|------|------|
| `n_atom_contacts ≤ 1500` | 上限 | 過度に巨大な界面を排除 |
| `n_res_contacts_a ≥ 10` | 下限 | 点接触（たまたま近い）を排除 |
| `n_res_contacts_b ≥ 10` | 下限 | 同上 |
| `1.6 ≤ min_atom_distance ≤ 4.5 Å` | 範囲 | クラッシュ（原子重複）や遠すぎを排除 |
| `contacts_per_residue ≤ 40` | 上限 | 過密界面（モデリングエラー等）を排除 |

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
---

## 7. BSA（ΔSASA）付与 → band 抽出

### 目的
- FreeSASA で BSA（total ΔSASA）を算出し、目的のBSA帯で候補を絞る。
- BSAとは
  - 複合体形成により「埋没する」表面積
  - 界面サイズの定量的指標

### 用いるスクリプト
- `scripts/add_bsa_with_freesasa.py`（sane/fixed 内包）

### 処理
1. FreeSASAライブラリでSASAを計算
2. 単独チェーンA、単独チェーンB、複合体ABそれぞれのSASAを算出
3. BSA = SASA(A) + SASA(B) - SASA(AB) を計算
4. **ΔSASA由来の界面残基集合** (`iface_sasa_*`) を同定
   - 複合体形成でSASAが減少した残基 = 界面残基

### 入力ファイル
- `data/interim/intact_pairs_with_pdb_contacts.tight.tsv`
- `data/raw/pdb_structures/*.cif.gz`

### 閾値（既定）
- BSA帯: `800–3200 Å²`
    - なぜこの範囲か
      - **800 Å²未満**: 小さすぎる界面（一過性・弱い相互作用）
      - **3200 Å²超**: 大きすぎる界面（結晶パッキング接触や、複数界面の混合の可能性）
      - この範囲は生物学的に意味のあるタンパク質-タンパク質界面のサイズとして一般的に受け入れられている
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

---

## 8. Top-K テンプレート選抜（MMR: 質×多様性）

## 目的
- 各PPIペアに対して、**質（スコア）と多様性（異なる構造コンフォメーション）**のバランスが取れたK個（デフォルト3）の代表テンプレートを選出する。


### 用いるスクリプト
- `scripts/select_template_topk_diverse.py`

### スコア計算 (`select_template_topk_diverse.py`):

1. ベーススコア（Quality）:
   - `dist_score`: 小さい距離ほど高スコア（1.6〜4.5 Åを線形変換）
   - `contacts_score`: 接触数のlog1p正規化
   - `bsa_score`: BSAのlog1p正規化

2. 多様性ペナルティ（Similarity）:
   - 界面Jaccard類似度: 界面残基集合の重なり（`iface_sasa_*` 優先）
   - 同一PDBペナルティ: 同じPDBからの選出にペナルティ
   - 同一Assemblyペナルティ: 同じ生物学的会合体からの選出にペナルティ

#### 選抜プロセス:
1. 最初は最高ベーススコアのテンプレートを選択
2. 以降は `λ × Quality − (1−λ) × (選抜済みとの最大類似度)` を最大化するものを選択
3. Jaccard距離 < 0.25 の類似しすぎるテンプレートは除外

#### なぜMMRか:
- 単純に「最良スコア上位K個」だと、同じPDBの同じような界面ばかりになる
- 構造予測等の下流タスクでは**多様な参照構造**があると精度が上がる
- 同一タンパク質が異なるコンフォメーションで結晶化されている場合、それらを複数選びたい

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
---

## 9. QC（代表の健全性・多様性・カバレッジ）

### 目的
- Top-K 出力の妥当性を自動チェック（不足列は上流から補完）。

### 用いるスクリプト
- `scripts/qc_representatives.py`

### 検査項目 (`qc_representatives.py`):
1. カバレッジ: 各ペアで何件の代表が選ばれたか（K未満のペアを検出）
2. 多様性: 
   - `iface_source` の内訳（SASA由来 vs 距離由来）
   - ペアあたりの異なるPDB比率
3. 健全性: 「疑わしい代表」の検出
   - `min_atom_distance > 4.5 Å` (遠すぎる)
   - `n_atom_contacts < 20` (接触が少なすぎる)

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




