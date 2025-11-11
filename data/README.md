# Data Directory

このディレクトリは、すべての外部データベースと、その中間生成物・処理済みデータを置く場所です。  
**生のデータファイル（外部DB由来のもの）は GitHub にコミットしません。**

ディレクトリ構成の想定は以下の通りです：

```text
data/
├── README.md       # このファイル
├── raw/            # ダウンロードしてきた生データ (Git管理外)
├── interim/        # 中間生成物 (Git管理外)
└── processed/      # 軽い最終テーブルなど (必要に応じてGit管理)

## IntAct
- File: `data/raw/intact_20251106.txt`
- Downloaded: 2025-11-06
- Notes: Full IntAct psimitab

## UniProt (Homo sapiens; all)
- File: `data/raw/uniprot_human_all_20251110.fasta`
- Downloaded: 2025-11-10
- Query: organism_id:9606 (reviewed:false)

## PDB SEQRES
- File: `data/raw/pdb_seqres_20251111.txt`
- Downloaded: 2025-11-11

## SIFTS PDB–UniProt mapping
- File: `data/raw/pdb_chain_uniprot_20251111.csv.gz`
- Downloaded: 2025-11-11