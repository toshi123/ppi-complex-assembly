#!/usr/bin/env python
"""
check_mapping_mmseqs_sifts_intact.py

目的:
- mmseqs 結果 (human_vs_pdb.m8)
- SIFTS PDB–UniProt マッピング (pdb_chain_uniprot_current.csv.gz)
- IntAct ヒト PPI ペア (intact_human_pairs.tsv)

がちゃんと「同じ ID 宇宙」でつながるか、ざっくり QC する。

出力:
- 標準出力に、以下のようなサマリを出す:
    - IntAct に出てくる UniProt ID の数
    - mmseqs の query 側 UniProt ID の数
    - その共通集合のサイズ
    - mmseqs の target チェーンのうち、SIFTS で UniProt が紐づく割合
    - 「IntAct のノードに属する UniProt × SIFTS で UniProt がとれる PDBチェーン」
      がどれくらい存在するか など
"""

import argparse
from pathlib import Path
import sys
import textwrap

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="QC for mapping between IntAct UniProt, mmseqs results, and SIFTS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            # デフォルトのファイルパスで実行
            python scripts/check_mapping_mmseqs_sifts_intact.py

            # パスを明示的に指定
            python scripts/check_mapping_mmseqs_sifts_intact.py \\
                --intact data/interim/intact_human_pairs.tsv \\
                --mmseqs data/interim/mmseqs/human_vs_pdb.m8 \\
                --sifts data/raw/pdb_chain_uniprot_current.csv.gz
            """
        ),
    )
    parser.add_argument(
        "--intact",
        type=str,
        default=None,
        help=(
            "IntAct ヒト PPI ペア TSV のパス "
            "(デフォルト: data/interim/intact_human_pairs.tsv)"
        ),
    )
    parser.add_argument(
        "--mmseqs",
        type=str,
        default=None,
        help=(
            "mmseqs 結果 human_vs_pdb.m8 のパス "
            "(デフォルト: data/interim/mmseqs/human_vs_pdb.m8)"
        ),
    )
    parser.add_argument(
        "--sifts",
        type=str,
        default=None,
        help=(
            "SIFTS mapping ファイル pdb_chain_uniprot_current.csv.gz のパス "
            "(デフォルト: data/raw/pdb_chain_uniprot_current.csv.gz)"
        ),
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    repo_root = Path(__file__).resolve().parents[1]

    intact_path = (
        repo_root / "data" / "interim" / "intact_human_pairs.tsv"
        if args.intact is None
        else (repo_root / args.intact if not Path(args.intact).is_absolute() else Path(args.intact))
    )

    mmseqs_path = (
        repo_root / "data" / "interim" / "mmseqs" / "human_vs_pdb.m8"
        if args.mmseqs is None
        else (repo_root / args.mmseqs if not Path(args.mmseqs).is_absolute() else Path(args.mmseqs))
    )

    sifts_path = (
        repo_root / "data" / "raw" / "pdb_chain_uniprot_current.csv.gz"
        if args.sifts is None
        else (repo_root / args.sifts if not Path(args.sifts).is_absolute() else Path(args.sifts))
    )

    return intact_path, mmseqs_path, sifts_path


def normalize_uniprot_from_query(q: str) -> str:
    """
    mmseqs の query ID から「素の UniProt アクセッション」を推定する。

    UniProt FASTA のヘッダが例えば:
        sp|Q9Y2X3|PROT_HUMAN ...
        tr|A0A123456|SOME_HUMAN ...

    のような場合、mmseqs の query は 'sp|Q9Y2X3|PROT_HUMAN' になるので、
    2番目のトークンを取る。

    もし '|' が含まれていなければ、そのまま返す。
    """
    if "|" in q:
        parts = q.split("|")
        if len(parts) >= 2:
            return parts[1]
    return q


def build_chain_label(pdb_id: str, chain_id: str) -> str:
    """
    SIFTS の pdb_id, chain_id から 'PDBID_CHAIN' 形式を作る。

    例:
        pdb_id = '1abc', chain_id = 'A' -> '1ABC_A'
    """
    if not isinstance(pdb_id, str) or not isinstance(chain_id, str):
        return None
    pdb_id = pdb_id.strip().upper()
    chain_id = chain_id.strip()
    if not pdb_id or not chain_id:
        return None
    return f"{pdb_id}_{chain_id}"


def main() -> None:
    args = parse_args()
    intact_path, mmseqs_path, sifts_path = resolve_paths(args)

    print(f"[INFO] IntAct pairs       : {intact_path}")
    print(f"[INFO] mmseqs results     : {mmseqs_path}")
    print(f"[INFO] SIFTS mapping      : {sifts_path}")

    # --- 1. IntAct ヒト PPI ペアの読み込み ---
    if not intact_path.exists():
        print(f"[ERROR] IntAct file not found: {intact_path}", file=sys.stderr)
        sys.exit(1)

    intact_df = pd.read_csv(intact_path, sep="\t")
    if not {"uniprot_a", "uniprot_b"}.issubset(intact_df.columns):
        print("[ERROR] IntAct TSV に 'uniprot_a', 'uniprot_b' 列が必要です。", file=sys.stderr)
        sys.exit(1)

    intact_uniprot_set = set(intact_df["uniprot_a"].astype(str)) | set(
        intact_df["uniprot_b"].astype(str)
    )

    print(f"[INFO] # IntAct pairs          : {len(intact_df)}")
    print(f"[INFO] # unique IntAct UniProt : {len(intact_uniprot_set)}")

    # --- 2. mmseqs human_vs_pdb.m8 の読み込み ---
    if not mmseqs_path.exists():
        print(f"[ERROR] mmseqs M8 file not found: {mmseqs_path}", file=sys.stderr)
        sys.exit(1)

    m8_cols = [
        "query",
        "target",
        "pident",
        "alnlen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bits",
    ]
    m8_df = pd.read_csv(mmseqs_path, sep="\t", header=None, names=m8_cols)

    print(f"[INFO] # mmseqs hits        : {len(m8_df)}")
    print(f"[INFO] # mmseqs queries     : {m8_df['query'].nunique()}")
    print(f"[INFO] # mmseqs targets     : {m8_df['target'].nunique()}")

    # query 側 UniProt ID の正規化
    m8_df["query_uniprot"] = m8_df["query"].astype(str).map(normalize_uniprot_from_query)

    n_unique_query_uniprot = m8_df["query_uniprot"].nunique()
    n_query_in_intact = len(set(m8_df["query_uniprot"]) & intact_uniprot_set)

    print(f"[INFO] # unique mmseqs query_uniprot : {n_unique_query_uniprot}")
    print(f"[INFO] # query_uniprot ∩ IntAct      : {n_query_in_intact}")

    # --- 3. SIFTS マッピングの読み込み ---
    if not sifts_path.exists():
        print(f"[ERROR] SIFTS file not found: {sifts_path}", file=sys.stderr)
        sys.exit(1)

    # 列名は SIFTS のフォーマットに依存するので、ここは必要に応じて調整してください。
    # 典型的な pdb_chain_uniprot.csv.gz では:
    #   pdb_id, chain_id, uniprot_acc, ...
    # コメント行 (# で始まる) をスキップして読み込む
    sifts_df = pd.read_csv(sifts_path, comment="#", low_memory=False)

    print(f"[INFO] SIFTS columns: {list(sifts_df.columns)}")

    # 典型的な SIFTS の列名は PDB, CHAIN, SP_PRIMARY などなので、
    # それらがあれば内部的に統一した名前にリネームする
    col_map = {}

    if {"pdb_id", "chain_id", "uniprot_acc"}.issubset(sifts_df.columns):
        # すでに想定どおりの列名の場合はそのまま使う
        col_map = {}
    elif {"PDB", "CHAIN", "SP_PRIMARY"}.issubset(sifts_df.columns):
        # よくある SIFTS 形式: PDB, CHAIN, SP_PRIMARY
        col_map = {
            "PDB": "pdb_id",
            "CHAIN": "chain_id",
            "SP_PRIMARY": "uniprot_acc",
        }
    else:
        print(
            "[ERROR] SIFTS CSV の列名が想定外です。\n"
            "  想定パターン:\n"
            "    - pdb_id, chain_id, uniprot_acc\n"
            "    - PDB, CHAIN, SP_PRIMARY\n"
            f"  実際の列: {list(sifts_df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    if col_map:
        sifts_df = sifts_df.rename(columns=col_map)

    expected_cols = {"pdb_id", "chain_id", "uniprot_acc"}
    if not expected_cols.issubset(sifts_df.columns):
        print(
            f"[ERROR] リネーム後も {expected_cols} が揃っていません。実際の列: {list(sifts_df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    sifts_df["chain_label"] = [
        build_chain_label(pdb_id, chain_id)
        for pdb_id, chain_id in zip(sifts_df["pdb_id"], sifts_df["chain_id"])
    ]
    sifts_df = sifts_df.dropna(subset=["chain_label"])

    # PDBチェーン→UniProt の対応表
    sifts_chain_uniprot = sifts_df[["chain_label", "uniprot_acc"]].drop_duplicates()
    print(f"[INFO] # SIFTS chain→UniProt rows : {len(sifts_chain_uniprot)}")
    print(f"[INFO] # unique SIFTS chain_label : {sifts_chain_uniprot['chain_label'].nunique()}")
    print(f"[INFO] # unique SIFTS UniProt     : {sifts_chain_uniprot['uniprot_acc'].nunique()}")

    sifts_df["chain_label"] = [
        build_chain_label(pdb_id, chain_id)
        for pdb_id, chain_id in zip(sifts_df["pdb_id"], sifts_df["chain_id"])
    ]
    sifts_df = sifts_df.dropna(subset=["chain_label"])

    # PDBチェーン→UniProt の対応表
    sifts_chain_uniprot = sifts_df[["chain_label", "uniprot_acc"]].drop_duplicates()
    print(f"[INFO] # SIFTS chain→UniProt rows : {len(sifts_chain_uniprot)}")
    print(f"[INFO] # unique SIFTS chain_label : {sifts_chain_uniprot['chain_label'].nunique()}")
    print(f"[INFO] # unique SIFTS UniProt     : {sifts_chain_uniprot['uniprot_acc'].nunique()}")

    # --- 4. mmseqs target と SIFTS chain_label のマージ ---
    m8_with_sifts = m8_df.merge(
        sifts_chain_uniprot,
        left_on="target",
        right_on="chain_label",
        how="left",
        indicator=True,
    )

    total_hits = len(m8_with_sifts)
    hits_with_sifts = (m8_with_sifts["_merge"] == "both").sum()

    print(f"[INFO] # mmseqs hits with SIFTS mapping : {hits_with_sifts} / {total_hits}")
    if total_hits > 0:
        print(f"[INFO]   fraction                      : {hits_with_sifts / total_hits:.3f}")

    # --- 5. 「IntAct のノード ↔ mmseqs ↔ SIFTS」がつながるかを見る ---
    # mmseqs query が IntAct に含まれていて、かつ SIFTS UniProt が存在するものを数える
    mask_query_in_intact = m8_with_sifts["query_uniprot"].isin(intact_uniprot_set)
    mask_has_sifts = m8_with_sifts["uniprot_acc"].notna()

    connected = m8_with_sifts[mask_query_in_intact & mask_has_sifts]

    n_connected_hits = len(connected)
    n_connected_queries = connected["query_uniprot"].nunique()
    n_connected_uniprot_sifts = connected["uniprot_acc"].nunique()

    print(f"[INFO] # hits with query∈IntAct & SIFTS UniProt : {n_connected_hits}")
    print(f"[INFO] # queries (IntAct) with SIFTS-mapped hits : {n_connected_queries}")
    print(f"[INFO] # SIFTS UniProt connected via mmseqs      : {n_connected_uniprot_sifts}")

    # 代表例を少しだけ表示
    print("\n[INFO] Example connected rows (head):")
    cols_show = [
        "query",
        "query_uniprot",
        "target",
        "uniprot_acc",
        "pident",
        "alnlen",
        "evalue",
    ]
    print(connected[cols_show].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
