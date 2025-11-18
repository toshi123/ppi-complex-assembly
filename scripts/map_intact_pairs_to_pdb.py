#!/usr/bin/env python
"""
map_intact_pairs_to_pdb.py

目的:
- IntAct ヒト PPI ペア (uniprot_a, uniprot_b)
- mmseqs2 の UniProt(human) vs PDBチェーン結果 (human_vs_pdb.m8)
- SIFTS PDB–UniProt 対応 (pdb_chain_uniprot_current.csv.gz)

を統合して、

  IntAct ペアごとに
    - A に対応する PDB チェーン
    - B に対応する PDB チェーン
    - かつ同じ PDB ID 内に存在するもの

を「複合体テンプレート候補」として列挙する。

出力:
- デフォルト: data/interim/intact_pairs_with_pdb_templates.tsv

列のイメージ:
  uniprot_a, uniprot_b,
  pdb_id,
  chain_id_a, chain_label_a, sifts_uniprot_a,
  pident_a, alnlen_a, evalue_a, bits_a,
  chain_id_b, chain_label_b, sifts_uniprot_b,
  pident_b, alnlen_b, evalue_b, bits_b
"""

import argparse
from pathlib import Path
import sys
import textwrap

import numpy as np
import pandas as pd


# ---------- ユーティリティ ----------

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
    if not isinstance(q, str):
        q = str(q)
    if "|" in q:
        parts = q.split("|")
        if len(parts) >= 2:
            return parts[1]
    return q


def build_chain_label(pdb_id: str, chain_id: str) -> str | None:
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


# ---------- 引数関連 ----------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map IntAct human PPI pairs to PDB complex templates via mmseqs2 + SIFTS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            # デフォルトパスで実行
            python scripts/map_intact_pairs_to_pdb.py

            # 閾値を少し厳しめにする例
            python scripts/map_intact_pairs_to_pdb.py \\
                --min-pident 0.30 \\
                --min-alnlen 30 \\
                --max-evalue 1e-5 \\
                --max-hits-per-uniprot 50

            # 出力パスを変える
            python scripts/map_intact_pairs_to_pdb.py \\
                --output data/interim/intact_pairs_with_pdb_templates_v2.tsv
            """
        ),
    )
    parser.add_argument(
        "--intact",
        type=str,
        default=None,
        help="IntAct ヒト PPI ペア TSV のパス (デフォルト: data/interim/intact_human_pairs.tsv)",
    )
    parser.add_argument(
        "--mmseqs",
        type=str,
        default=None,
        help="mmseqs 結果 human_vs_pdb.m8 のパス (デフォルト: data/interim/mmseqs/human_vs_pdb.m8)",
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
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help=(
            "IntAct ペアごとの PDB テンプレートリストの出力パス。"
            "デフォルト: data/interim/intact_pairs_with_pdb_templates.tsv"
        ),
    )
    parser.add_argument(
        "--min-pident",
        type=float,
        default=0.25,
        help="mmseqs ヒットの最小 sequence identity (0〜1). デフォルト: 0.25",
    )
    parser.add_argument(
        "--min-alnlen",
        type=int,
        default=20,
        help="mmseqs ヒットの最小アラインメント長。デフォルト: 20",
    )
    parser.add_argument(
        "--max-evalue",
        type=float,
        default=1e-3,
        help="mmseqs ヒットの最大 e-value。デフォルト: 1e-3",
    )
    parser.add_argument(
        "--max-hits-per-uniprot",
        type=int,
        default=50,
        help="1 UniProt あたりの PDB チェーン候補数の上限。デフォルト: 50",
    )
    parser.add_argument(
        "--allow-self-chain",
        action="store_true",
        help=(
            "uniprot_a == uniprot_b の自己相互作用で、"
            "同じ chain_id を A/B 両方に割り当てるケースも残したい場合に指定。"
            "デフォルトでは、uniprot_a == uniprot_b かつ chain_id_a == chain_id_b の行は除外される。"
        ),
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    repo_root = Path(__file__).resolve().parents[1]

    if args.intact is None:
        intact_path = repo_root / "data" / "interim" / "intact_human_pairs.tsv"
    else:
        intact_path = Path(args.intact)
        if not intact_path.is_absolute():
            intact_path = (repo_root / intact_path).resolve()

    if args.mmseqs is None:
        mmseqs_path = repo_root / "data" / "interim" / "mmseqs" / "human_vs_pdb.m8"
    else:
        mmseqs_path = Path(args.mmseqs)
        if not mmseqs_path.is_absolute():
            mmseqs_path = (repo_root / mmseqs_path).resolve()

    if args.sifts is None:
        sifts_path = repo_root / "data" / "raw" / "pdb_chain_uniprot_current.csv.gz"
    else:
        sifts_path = Path(args.sifts)
        if not sifts_path.is_absolute():
            sifts_path = (repo_root / sifts_path).resolve()

    if args.output is None:
        output_path = repo_root / "data" / "interim" / "intact_pairs_with_pdb_templates.tsv"
    else:
        output_path = Path(args.output)
        if not output_path.is_absolute():
            output_path = (repo_root / output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    return intact_path, mmseqs_path, sifts_path, output_path


# ---------- SIFTS 読み込み ----------

def load_sifts_chain_uniprot(sifts_path: Path) -> pd.DataFrame:
    """
    SIFTS pdb_chain_uniprot_current.csv.gz を読み込み、

        chain_label (例: 1ABC_A)
        pdb_id      (例: 1ABC)
        chain_id    (例: A)
        uniprot_acc (UniProt ID)

    の 4列を持つ DataFrame を返す。
    """
    if not sifts_path.exists():
        print(f"[ERROR] SIFTS file not found: {sifts_path}", file=sys.stderr)
        sys.exit(1)

    sifts_df = pd.read_csv(sifts_path, comment="#", low_memory=False)

    print(f"[INFO] SIFTS columns: {list(sifts_df.columns)}")

    # 列名のパターンを判定
    col_map = {}
    if {"pdb_id", "chain_id", "uniprot_acc"}.issubset(sifts_df.columns):
        # 既に期待される名前のとき
        col_map = {}
    elif {"PDB", "CHAIN", "SP_PRIMARY"}.issubset(sifts_df.columns):
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
            f"[ERROR] リネーム後も {expected_cols} が揃っていません。"
            f" 実際の列: {list(sifts_df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    sifts_df["chain_label"] = [
        build_chain_label(pdb_id, chain_id)
        for pdb_id, chain_id in zip(sifts_df["pdb_id"], sifts_df["chain_id"])
    ]
    sifts_df = sifts_df.dropna(subset=["chain_label"])

    chain_uniprot = (
        sifts_df[["chain_label", "pdb_id", "chain_id", "uniprot_acc"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    print(f"[INFO] # SIFTS chain→UniProt rows : {len(chain_uniprot)}")
    print(f"[INFO] # unique SIFTS chain_label : {chain_uniprot['chain_label'].nunique()}")
    print(f"[INFO] # unique SIFTS UniProt     : {chain_uniprot['uniprot_acc'].nunique()}")

    return chain_uniprot


# ---------- mmseqs 読み込み & IntAct に関係するヒット抽出 ----------

def load_intact_pairs(intact_path: Path) -> pd.DataFrame:
    if not intact_path.exists():
        print(f"[ERROR] IntAct file not found: {intact_path}", file=sys.stderr)
        sys.exit(1)
    df = pd.read_csv(intact_path, sep="\t")
    if not {"uniprot_a", "uniprot_b"}.issubset(df.columns):
        print(
            "[ERROR] IntAct TSV に 'uniprot_a', 'uniprot_b' 列が必要です。",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"[INFO] # IntAct pairs          : {len(df)}")
    print(f"[INFO] # unique IntAct UniProt : {len(set(df['uniprot_a']) | set(df['uniprot_b']))}")
    return df


def load_mmseqs_hits_for_intact(
    mmseqs_path: Path,
    intact_uniprot_set: set[str],
    chain_uniprot: pd.DataFrame,
    min_pident: float,
    min_alnlen: int,
    max_evalue: float,
) -> pd.DataFrame:
    """
    mmseqs m8 ファイルをチャンクで読み込み、
    - query_uniprot が IntAct のノードに含まれる
    - pident, alnlen, evalue が閾値を満たす
    - target が SIFTS chain_label とマッチする
    ヒットだけを集めて返す。
    """
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

    print(f"[INFO] Loading mmseqs hits (chunked) from: {mmseqs_path}")
    print(
        f"[INFO] Filters: min_pident={min_pident}, "
        f"min_alnlen={min_alnlen}, max_evalue={max_evalue}"
    )

    chain_uniprot = chain_uniprot.copy()
    # join用に index を貼っておくと若干速い
    chain_uniprot = chain_uniprot.set_index("chain_label")

    chunks = []
    chunk_iter = pd.read_csv(
        mmseqs_path,
        sep="\t",
        header=None,
        names=m8_cols,
        chunksize=1_000_000,
    )

    total_rows = 0
    kept_rows = 0

    for i, chunk in enumerate(chunk_iter, start=1):
        total_rows += len(chunk)
        # query_uniprot に正規化
        chunk["query_uniprot"] = chunk["query"].map(normalize_uniprot_from_query)

        # IntAct に出てくるノードだけに絞る
        mask_intact = chunk["query_uniprot"].isin(intact_uniprot_set)

        # pident, alnlen, evalue でフィルタ
        mask_quality = (
            (chunk["pident"] >= min_pident * 100.0)  # pident は百分率と想定
            & (chunk["alnlen"] >= min_alnlen)
            & (chunk["evalue"] <= max_evalue)
        )

        sub = chunk[mask_intact & mask_quality].copy()
        if sub.empty:
            continue

        # SIFTS の chain_label によるマッピング
        # mmseqs の target は "1ABC_A" の形式を想定
        sub = sub.join(chain_uniprot, on="target", how="inner")

        if sub.empty:
            continue

        # mmseqs の target は SIFTS の chain_label と同じ表記なので、
        # 明示的に chain_label 列を付けておく
        sub["chain_label"] = sub["target"]

        kept_rows += len(sub)
        chunks.append(sub)

        print(
            f"[INFO]  chunk {i:3d}: read {len(chunk):7d} rows, "
            f"kept {len(sub):7d} (total kept: {kept_rows})"
        )

    if not chunks:
        print("[ERROR] No mmseqs hits remained after filtering.", file=sys.stderr)
        sys.exit(1)

    hits = pd.concat(chunks, ignore_index=True)
    print(f"[INFO] Total mmseqs rows read : {total_rows}")
    print(f"[INFO] Total mmseqs hits kept: {len(hits)}")

    # インデックスを戻しておく
    hits = hits.reset_index(drop=False).rename(columns={"index": "sifts_row_index"})

    return hits


def build_uniprot_to_chain_table(
    hits: pd.DataFrame,
    max_hits_per_uniprot: int,
) -> pd.DataFrame:
    """
    mmseqs + SIFTS のマージ結果 hits から、

        query_uniprot -> (pdb_id, chain_id, chain_label, uniprot_acc, pident, alnlen, evalue, bits)

    の対応を作り、
    - 同一 (query_uniprot, pdb_id, chain_id) で重複する行は
      最も evalue が良く bits の高いものだけ残す
    - さらに各 query_uniprot あたり max_hits_per_uniprot までに制限

    したテーブルを返す。
    """
    # まず (query_uniprot, pdb_id, chain_id) ごとにベストヒットを残す
    hits = hits.copy()
    # evalue が小さい順、bits が大きい順にソート
    hits["neg_bits"] = -hits["bits"]
    hits = hits.sort_values(
        ["query_uniprot", "pdb_id", "chain_id", "evalue", "neg_bits"],
        ascending=[True, True, True, True, True],
    )
    hits = hits.drop_duplicates(
        subset=["query_uniprot", "pdb_id", "chain_id"], keep="first"
    )

    # 各 query_uniprot あたり max_hits_per_uniprot まで
    hits = hits.sort_values(
        ["query_uniprot", "evalue", "neg_bits"],
        ascending=[True, True, True],
    )
    hits = hits.groupby("query_uniprot", as_index=False).head(max_hits_per_uniprot)

    # 必要な列だけ残す＆整理
    cols = [
        "query_uniprot",
        "pdb_id",
        "chain_id",
        "chain_label",
        "uniprot_acc",
        "pident",
        "alnlen",
        "evalue",
        "bits",
    ]
    hits = hits[cols].reset_index(drop=True)

    print(f"[INFO] # rows in uniprot→chain table: {len(hits)}")
    print(f"[INFO] # unique UniProt (with chains): {hits['query_uniprot'].nunique()}")

    return hits


# ---------- IntAct ペア × PDB テンプレート構築 ----------

def build_intact_pair_templates(
    intact_df: pd.DataFrame,
    uniprot_chain_table: pd.DataFrame,
    allow_self_chain: bool,
) -> pd.DataFrame:
    """
    IntAct ペア (uniprot_a, uniprot_b) と、
    UniProt→PDBチェーン対応 uniprot_chain_table を使って、

    IntAct ペアごとに「同じ PDB 内に chain_a と chain_b の両方が存在する」テンプレートを列挙する。
    """
    # uniprot_chain_table を A/B 用に複製＆リネーム
    a_map = uniprot_chain_table.rename(
        columns={
            "query_uniprot": "uniprot_a",
            "pdb_id": "pdb_id",
            "chain_id": "chain_id_a",
            "chain_label": "chain_label_a",
            "uniprot_acc": "sifts_uniprot_a",
            "pident": "pident_a",
            "alnlen": "alnlen_a",
            "evalue": "evalue_a",
            "bits": "bits_a",
        }
    )

    b_map = uniprot_chain_table.rename(
        columns={
            "query_uniprot": "uniprot_b",
            "pdb_id": "pdb_id",
            "chain_id": "chain_id_b",
            "chain_label": "chain_label_b",
            "uniprot_acc": "sifts_uniprot_b",
            "pident": "pident_b",
            "alnlen": "alnlen_b",
            "evalue": "evalue_b",
            "bits": "bits_b",
        }
    )

    print("[INFO] Merging IntAct pairs with PDB chains for A-side ...")
    df_a = intact_df.merge(a_map, on="uniprot_a", how="inner")
    print(f"[INFO]  After A-side merge: {len(df_a)} rows")

    print("[INFO] Merging with PDB chains for B-side (same PDB ID) ...")
    df_ab = df_a.merge(b_map, on=["uniprot_b", "pdb_id"], how="inner")
    print(f"[INFO]  After B-side merge: {len(df_ab)} rows")

    # 自己相互作用 (uniprot_a == uniprot_b) で、同一 chain_id を A/B 両方に割り当てた行は
    # デフォルトでは除外する（同じ鎖を2回使うのは実際の複合体にならないため）
    if not allow_self_chain:
        mask_bad_self = (
            (df_ab["uniprot_a"] == df_ab["uniprot_b"])
            & (df_ab["chain_id_a"] == df_ab["chain_id_b"])
        )
        n_bad = mask_bad_self.sum()
        if n_bad > 0:
            print(
                f"[INFO] Removing {n_bad} rows where uniprot_a == uniprot_b "
                "and chain_id_a == chain_id_b (self-chain)."
            )
            df_ab = df_ab[~mask_bad_self]

    # 重複の整理（同じテンプレートが複数回出てくるのを防ぐ）
    df_ab = df_ab.drop_duplicates(
        subset=[
            "uniprot_a",
            "uniprot_b",
            "pdb_id",
            "chain_id_a",
            "chain_id_b",
        ]
    ).reset_index(drop=True)

    print(f"[INFO] # unique IntAct pair × PDB template rows: {len(df_ab)}")
    print(
        f"[INFO] # IntAct pairs with ≥1 template: "
        f"{df_ab[['uniprot_a','uniprot_b']].drop_duplicates().shape[0]}"
    )

    # 列の順序を整理
    cols_out = [
        "uniprot_a",
        "uniprot_b",
        "pdb_id",
        "chain_id_a",
        "chain_label_a",
        "sifts_uniprot_a",
        "pident_a",
        "alnlen_a",
        "evalue_a",
        "bits_a",
        "chain_id_b",
        "chain_label_b",
        "sifts_uniprot_b",
        "pident_b",
        "alnlen_b",
        "evalue_b",
        "bits_b",
    ]
    df_ab = df_ab[cols_out]

    # 同一 PDB chain を A/B に割り当てた行を除外
    mask_same_chain = df_ab["chain_id_a"] == df_ab["chain_id_b"]
    df_ab = df_ab[~mask_same_chain].copy()

    return df_ab


# ---------- main ----------

def main() -> None:
    args = parse_args()
    intact_path, mmseqs_path, sifts_path, output_path = resolve_paths(args)

    print(f"[INFO] IntAct pairs       : {intact_path}")
    print(f"[INFO] mmseqs results     : {mmseqs_path}")
    print(f"[INFO] SIFTS mapping      : {sifts_path}")
    print(f"[INFO] Output             : {output_path}")

    # 1. IntAct ペア
    intact_df = load_intact_pairs(intact_path)
    intact_uniprot_set = set(intact_df["uniprot_a"]) | set(intact_df["uniprot_b"])

    # 2. SIFTS (chain_label -> pdb_id, chain_id, uniprot_acc)
    chain_uniprot = load_sifts_chain_uniprot(sifts_path)

    # 3. mmseqs 結果から IntAct ノードに関係するヒットだけ抽出
    hits = load_mmseqs_hits_for_intact(
        mmseqs_path=mmseqs_path,
        intact_uniprot_set=intact_uniprot_set,
        chain_uniprot=chain_uniprot,
        min_pident=args.min_pident,
        min_alnlen=args.min_alnlen,
        max_evalue=args.max_evalue,
    )

    # 4. UniProt -> PDBチェーン対応表を構築
    uniprot_chain_table = build_uniprot_to_chain_table(
        hits,
        max_hits_per_uniprot=args.max_hits_per_uniprot,
    )

    # 5. IntAct ペアごとの PDB テンプレート候補を列挙
    pair_templates = build_intact_pair_templates(
        intact_df=intact_df,
        uniprot_chain_table=uniprot_chain_table,
        allow_self_chain=args.allow_self_chain,
    )

    # 6. TSV に出力
    pair_templates.to_csv(output_path, sep="\t", index=False)
    print(f"[INFO] Written: {output_path}")


if __name__ == "__main__":
    main()
