#!/usr/bin/env python
# scripts/qc_representatives.py
"""
Top-K 代表テンプレートのQCスクリプト

主な機能
- 列の自動補完: min_atom_distance / n_atom_contacts / bsa_total が無ければ
  指定した upstream TSV 群から (uniprot_a,b,pdb_id,chain_id_a,b) キーでマージ補完
- カバレッジ: 各 IntAct ペアの選抜件数 (n_sel) 分布と <K の件数
- 多様性: iface_source (sasa/dist) の内訳、pairごとの distinct-PDB 比率
- 代表の健全性: 距離>max_dist or 接触<min_contacts の “疑わしい代表” を検出
- 要約の保存: summary.json / pairs_counts.tsv / suspicious.tsv などを --out-prefix で保存

使い方例:
  python scripts/qc_representatives.py \
    --topk data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
    --k 3 \
    --upstream data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.fixed.tsv \
               data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv \
               data/interim/intact_pairs_with_pdb_contacts.tight_bsa.sane.tsv \
               data/interim/intact_pairs_with_pdb_contacts.tight_bsa.fixed.tsv \
               data/interim/intact_pairs_with_pdb_contacts.tight.tsv \
               data/interim/intact_pairs_with_pdb_templates.tsv \
    --max-dist 4.5 --min-contacts 20 \
    --out-prefix data/interim/qc_topk

出力物 (例):
  data/interim/qc_topk.summary.json
  data/interim/qc_topk.pairs_counts.tsv
  data/interim/qc_topk.suspicious.tsv
"""

import argparse
import json
from pathlib import Path
from typing import List

import pandas as pd


NEED_COLS = ["min_atom_distance", "n_atom_contacts", "bsa_total"]
KEY_COLS  = ["uniprot_a","uniprot_b","pdb_id","chain_id_a","chain_id_b"]

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--topk", required=True, help="Top-K TSV (代表出力)")
    ap.add_argument("--upstream", nargs="*", default=[],
                    help="不足列の補完に使うTSV群（上から優先）")
    ap.add_argument("--k", type=int, default=3, help="K（各ペアの目標件数）")
    ap.add_argument("--max-dist", type=float, default=4.5,
                    help="疑わしい代表: min_atom_distance > この値")
    ap.add_argument("--min-contacts", type=int, default=20,
                    help="疑わしい代表: n_atom_contacts < この値")
    ap.add_argument("--out-prefix", type=str, default="",
                    help="保存プレフィックス（空なら保存しない）")
    return ap.parse_args()


def _safe_read_tsv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        return None


def backfill_columns(df: pd.DataFrame, upstream_paths: List[Path]) -> pd.DataFrame:
    """TopK df に必要列が無ければ upstream から補完（m:1でマージ）。"""
    missing = [c for c in NEED_COLS if c not in df.columns]
    if not missing:
        return df

    for p in upstream_paths:
        src = _safe_read_tsv(p)
        if src is None:
            continue
        if not set(KEY_COLS).issubset(src.columns):
            continue
        have = [c for c in missing if c in src.columns]
        if not have:
            continue

        add = src[KEY_COLS + have].drop_duplicates(KEY_COLS)
        df  = df.merge(add, on=KEY_COLS, how="left", validate="m:1")
        # まだ欠けている列を更新
        missing = [c for c in NEED_COLS if c not in df.columns]
        if not missing:
            break

    # なお欠けている列は NaN のダミーを作る（QCは走らせる）
    for c in missing:
        df[c] = pd.NA
    return df


def main():
    args = parse_args()
    topk_path = Path(args.topk)
    df = _safe_read_tsv(topk_path)
    if df is None or df.empty:
        raise SystemExit(f"[ERROR] cannot read topk file or empty: {topk_path}")

    # 入力健全性確認（キー列）
    if not set(KEY_COLS).issubset(df.columns):
        raise SystemExit(f"[ERROR] missing key columns in topk: need {KEY_COLS}, have {list(df.columns)}")

    # 不足列の補完
    upstream_paths = [Path(p) for p in (args.upstream or [])]
    df = backfill_columns(df, upstream_paths)

    # 1) カバレッジ
    pairs = df.groupby(["uniprot_a","uniprot_b"]).size().rename("n_sel").reset_index()
    n_pairs_total = pairs.shape[0]
    lt_k = (pairs["n_sel"] < args.k).sum()
    n_sel_count = pairs["n_sel"].value_counts().sort_index().to_dict()

    # 2) 多様性: iface_source, distinct-PDB 比率
    iface_counts = {}
    if "iface_source" in df.columns:
        iface_counts = df["iface_source"].value_counts(dropna=False).to_dict()

    ratio = None
    try:
        per_pair_nuniq = (df.groupby(["uniprot_a","uniprot_b"])["pdb_id"].nunique())
        per_pair_size  = (df.groupby(["uniprot_a","uniprot_b"]).size())
        ratio = float((per_pair_nuniq / per_pair_size).mean())
    except Exception:
        ratio = None

    # 3) 代表の健全性: suspicious
    md = pd.to_numeric(df.get("min_atom_distance"), errors="coerce")
    nc = pd.to_numeric(df.get("n_atom_contacts"), errors="coerce")
    suspicious_mask = pd.Series([False]*len(df))
    if "min_atom_distance" in df.columns:
        suspicious_mask |= (md > float(args.max_dist))
    if "n_atom_contacts" in df.columns:
        suspicious_mask |= (nc < int(args.min_contacts))
    suspicious = df[suspicious_mask].copy()

    # 4) 基本統計（ある列だけ）
    stats_cols = [c for c in ["bsa_total","final_score","min_atom_distance","n_atom_contacts"] if c in df.columns]
    stats_desc = df[stats_cols].describe().to_dict() if stats_cols else {}

    # 5) 出力
    summary = {
        "input_rows": int(len(df)),
        "n_pairs": int(n_pairs_total),
        "k": int(args.k),
        "n_sel_hist": n_sel_count,
        "pairs_with_<k": int(lt_k),
        "iface_source_counts": iface_counts,
        "avg_distinct_PDB_per_pair_ratio": ratio,
        "suspicious_count": int(len(suspicious)),
        "stats": stats_desc,
        "thresholds": {
            "max_distance": float(args.max_dist),
            "min_contacts": int(args.min_contacts),
        },
    }

    # ログ表示
    print("=== QC Summary ===")
    print(json.dumps(summary, indent=2, ensure_ascii=False))

    # 保存
    if args.out_prefix:
        out_prefix = Path(args.out_prefix)
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

        # サマリ
        (out_prefix.with_suffix(".summary.json")).write_text(
            json.dumps(summary, indent=2, ensure_ascii=False)
        )

        # ペアごとの選抜件数
        pairs.sort_values(["n_sel","uniprot_a","uniprot_b"], inplace=True)
        pairs.to_csv(out_prefix.with_suffix(".pairs_counts.tsv"), sep="\t", index=False)

        # 疑わしい代表
        if len(suspicious):
            keep_cols = list({*KEY_COLS, "rank", "min_atom_distance", "n_atom_contacts", "bsa_total", "final_score"} & set(df.columns))
            suspicious = suspicious[keep_cols] if keep_cols else suspicious
            suspicious.sort_values(["min_atom_distance","n_atom_contacts"], ascending=[False, True], inplace=True)
            suspicious.to_csv(out_prefix.with_suffix(".suspicious.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
