#!/usr/bin/env python
import argparse
from pathlib import Path
import sys
import pandas as pd

PAIR_KEY = ["pdb_id", "chain_id_a", "chain_id_b"]

# 列名候補辞書
UA_CAND = ["uniprot_a","UniProt_A","A_uniprot","uniprot_id_a",
           "uniprotA","uniprot_acc_a","uniprotA_acc","query_uniprot","A_acc"]
UB_CAND = ["uniprot_b","UniProt_B","B_uniprot","uniprot_id_b",
           "uniprotB","uniprot_acc_b","uniprotB_acc","target_uniprot","B_acc"]

def pick(df: pd.DataFrame, names):
    for n in names:
        if n in df.columns:
            return n
    return None

def need_keys(df: pd.DataFrame):
    missing = [k for k in PAIR_KEY if k not in df.columns]
    if missing:
        raise SystemExit(f"[ERROR] missing key columns {missing} in this table. "
                         f"Have: {list(df.columns)}")

def load_df(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", required=True, help="band TSV needing uniprot columns")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--candidates", nargs="+", required=True,
                    help="upstream TSVs to source uniprot_a/b from (ordered)")
    args = ap.parse_args()

    band_p = Path(args.band)
    out_p  = Path(args.out)
    need_srcs = [Path(p) for p in args.candidates]

    band = load_df(band_p)
    need_keys(band)

    # 既に uniprot が band にあればそのまま保存
    ua = pick(band, UA_CAND)
    ub = pick(band, UB_CAND)
    if ua and ub:
        # 標準名に正規化
        band = band.rename(columns={ua:"uniprot_a", ub:"uniprot_b"})
        band.to_csv(out_p, sep="\t", index=False)
        print(f"[INFO] band already had uniprot columns -> normalized & written: {out_p}")
        return

    # 上流候補から復元
    merged = None
    for srcp in need_srcs:
        if not srcp.exists():
            continue
        src = load_df(srcp)
        try:
            need_keys(src)
        except SystemExit:
            # 上流が key を欠いていたらスキップ
            continue

        sua = pick(src, UA_CAND)
        sub = pick(src, UB_CAND)
        if not (sua and sub):
            continue

        add = src[PAIR_KEY + [sua, sub]].drop_duplicates(PAIR_KEY).rename(
            columns={sua:"uniprot_a", sub:"uniprot_b"}
        )

        tmp = band.merge(add, on=PAIR_KEY, how="left", validate="m:1")
        missing = tmp["uniprot_a"].isna().sum() + tmp["uniprot_b"].isna().sum()
        print(f"[INFO] tried {srcp.name}: filled={len(tmp)-missing//2} rows, missing_fields={missing}")

        # 採用条件：少なくとも一部でも埋まったら採用し、残欠損があれば次の候補で上書きマージ
        if merged is None:
            merged = tmp
        else:
            # 既存欠損に限って新値で埋める
            for col in ("uniprot_a","uniprot_b"):
                mask = merged[col].isna() & tmp[col].notna()
                merged.loc[mask, col] = tmp.loc[mask, col]

        # 両列がすべて埋まったら終了
        if merged["uniprot_a"].notna().all() and merged["uniprot_b"].notna().all():
            break

    if merged is None:
        raise SystemExit("[ERROR] could not find any upstream file with uniprot columns.")

    # 最終チェック
    still_miss = merged["uniprot_a"].isna().sum() + merged["uniprot_b"].isna().sum()
    if still_miss > 0:
        print(f"[WARN] some rows still missing uniprot columns: {still_miss} fields total", file=sys.stderr)

    out_p.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_p, sep="\t", index=False)
    print(f"[INFO] written: {out_p}  rows={len(merged)}")

if __name__ == "__main__":
    main()
