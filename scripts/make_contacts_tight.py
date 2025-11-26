#!/usr/bin/env python
# scripts/make_contacts_tight.py
"""
make_contacts_tight.py

filter_pdb_contacts_parallel.py の出力（contacts TSV）を読み込み、
“tight” なしきい値でノイズを間引いたTSVを作成します。

既定のフィルタ（BSA未導入の暫定版）:
  - n_atom_contacts <= 1500          # 過剰に巨大な界面を抑制
  - n_res_contacts_a >= 10
  - n_res_contacts_b >= 10           # 点接触を除外
  - 1.6 <= min_atom_distance <= 4.5  # 異常に近すぎ/遠すぎを除外
  - contacts_per_res <= 40           # 過密すぎる界面を抑制
  - （任意）--only-hetero で同一タンパク（homodimer）を除外

必要に応じてコマンドラインで閾値を変更できます。
"""

import argparse
import sys
from pathlib import Path
import pandas as pd

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  required=True,
                    help="contacts TSV（filter_pdb_contacts_parallel.py の出力）")
    ap.add_argument("--output", default=None,
                    help="出力パス（省略時は <input>.tight.tsv）")
    ap.add_argument("--max-atom-contacts", type=int, default=1500)
    ap.add_argument("--min-res-a", type=int, default=10)
    ap.add_argument("--min-res-b", type=int, default=10)
    ap.add_argument("--min-distance", type=float, default=1.6)
    ap.add_argument("--max-distance", type=float, default=4.5)
    ap.add_argument("--max-contacts-per-res", type=float, default=40.0)
    ap.add_argument("--only-hetero", action="store_true",
                    help="uniprot_a != uniprot_b のみ残す（自己相互作用を除外）")
    return ap.parse_args()

def main():
    args = parse_args()
    inp = Path(args.input)
    if not inp.exists():
        print(f"[ERROR] input not found: {inp}", file=sys.stderr); sys.exit(1)
    out = Path(args.output) if args.output else inp.with_suffix(".tight.tsv")

    df = pd.read_csv(inp, sep="\t")
    n0 = len(df)

    # 必須列のざっくりチェック
    need = {"uniprot_a","uniprot_b","pdb_id","chain_id_a","chain_id_b",
            "min_atom_distance","n_atom_contacts","n_res_contacts_a","n_res_contacts_b"}
    miss = need - set(df.columns)
    if miss:
        print(f"[ERROR] missing columns: {sorted(miss)}", file=sys.stderr); sys.exit(1)

    # 過密度
    denom = (df["n_res_contacts_a"] + df["n_res_contacts_b"]).clip(lower=1)
    df["contacts_per_res"] = df["n_atom_contacts"] / denom

    # フィルタ適用
    mask = (
        (df["n_atom_contacts"] <= args.max_atom_contacts) &
        (df["n_res_contacts_a"] >= args.min_res_a) &
        (df["n_res_contacts_b"] >= args.min_res_b) &
        (df["min_atom_distance"] >= args.min_distance) &
        (df["min_atom_distance"] <= args.max_distance) &
        (df["contacts_per_res"] <= args.max_contacts_per_res)
    )
    if args.only_hetero:
        mask &= (df["uniprot_a"] != df["uniprot_b"])

    kept = df[mask].copy()
    kept.to_csv(out, sep="\t", index=False)

    # 簡易サマリ
    def uniq_pairs(x): return x[["uniprot_a","uniprot_b"]].drop_duplicates().shape[0]
    print("[INFO] thresholds:",
          f"max_atom_contacts<={args.max_atom_contacts},",
          f"min_res_a>={args.min_res_a}, min_res_b>={args.min_res_b},",
          f"{args.min_distance}<=min_atom_distance<={args.max_distance},",
          f"contacts_per_res<={args.max_contacts_per_res},",
          f"only_hetero={args.only_hetero}")
    print(f"[INFO] input rows  : {n0}")
    print(f"[INFO] kept rows   : {len(kept)}")
    print(f"[INFO] uniq pairs  : {uniq_pairs(kept)}")
    print(f"[INFO] uniq PDB ids: {kept['pdb_id'].nunique()}")
    print(f"[INFO] written     : {out}")

if __name__ == "__main__":
    main()
