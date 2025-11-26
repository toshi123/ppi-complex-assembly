#!/usr/bin/env python
# scripts/select_template_topk_diverse.py
"""
多様性考慮 Top-K 代表テンプレート選出

入力: filter_pdb_contacts_parallel.py の出力TSV
 必須列:
   - uniprot_a, uniprot_b, pdb_id, chain_id_a, chain_id_b
   - min_atom_distance, n_atom_contacts
   - iface_res_a, iface_res_b   # カンマ区切りの残基ID

任意列（あれば使う）:
   - has_ligand (bool), ligand_class (str)
   - resol, pident_mean, assembly_id など

特徴:
 - base_score = (-min_atom_distance) + log1p(n_atom_contacts) (+任意の副次スコア)
 - MMR（Maximal Marginal Relevance）で多様性を確保
   類似度: インターフェース残基Jaccard + (任意) リガンド状態一致, 同一PDB/assembly

使い方:
  python scripts/select_template_topk_diverse.py \
    --input  data/interim/intact_pairs_with_pdb_contacts.tsv \
    --output data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
    --k 3 --lambda 0.5 --min-iface-diff 0.25
"""

import argparse
import math
import pandas as pd
import numpy as np

def jaccard(a:set, b:set) -> float:
    if not a and not b: return 1.0
    if not a or not b:  return 0.0
    return len(a & b) / len(a | b)

def parse_set(s):
    if pd.isna(s) or s == "": return set()
    return set(str(s).split(","))

def compute_base_score(row, w_dist=1.0, w_contacts=1.0, w_pident=0.0, w_resol=0.0):
    dist = float(row.get("min_atom_distance", 10.0))
    ncnt = float(row.get("n_atom_contacts", 0.0))
    pid  = float(row.get("pident_mean", 0.0)) if not pd.isna(row.get("pident_mean", np.nan)) else 0.0
    resol = float(row.get("resol", 0.0)) if not pd.isna(row.get("resol", np.nan)) else 0.0
    # 小さい距離＆多い接触を優先。分解能は小さいほど良いのでマイナス符号。
    score = w_dist * (-min(dist, 10.0)) + w_contacts * math.log1p(max(ncnt, 0.0)) \
            + w_pident * pid - w_resol * resol
    return score

def mmr_select(rows, K=3, lam=0.5,
               w_iface=1.0, w_lig=0.3, w_same_pdb=0.2, w_same_asm=0.2,
               min_iface_diff=0.0):
    """
    rows: list[dict] 同一 (uniprot_a, uniprot_b) の候補
    lam: 多様性重み（大きいほど多様性重視）
    min_iface_diff: 0.0〜1.0。選出済みとJaccardが高すぎる候補を弾くための下限距離（=1-類似度）。
                    例: 0.25 を指定 → 類似度 0.75 超は落とす（インターフェースがほぼ同じなら弾く）
    """
    # base_score で1つ目
    rows = rows.copy()
    rows.sort(key=lambda r: r["base_score"], reverse=True)
    chosen = []
    if rows:
        chosen.append(rows.pop(0))

    # 2つ目以降: base_score - lam * sim_max を最大化
    def lig_key(r):
        return (r.get("has_ligand", False), r.get("ligand_class", "NA"))

    while len(chosen) < K and rows:
        best = None
        best_val = -1e9
        for c in rows:
            # 類似度（最大）を計算
            sim_max = 0.0
            too_similar = False
            for s in chosen:
                ia = jaccard(c["iface_a"], s["iface_a"])
                ib = jaccard(c["iface_b"], s["iface_b"])
                iface_sim = (ia + ib) / 2.0
                # 類似しすぎを強制的に弾く（多様性確保）
                if min_iface_diff > 0.0 and (1.0 - iface_sim) < min_iface_diff:
                    too_similar = True
                    break
                lig_sim = 1.0 if lig_key(c) == lig_key(s) else 0.0
                same_pdb = 1.0 if str(c["pdb_id"]).upper()==str(s["pdb_id"]).upper() else 0.0
                same_asm = 1.0 if c.get("assembly_id") and c.get("assembly_id")==s.get("assembly_id") else 0.0
                sim = w_iface*iface_sim + w_lig*lig_sim + w_same_pdb*same_pdb + w_same_asm*same_asm
                sim_max = max(sim_max, sim)
            if too_similar:
                continue
            val = c["base_score"] - lam * sim_max
            if val > best_val:
                best_val, best = val, c
        if best is None:
            break
        chosen.append(best)
        rows.remove(best)
    return chosen

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="contacts TSV")
    ap.add_argument("--output", required=True, help="topK diverse TSV")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--lambda", dest="lam", type=float, default=0.5)
    ap.add_argument("--min-iface-diff", type=float, default=0.0, help="0〜1。大きいほど多様性を強制")
    # base_score の重み
    ap.add_argument("--w-dist", type=float, default=1.0)
    ap.add_argument("--w-contacts", type=float, default=1.0)
    ap.add_argument("--w-pident", type=float, default=0.0)
    ap.add_argument("--w-resol", type=float, default=0.0)
    # 類似度の重み
    ap.add_argument("--w-iface", type=float, default=1.0)
    ap.add_argument("--w-lig", type=float, default=0.3)
    ap.add_argument("--w-same-pdb", type=float, default=0.2)
    ap.add_argument("--w-same-asm", type=float, default=0.2)
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # 必須列チェック
    need = {"uniprot_a","uniprot_b","pdb_id","chain_id_a","chain_id_b",
            "min_atom_distance","n_atom_contacts","iface_res_a","iface_res_b"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"Missing required columns: {miss}")

    # セット化
    df = df.copy()
    df["iface_a"] = df["iface_res_a"].map(parse_set)
    df["iface_b"] = df["iface_res_b"].map(parse_set)

    # base_score
    df["base_score"] = df.apply(
        lambda r: compute_base_score(
            r, w_dist=args.w_dist, w_contacts=args.w_contacts,
            w_pident=args.w_pident, w_resol=args.w_resol
        ),
        axis=1
    )

    out_rows = []
    for (ua, ub), g in df.groupby(["uniprot_a","uniprot_b"], sort=False):
        rows = g.to_dict(orient="records")
        chosen = mmr_select(
            rows, K=args.k, lam=args.lam,
            w_iface=args.w_iface, w_lig=args.w_lig,
            w_same_pdb=args.w_same_pdb, w_same_asm=args.w_same_asm,
            min_iface_diff=args.min_iface_diff
        )
        out_rows.extend(chosen)

    out_df = pd.DataFrame(out_rows)
    # 並べ替え（対称重複を避けるため一応正規化キー付け、必要なら解除）
    out_df["sort_key"] = out_df["uniprot_a"] + "_" + out_df["uniprot_b"] + "_" + out_df["pdb_id"].astype(str)
    out_df = out_df.sort_values(["uniprot_a","uniprot_b","base_score"], ascending=[True,True,False]).drop(columns=["sort_key"])
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] written: {args.output}  rows={len(out_df)}")

if __name__ == "__main__":
    main()
