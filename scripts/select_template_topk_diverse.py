#!/usr/bin/env python
# scripts/select_template_topk_diverse.py
"""
各 IntAct ペア (uniprot_a, uniprot_b) について、テンプレート候補行から
MMR(Maximal Marginal Relevance)で Top-K を選抜する。

特徴:
- 界面集合は --use-sasa-iface で ΔSASA 由来 (iface_sasa_*) を優先
  (欠損時は距離由来 iface_res_* へフォールバック)
- ベーススコア = w_dist*DistScore + w_contacts*ContactScore + w_bsa*BsaScore
- 多様性: Jaccard類似度を用い MMR で選抜（同一PDB/Assembly/鎖組の小ペナ追加）
- 出力にスコア内訳と選抜順位を付与

必須列:
  uniprot_a, uniprot_b, pdb_id, chain_id_a, chain_id_b
推奨列:
  min_atom_distance, n_atom_contacts
  bsa_total (w_bsa>0なら)
  iface_sasa_a, iface_sasa_b  または  iface_res_a, iface_res_b
任意列:
  assembly_id（あれば多様性ペナルティに使用）
  has_ligand（0/1, True/False; w_lig で軽い多様性制御に使用可）

使い方(例):
  python scripts/select_template_topk_diverse.py \
    --input  data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv \
    --output data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
    --k 3 --lambda 0.5 --min-iface-diff 0.25 \
    --use-sasa-iface \
    --w-dist 1.0 --w-contacts 1.0 --w-bsa 0.5 \
    --w-iface 1.0 --w-lig 0.3 --w-same-pdb 0.2 --w-same-asm 0.2
"""

import argparse
import math
from pathlib import Path
from typing import List, Set, Tuple, Optional

import numpy as np
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--k", type=int, default=3, help="Top-K per pair")
    ap.add_argument("--lambda", dest="lam", type=float, default=0.5,
                    help="MMR lambda (0=only diversity, 1=only quality)")
    ap.add_argument("--min-iface-diff", type=float, default=0.25,
                    help="同一ペア内で別テンプレートとして扱う最小Jaccard距離(=1-Jaccard)")
    ap.add_argument("--use-sasa-iface", action="store_true",
                    help="ΔSASA由来 iface_sasa_* を優先使用（欠損時は距離由来にフォールバック）")

    # base score weights
    ap.add_argument("--w-dist", type=float, default=1.0,
                    help="小さい距離ほど高スコア。1.6–4.5Åを中心に正規化")
    ap.add_argument("--w-contacts", type=float, default=1.0,
                    help="n_atom_contacts の log1p 正規化スコア重み")
    ap.add_argument("--w-bsa", type=float, default=0.0,
                    help="bsa_total の log1p 正規化スコア重み（ΔSASA付きデータ推奨）")

    # diversity penalty weights
    ap.add_argument("--w-iface", type=float, default=1.0,
                    help="界面集合のJaccard類似ペナルティ重み")
    ap.add_argument("--w-lig", type=float, default=0.3,
                    help="同一リガンド有無の軽ペナルティ（has_ligand列があれば使用）")
    ap.add_argument("--w-same-pdb", type=float, default=0.2,
                    help="同一PDBの軽ペナルティ")
    ap.add_argument("--w-same-asm", type=float, default=0.2,
                    help="同一assembly_idの軽ペナルティ（列がある場合のみ）")
    return ap.parse_args()


# ------------------------ helpers ------------------------

def _safe_set_from_str(s: str) -> Set[str]:
    if pd.isna(s) or s == "":
        return set()
    return set(str(s).split(","))


def _jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return 0.0 if union == 0 else inter / union


def _norm01(series: pd.Series) -> pd.Series:
    # robust min-max with fallback
    s = pd.to_numeric(series, errors="coerce")
    vmin, vmax = s.min(), s.max()
    if pd.isna(vmin) or pd.isna(vmax) or vmax <= vmin:
        return pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    return (s - vmin) / (vmax - vmin)


def _dist_to_score(d: pd.Series) -> pd.Series:
    """
    小さいほどよい距離を [0,1] のスコアに変換。
    目安: 1.6Å(最良)〜4.5Å(許容下限)で線形に落ちるクリップ。
    """
    x = pd.to_numeric(d, errors="coerce")
    s = 1.0 - (x - 1.6) / (4.5 - 1.6)
    s = s.clip(lower=0.0, upper=1.0).fillna(0.0)
    return s


def _log1p_norm(x: pd.Series) -> pd.Series:
    s = pd.to_numeric(x, errors="coerce").fillna(0)
    s = np.log1p(s)
    return _norm01(pd.Series(s, index=x.index))


def _choose_iface_sets(row, use_sasa: bool) -> Tuple[Set[str], Set[str], str]:
    src = "dist"
    if use_sasa:
        ia = _safe_set_from_str(row.get("iface_sasa_a", ""))
        ib = _safe_set_from_str(row.get("iface_sasa_b", ""))
        if ia or ib:
            return ia, ib, "sasa"
        # フォールバック
    ia = _safe_set_from_str(row.get("iface_res_a", ""))
    ib = _safe_set_from_str(row.get("iface_res_b", ""))
    return ia, ib, src


def _mmr_select(group_df: pd.DataFrame,
                k: int,
                lam: float,
                min_iface_diff: float,
                w_iface: float,
                w_lig: float,
                w_same_pdb: float,
                w_same_asm: float,
                use_sasa_iface: bool) -> pd.DataFrame:
    """
    1ペアの候補から Top-K をMMRで選ぶ。group_df は同一 (uniprot_a, uniprot_b) のみ。
    返り値は選抜行のみ（rank/score列付き）。
    """
    if group_df.empty:
        return group_df

    # 事前に界面集合とソースを作る
    iface_a, iface_b, srcs = [], [], []
    for _, r in group_df.iterrows():
        ia, ib, src = _choose_iface_sets(r, use_sasa_iface)
        iface_a.append(ia); iface_b.append(ib); srcs.append(src)

    G = group_df.copy()
    G["_iface_a"] = iface_a
    G["_iface_b"] = iface_b
    G["_iface_src"] = srcs

    # 類似度（Jaccard）はA/B両側の平均で代表
    def iface_sim(i, j) -> float:
        ja = _jaccard(G["_iface_a"].iat[i], G["_iface_a"].iat[j])
        jb = _jaccard(G["_iface_b"].iat[i], G["_iface_b"].iat[j])
        return 0.5 * (ja + jb)

    # 多様性“追加”ペナルティ（同一PDB/同一Asm/同一Lig）
    def extra_pen(i, j) -> float:
        pen = 0.0
        if w_same_pdb != 0:
            pen += w_same_pdb * (G["pdb_id"].iat[i] == G["pdb_id"].iat[j])
        if w_same_asm != 0 and "assembly_id" in G.columns:
            a = G["assembly_id"].iat[i]; b = G["assembly_id"].iat[j]
            pen += w_same_asm * (pd.notna(a) and pd.notna(b) and a == b)
        if w_lig != 0 and "has_ligand" in G.columns:
            a = str(G["has_ligand"].iat[i]); b = str(G["has_ligand"].iat[j])
            pen += w_lig * (a == b and a != "nan")
        return float(pen)

    # MMR本体
    selected_idx: List[int] = []
    avail = list(range(len(G)))

    # 先頭は base_score 最大を選ぶ
    first = int(np.argmax(G["base_score"].values))
    selected_idx.append(first)
    avail.remove(first)

    while avail and len(selected_idx) < k:
        best_cand = None
        best_score = -1e9
        for idx in avail:
            # 類似度（=1-距離）と追加ペナ
            sim = 0.0
            for j in selected_idx:
                sim += w_iface * iface_sim(idx, j) + extra_pen(idx, j)
            if selected_idx:
                sim /= len(selected_idx)
            # MMR: λ*quality − (1−λ)*similarity
            mmr = lam * G["base_score"].iat[idx] - (1.0 - lam) * sim
            if mmr > best_score:
                best_score = mmr
                best_cand = idx
        selected_idx.append(best_cand)
        avail.remove(best_cand)

    S = G.iloc[selected_idx].copy()
    # rankと最終スコア（再計算して付与）
    ranks, finals, divs, src = [], [], [], []
    for rpos, ridx in enumerate(selected_idx, start=1):
        sim = 0.0
        for j in selected_idx[:rpos-1]:
            sim += w_iface * iface_sim(ridx, j) + extra_pen(ridx, j)
        if rpos > 1:
            sim /= (rpos - 1)
        mmr = lam * G["base_score"].iat[ridx] - (1.0 - lam) * sim
        ranks.append(rpos); finals.append(mmr); divs.append(sim); src.append(G["_iface_src"].iat[ridx])

    S["rank"] = ranks
    S["final_score"] = finals
    S["div_penalty"] = divs
    S["iface_source"] = src

    # オプション: ほぼ同一界面（Jaccard距離 < min_iface_diff）は同一扱い → 下位を落とす
    keep_mask = [True]*len(S)
    for i in range(len(S)):
        if not keep_mask[i]:
            continue
        for j in range(i+1, len(S)):
            if not keep_mask[j]:
                continue
            # 距離 = 1 - 類似度
            ja = _jaccard(S["_iface_a"].iat[i], S["_iface_a"].iat[j])
            jb = _jaccard(S["_iface_b"].iat[i], S["_iface_b"].iat[j])
            dist = 1.0 - 0.5*(ja+jb)
            if dist < min_iface_diff:
                # 類似過多 → 下位（スコア低い方）を落とす
                if S["final_score"].iat[i] >= S["final_score"].iat[j]:
                    keep_mask[j] = False
                else:
                    keep_mask[i] = False
                    break
    S = S.loc[keep_mask].sort_values(["final_score","base_score"], ascending=[False,False])
    return S.head(k)


def main():
    args = parse_args()
    inp = Path(args.input)
    outp = Path(args.output)

    df = pd.read_csv(inp, sep="\t", low_memory=False)
    need = {"uniprot_a","uniprot_b","pdb_id","chain_id_a","chain_id_b"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"[ERROR] missing columns: {sorted(miss)}")

    # ベーススコアの構築
    # 距離
    dist_score = _dist_to_score(df.get("min_atom_distance", pd.Series([np.nan]*len(df))))
    dist_score *= float(args.w_dist)

    # 接触
    contacts_score = _log1p_norm(df.get("n_atom_contacts", pd.Series([0]*len(df))))
    contacts_score *= float(args.w_contacts)

    # BSA
    bsa_col = df.get("bsa_total", pd.Series([np.nan]*len(df)))
    bsa_score = _log1p_norm(bsa_col.fillna(0))
    bsa_score *= float(args.w_bsa)

    base = dist_score.add(contacts_score, fill_value=0).add(bsa_score, fill_value=0)
    df["base_score"] = base.fillna(0.0)

    # 1ペアごとのMMR選抜
    groups = df.groupby(["uniprot_a","uniprot_b"], sort=False)
    picks = []
    for (ua, ub), g in groups:
        sel = _mmr_select(
            g, k=int(args.k), lam=float(args.lam),
            min_iface_diff=float(args.min_iface_diff),
            w_iface=float(args.w_iface),
            w_lig=float(args.w_lig),
            w_same_pdb=float(args.w_same_pdb),
            w_same_asm=float(args.w_same_asm),
            use_sasa_iface=bool(args.use_sasa_iface)
        )
        picks.append(sel)

    if picks:
        out = pd.concat(picks, axis=0, ignore_index=True)
    else:
        out = df.iloc[0:0].copy()

    # 内部列を整える
    drop_tmp = ["_iface_a","_iface_b","_iface_src"]
    for c in drop_tmp:
        if c in out.columns:
            out.drop(columns=[c], inplace=True)

    out.sort_values(["uniprot_a","uniprot_b","rank"], inplace=True)
    outp.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outp, sep="\t", index=False)

    # ログ風の出力
    n_pairs = df[["uniprot_a","uniprot_b"]].drop_duplicates().shape[0]
    n_out_pairs = out[["uniprot_a","uniprot_b"]].drop_duplicates().shape[0]
    print(f"[INFO] input rows: {len(df)}  pairs: {n_pairs}")
    print(f"[INFO] written: {outp}  rows: {len(out)}  pairs: {n_out_pairs}")
    if "bsa_total" in out.columns:
        print(out[["bsa_total"]].describe())


if __name__ == "__main__":
    main()
