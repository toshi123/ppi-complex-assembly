#!/usr/bin/env python
"""
make_interface_json.py (upstream-aware)

TopKやbandなど、入力TSVに不足しがちな列（iface_res_*, min_atom_distance, n_atom_contacts）を
--upstream で与えた上流TSV群からキー結合で補完し、インターフェースJSONLを出力する。

使い方の例:
python scripts/make_interface_json.py \
  --input  data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
  --upstream data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv \
            data/interim/intact_pairs_with_pdb_contacts.tight_bsa.tsv \
            data/interim/intact_pairs_with_pdb_contacts.tight.tsv \
  --sifts   data/raw/pdb_chain_uniprot_current.csv.gz \
  --fasta   data/raw/uniprot_human_all_20251110.fasta \
  --output  data/processed/interfaces/interfaces.jsonl \
  --log-every 1000
"""

import argparse
import gzip
import json
import logging
from pathlib import Path
from collections import defaultdict
import pandas as pd

AA_SETS = {
    "hydrophobic": set("AILMVFWY"),
    "polar": set("STNQCY"),
    "charged_pos": set("KRH"),
    "charged_neg": set("DE"),
    "aromatic": set("FWY"),
}

NEEDED_BASE = {"uniprot_a","uniprot_b","pdb_id","chain_id_a","chain_id_b"}
NEEDED_DERIVED = {"bsa_total","n_atom_contacts","min_atom_distance"}
NEEDED_IFACE = {"iface_res_a","iface_res_b"}
CANDIDATE_IFACE_A = ["iface_res_a", "iface_sasa_res_a"]
CANDIDATE_IFACE_B = ["iface_res_b", "iface_sasa_res_b"]

def parse_args():
    p = argparse.ArgumentParser(description="Build interface JSONL; supplement missing cols from upstream TSVs.")
    p.add_argument("--input", required=True, help="Main TSV (e.g., topk_diverse.tsv or tight_bsa.band.tsv)")
    p.add_argument("--upstream", nargs="*", default=[], help="Zero or more TSVs to left-merge for missing columns (priority = order)")
    p.add_argument("--sifts", required=True, help="SIFTS chain mapping CSV (pdb_chain_uniprot_current.csv.gz)")
    p.add_argument("--fasta", required=True, help="UniProt FASTA (human)")
    p.add_argument("--output", default="data/processed/interfaces/interfaces.jsonl", help="Output JSONL")
    p.add_argument("--limit", type=int, default=None, help="Process first N rows of INPUT")
    p.add_argument("--log-every", type=int, default=2000)
    p.add_argument("--pident-threshold", type=float, default=90.0,
                   help="If pident >= threshold, use sifts_uniprot for mapping even when uniprot differs (default: 90.0)")
    return p.parse_args()


# Candidate column names for sifts_uniprot and pident
CANDIDATE_SIFTS_A = ["sifts_uniprot_a", "sifts_uniprot_a_x", "sifts_uniprot_a_y"]
CANDIDATE_SIFTS_B = ["sifts_uniprot_b", "sifts_uniprot_b_x", "sifts_uniprot_b_y"]
CANDIDATE_PIDENT_A = ["pident_a", "pident_a_x", "pident_a_y"]
CANDIDATE_PIDENT_B = ["pident_b", "pident_b_x", "pident_b_y"]

def extract_uniprot_ac(s: str) -> str:
    """Extract UniProt accession from formats like 'sp|P12345|NAME' or 'tr|A0A024QYR6|NAME' -> 'P12345' or 'A0A024QYR6'"""
    s = str(s).strip()
    if '|' in s:
        parts = s.split('|')
        if len(parts) >= 2:
            return parts[1]
    return s

def load_uniprot_fasta(fasta_path: str) -> dict:
    seqs = {}
    def _open(path):
        return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")
    with _open(fasta_path) as f:
        ac = None
        buf = []
        for line in f:
            if line.startswith(">"):
                if ac and buf:
                    seqs[ac] = "".join(buf).replace("\n","").strip()
                header = line[1:].strip()
                # Extract accession from 'sp|P12345|NAME' or 'tr|A0A024QYR6|NAME' format
                raw_ac = header.split()[0]
                ac = extract_uniprot_ac(raw_ac)
                buf = []
            else:
                buf.append(line.strip())
        if ac and buf:
            seqs[ac] = "".join(buf).replace("\n","").strip()
    return seqs

def normalize_res_list(x):
    if pd.isna(x):
        return []
    s = str(x).replace(";", ",").replace(" ", ",")
    toks = [t for t in s.split(",") if t.strip() != ""]
    out = []
    for t in toks:
        t = t.strip()
        num = ""
        for ch in t:
            if ch.isdigit():
                num += ch
            else:
                break
        if num:
            out.append(int(num))
    return sorted(set(out))

def build_sifts_index(sifts_df: pd.DataFrame):
    idx = defaultdict(list)
    needed = {"PDB","CHAIN","SP_PRIMARY","PDB_BEG","PDB_END","SP_BEG","SP_END"}
    if not needed.issubset(set(sifts_df.columns)):
        raise ValueError(f"SIFTS columns missing. Need {needed}, got {list(sifts_df.columns)}")
    # Convert numeric columns (None string -> NaN)
    for col in ["RES_BEG","RES_END","PDB_BEG","PDB_END","SP_BEG","SP_END"]:
        if col in sifts_df.columns:
            sifts_df[col] = pd.to_numeric(sifts_df[col], errors="coerce")
    # Fallback: use RES_BEG/RES_END if PDB_BEG/PDB_END are missing
    if "RES_BEG" in sifts_df.columns:
        sifts_df["PDB_BEG"] = sifts_df["PDB_BEG"].fillna(sifts_df["RES_BEG"])
    if "RES_END" in sifts_df.columns:
        sifts_df["PDB_END"] = sifts_df["PDB_END"].fillna(sifts_df["RES_END"])
    # Only drop rows missing SP_BEG/SP_END (required for UniProt mapping)
    sifts_df = sifts_df.dropna(subset=["PDB","CHAIN","SP_PRIMARY","SP_BEG","SP_END","PDB_BEG","PDB_END"])
    sifts_df[["PDB_BEG","PDB_END","SP_BEG","SP_END"]] = sifts_df[["PDB_BEG","PDB_END","SP_BEG","SP_END"]].astype(int)
    for _, r in sifts_df.iterrows():
        key = (str(r["PDB"]).lower(), str(r["CHAIN"]), str(r["SP_PRIMARY"]))
        idx[key].append((r["PDB_BEG"], r["PDB_END"], r["SP_BEG"], r["SP_END"]))
    for k in idx:
        idx[k].sort(key=lambda t: t[0])
    return idx

def pdb_to_sp_pos(pdb_pos: int, ranges):
    for PDB_BEG, PDB_END, SP_BEG, SP_END in ranges:
        if PDB_BEG <= pdb_pos <= PDB_END:
            sp = SP_BEG + (pdb_pos - PDB_BEG)
            if sp <= SP_END:
                return sp
    return None

def aa_stats(aa_list):
    n = len(aa_list)
    if n == 0:
        return {"n": 0, "hydrophobic_frac": 0.0, "polar_frac": 0.0,
                "charged_frac": 0.0, "aromatic_frac": 0.0, "net_charge": 0}
    hyd = sum(a in AA_SETS["hydrophobic"] for a in aa_list)
    pol = sum(a in AA_SETS["polar"] for a in aa_list)
    pos = sum(a in AA_SETS["charged_pos"] for a in aa_list)
    neg = sum(a in AA_SETS["charged_neg"] for a in aa_list)
    aro = sum(a in AA_SETS["aromatic"] for a in aa_list)
    return {
        "n": n,
        "hydrophobic_frac": hyd / n,
        "polar_frac": pol / n,
        "charged_frac": (pos + neg) / n,
        "aromatic_frac": aro / n,
        "net_charge": pos - neg
    }

def pick_first_present(cols, df):
    for c in cols:
        if c in df.columns:
            return c
    return None

def dedupe_upstream(df_u, merge_keys):
    # 代表値の優先序（n_atom_contacts desc → bsa_total desc → min_atom_distance asc）
    sort_cols = []
    if "n_atom_contacts" in df_u.columns:
        sort_cols.append(("n_atom_contacts", False))
    if "bsa_total" in df_u.columns:
        sort_cols.append(("bsa_total", False))
    if "min_atom_distance" in df_u.columns:
        sort_cols.append(("min_atom_distance", True))
    if sort_cols:
        by = [c for c,_ in sort_cols]
        asc = [a for _,a in sort_cols]
        df_u = df_u.sort_values(by=by, ascending=asc)
    return df_u.drop_duplicates(subset=merge_keys, keep="first")

def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    inp = Path(args.input)
    logging.info(f"Reading input: {inp}")
    df = pd.read_csv(inp, sep="\t")

    if args.limit:
        df = df.head(args.limit)

    # 決定するマージキー
    merge_keys = list(NEEDED_BASE)
    if "assembly_id" in df.columns:
        merge_keys.append("assembly_id")

    # 上流から不足列を補完
    missing_now = set()
    # 先に、どのiface列名があるか確認（入力側）
    iface_a_col = pick_first_present(CANDIDATE_IFACE_A, df)
    iface_b_col = pick_first_present(CANDIDATE_IFACE_B, df)

    needed_cols = set(merge_keys) | {"bsa_total","min_atom_distance","n_atom_contacts"} | set(CANDIDATE_IFACE_A) | set(CANDIDATE_IFACE_B)
    for up in args.upstream:
        up_path = Path(up)
        if not up_path.exists():
            logging.warning(f"[WARN] upstream not found: {up}")
            continue
        u = pd.read_csv(up_path, sep="\t")
        # マージキーが揃っていない上流は飛ばす
        u_keys = list(NEEDED_BASE)
        if "assembly_id" in df.columns and "assembly_id" in u.columns:
            u_keys.append("assembly_id")
        if not set(u_keys).issubset(set(u.columns)):
            logging.warning(f"[WARN] upstream missing merge keys {u_keys}: {up}")
            continue
        # 必要になりそうな列だけ残す
        keep_cols = [c for c in u.columns if (c in needed_cols)]
        u = u[u_keys + [c for c in keep_cols if c not in u_keys]].copy()
        # 重複キーは代表を1つに
        u = dedupe_upstream(u, u_keys)
        # left merge（不足列だけを埋める）
        before_cols = set(df.columns)
        df = df.merge(u, how="left", on=u_keys, suffixes=("", "_u"))
        # 欲しい列を確定（既にあるならそれを使い、無ければ *_u を採用）
        for col in ["bsa_total","min_atom_distance","n_atom_contacts"] + CANDIDATE_IFACE_A + CANDIDATE_IFACE_B:
            if col not in df.columns and f"{col}_u" in df.columns:
                df[col] = df[f"{col}_u"]
            elif col in df.columns and f"{col}_u" in df.columns:
                # 既存がNaNで upstream に値がある場合は埋める
                df[col] = df[col].where(df[col].notna(), df[f"{col}_u"])
        # *_u を掃除
        drop_cols = [c for c in df.columns if c.endswith("_u")]
        if drop_cols:
            df = df.drop(columns=drop_cols)

        # 更新された iface 列名の再解決
        iface_a_col = pick_first_present(CANDIDATE_IFACE_A, df)
        iface_b_col = pick_first_present(CANDIDATE_IFACE_B, df)

    # 最終チェック
    missing_base = NEEDED_BASE - set(df.columns)
    if missing_base:
        raise SystemExit(f"[ERROR] input missing merge keys: {sorted(missing_base)}")
    still_missing = []
    for col in ["min_atom_distance","n_atom_contacts"]:
        if col not in df.columns:
            still_missing.append(col)
    if not iface_a_col or not iface_b_col:
        still_missing += ["iface_res_a/iface_sasa_res_a", "iface_res_b/iface_sasa_res_b"]
    if still_missing:
        raise SystemExit(f"[ERROR] still missing columns after upstream merge: {sorted(still_missing)}.\n"
                         f"Try adding band/bsa/tight TSVs to --upstream in this order.")

    # SIFTS / FASTA 読み込み
    logging.info(f"Reading SIFTS: {args.sifts}")
    sifts_df = pd.read_csv(args.sifts, comment='#', low_memory=False)
    sifts_idx = build_sifts_index(sifts_df)

    logging.info(f"Reading UniProt FASTA: {args.fasta}")
    up_seqs = load_uniprot_fasta(args.fasta)

    # Resolve sifts_uniprot and pident column names
    sifts_a_col = pick_first_present(CANDIDATE_SIFTS_A, df)
    sifts_b_col = pick_first_present(CANDIDATE_SIFTS_B, df)
    pident_a_col = pick_first_present(CANDIDATE_PIDENT_A, df)
    pident_b_col = pick_first_present(CANDIDATE_PIDENT_B, df)

    if sifts_a_col and pident_a_col:
        logging.info(f"Using sifts columns: {sifts_a_col}, {sifts_b_col}")
        logging.info(f"Using pident columns: {pident_a_col}, {pident_b_col}")
        logging.info(f"pident threshold for fallback: {args.pident_threshold}%")
    else:
        logging.warning("sifts_uniprot or pident columns not found; fallback mapping disabled")

    # JSONL 出力
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_rows = len(df)
    n_written = 0
    n_missing_map = 0
    n_missing_seq = 0
    n_exact_match = 0
    n_sifts_fallback = 0

    with open(out, "w") as fw:
        for i, r in df.iterrows():
            ua = str(r["uniprot_a"]); ub = str(r["uniprot_b"])
            pdb = str(r["pdb_id"]).lower()
            ca  = str(r["chain_id_a"]); cb  = str(r["chain_id_b"])

            # Get sifts_uniprot and pident values
            sifts_ua = str(r[sifts_a_col]) if sifts_a_col and pd.notna(r.get(sifts_a_col)) else None
            sifts_ub = str(r[sifts_b_col]) if sifts_b_col and pd.notna(r.get(sifts_b_col)) else None
            pident_a = float(r[pident_a_col]) if pident_a_col and pd.notna(r.get(pident_a_col)) else None
            pident_b = float(r[pident_b_col]) if pident_b_col and pd.notna(r.get(pident_b_col)) else None

            # iface 残基
            ra = normalize_res_list(r[iface_a_col]) if iface_a_col in df.columns else []
            rb = normalize_res_list(r[iface_b_col]) if iface_b_col in df.columns else []

            # Determine which UniProt AC to use for SIFTS lookup (side A)
            is_exact_a = (ua == sifts_ua) if sifts_ua else True
            mapping_used_a = "none"
            key_a = (pdb, ca, ua)
            ranges_a = sifts_idx.get(key_a, [])
            if ranges_a:
                mapping_used_a = "exact"
                n_exact_match += 1
            elif sifts_ua and pident_a is not None and pident_a >= args.pident_threshold:
                # Fallback to sifts_uniprot
                key_a_fallback = (pdb, ca, sifts_ua)
                ranges_a = sifts_idx.get(key_a_fallback, [])
                if ranges_a:
                    mapping_used_a = "sifts_fallback"
                    n_sifts_fallback += 1

            # Determine which UniProt AC to use for SIFTS lookup (side B)
            is_exact_b = (ub == sifts_ub) if sifts_ub else True
            mapping_used_b = "none"
            key_b = (pdb, cb, ub)
            ranges_b = sifts_idx.get(key_b, [])
            if ranges_b:
                mapping_used_b = "exact"
            elif sifts_ub and pident_b is not None and pident_b >= args.pident_threshold:
                # Fallback to sifts_uniprot
                key_b_fallback = (pdb, cb, sifts_ub)
                ranges_b = sifts_idx.get(key_b_fallback, [])
                if ranges_b:
                    mapping_used_b = "sifts_fallback"

            if not ranges_a or not ranges_b:
                n_missing_map += 1

            # UniProt 配列 - use sifts_uniprot for sequence lookup if fallback was used
            seq_ac_a = sifts_ua if mapping_used_a == "sifts_fallback" and sifts_ua else ua
            seq_ac_b = sifts_ub if mapping_used_b == "sifts_fallback" and sifts_ub else ub
            seq_a = up_seqs.get(seq_ac_a)
            seq_b = up_seqs.get(seq_ac_b)
            if seq_a is None or seq_b is None:
                n_missing_seq += 1

            def map_side(res_list, ranges, seq):
                residues = []
                aa_list = []
                for ppos in res_list:
                    sp = pdb_to_sp_pos(ppos, ranges) if ranges else None
                    aa = None
                    if sp and seq and 1 <= sp <= len(seq):
                        aa = seq[sp-1]
                    residues.append({"pdb_pos": ppos, "uniprot_pos": sp, "aa": aa})
                    if aa:
                        aa_list.append(aa)
                return residues, aa_list

            residues_a, aa_list_a = map_side(ra, ranges_a, seq_a)
            residues_b, aa_list_b = map_side(rb, ranges_b, seq_b)

            stats_a = aa_stats(aa_list_a)
            stats_b = aa_stats(aa_list_b)

            rec = {
                "pair_id": f"{ua}|{ub}|{pdb}|{ca}|{cb}",
                "uniprot": {"a": ua, "b": ub},
                "source": {
                    "pdb_id": pdb,
                    "chain_id_a": ca,
                    "chain_id_b": cb,
                    "assembly_id": r.get("assembly_id") if "assembly_id" in df.columns else None
                },
                "mapping_info": {
                    "sifts_uniprot_a": sifts_ua,
                    "sifts_uniprot_b": sifts_ub,
                    "pident_a": pident_a,
                    "pident_b": pident_b,
                    "is_exact_match_a": is_exact_a,
                    "is_exact_match_b": is_exact_b,
                    "mapping_used_a": mapping_used_a,
                    "mapping_used_b": mapping_used_b,
                    "pident_threshold": args.pident_threshold
                },
                "metrics": {
                    "bsa_total": float(r["bsa_total"]) if "bsa_total" in df.columns and pd.notna(r["bsa_total"]) else None,
                    "n_atom_contacts": int(r["n_atom_contacts"]) if "n_atom_contacts" in df.columns and pd.notna(r["n_atom_contacts"]) else None,
                    "min_atom_distance": float(r["min_atom_distance"]) if "min_atom_distance" in df.columns and pd.notna(r["min_atom_distance"]) else None
                },
                "interface": {
                    "residues_a": residues_a,
                    "residues_b": residues_b,
                    "stats_a": stats_a,
                    "stats_b": stats_b
                },
                "provenance": {
                    "input_file": str(inp),
                    "row_index": int(i),
                    "params": {
                        "sifts_file": str(args.sifts),
                        "uniprot_fasta": str(args.fasta),
                        "upstream": args.upstream,
                        "pident_threshold": args.pident_threshold
                    },
                    "version": "0.3.0"
                }
            }
            fw.write(json.dumps(rec, ensure_ascii=False) + "\n")
            n_written += 1
            if (i+1) % args.log_every == 0:
                logging.info(f"processed {i+1}/{n_rows}")

    logging.info(f"written: {out} rows={n_written}")
    logging.info(f"exact match mapping (side A): {n_exact_match}")
    logging.info(f"sifts fallback mapping (side A, pident >= {args.pident_threshold}%): {n_sifts_fallback}")
    logging.info(f"missing SIFTS map (any-side): {n_missing_map}")
    logging.info(f"missing UniProt sequence (any-side): {n_missing_seq}")

if __name__ == "__main__":
    main()