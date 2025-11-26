#!/usr/bin/env python
# scripts/add_bsa_with_freesasa.py
"""
TopK多様性ファイル（または tight/contacts ファイル）に
BSA (Buried Surface Area) と ΔSASA 由来の界面残基集合を付与する。

入力: TSV（必須列）
  - pdb_id, chain_id_a, chain_id_b
  - （任意）iface_res_a/b があってもOK（上書きはしない）

出力: 入力列 + 追加列
  - bsa_total                # Å^2
  - sasa_a, sasa_b, sasa_ab  # 参考
  - iface_sasa_a, iface_sasa_b  # ΔSASA>0 の残基ID（カンマ区切り）
  - n_iface_sasa_a, n_iface_sasa_b

使い方:
  python scripts/add_bsa_with_freesasa.py \
    --input  data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
    --output data/interim/intact_pairs_with_pdb_contacts_topk_diverse_bsa.tsv \
    --pdb-dir data/raw/pdb_structures --n-proc 6
"""

import argparse, math, sys, gzip
from pathlib import Path
from multiprocessing import Pool, cpu_count
import pandas as pd
import gemmi, freesasa

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--pdb-dir", default="data/raw/pdb_structures")
    ap.add_argument("--n-proc", type=int, default=max(1, cpu_count()-1))
    ap.add_argument("--assembly", action="store_true", help="生物学的会合体（最初のassembly）で評価")
    return ap.parse_args()

def _read_structure_gz(path: Path) -> gemmi.Structure:
    return gemmi.read_structure(str(path))

def _write_chains_to_pdb(model: gemmi.Model, chains:set, out_pdb: Path):
    st = gemmi.Structure()
    st.name = "sel"
    m = gemmi.Model("1")
    for ch in model:
        if ch.name in chains:
            m.add_chain(ch)
    st.add_model(m)
    st.write_minimal_pdb(str(out_pdb))

def _sasa_total(pdb_path: Path) -> float:
    s = freesasa.Structure(str(pdb_path))
    r = freesasa.calc(s)
    return float(r.totalArea())

def _residue_deltasasa(pdb_path: Path, tag:str):
    """
    RSAフォーマットを一度出さず、freesasa.Result から per-atom を足し上げても良いが、
    簡便のため Structure を再利用して residue 単位で集計。
    ここでは“原始的に”PDBのresSeq+icodeをIDとして返す。
    """
    s = freesasa.Structure(str(pdb_path))
    r = freesasa.calc(s)
    res_sasa = {}
    for i in range(s.nAtoms()):
        atom = s.atomName(i).strip()
        res  = s.residueNumber(i)  # 例: "123" or "123A"
        res_sasa[res] = res_sasa.get(res, 0.0) + r.atomArea(i)
    return res_sasa  # dict: { "123" : area, ... }

def _resid_set_from_delta(before:dict, after:dict):
    # ΔSASA = SASA(単独) - SASA(複合体) > 0
    ids = set(before.keys()) | set(after.keys())
    iface = set()
    for k in ids:
        if before.get(k,0.0) - after.get(k,0.0) > 1e-6:
            iface.add(k)
    return iface

def _job(row, pdb_dir:Path, use_assembly:bool):
    pdb_id = str(row.pdb_id).lower()
    ca = str(row.chain_id_a); cb = str(row.chain_id_b)
    # 構造ファイルの候補
    gz = None
    for ext in (".bcif.gz",".cif.gz",".cif",".bcif"):
        p = pdb_dir / f"{pdb_id}{ext}"
        if p.exists():
            gz = p; break
    if gz is None:
        return None

    st = _read_structure_gz(gz)
    st.remove_waters()
    model = None
    if use_assembly and len(st.assemblies)>0:
        try:
            st_asm = st.make_assembly(st.assemblies[0].id)
            model = st_asm[0]
        except Exception:
            model = st[0]
    else:
        model = st[0]

    names = {c.name for c in model}
    if ca not in names or cb not in names:
        return None

    # 一時PDBを書き出してFreeSASAで測る
    outdir = Path(".tmp_bsa"); outdir.mkdir(exist_ok=True)
    complex_pdb = outdir / f"{pdb_id}_{ca}{cb}_AB.pdb"
    a_pdb       = outdir / f"{pdb_id}_{ca}_A.pdb"
    b_pdb       = outdir / f"{pdb_id}_{cb}_B.pdb"

    _write_chains_to_pdb(model, {ca,cb}, complex_pdb)
    _write_chains_to_pdb(model, {ca}, a_pdb)
    _write_chains_to_pdb(model, {cb}, b_pdb)

    try:
        sasa_a  = _sasa_total(a_pdb)
        sasa_b  = _sasa_total(b_pdb)
        sasa_ab = _sasa_total(complex_pdb)
        bsa     = sasa_a + sasa_b - sasa_ab

        # ΔSASAで界面残基集合を作る
        sa_a = _residue_deltasasa(a_pdb, "A")
        sa_b = _residue_deltasasa(b_pdb, "B")
        sab  = _residue_deltasasa(complex_pdb, "AB")

        iface_a = _resid_set_from_delta(sa_a, sab)
        iface_b = _resid_set_from_delta(sa_b, sab)

        return {
            "pdb_id": pdb_id, "chain_id_a": ca, "chain_id_b": cb,
            "sasa_a": sasa_a, "sasa_b": sasa_b, "sasa_ab": sasa_ab,
            "bsa_total": bsa,
            "iface_sasa_a": ",".join(sorted(list(iface_a))),
            "iface_sasa_b": ",".join(sorted(list(iface_b))),
            "n_iface_sasa_a": len(iface_a),
            "n_iface_sasa_b": len(iface_b),
        }
    finally:
        # 好みで中間PDBを残したいならコメントアウト
        try:
            complex_pdb.unlink(); a_pdb.unlink(); b_pdb.unlink()
        except Exception:
            pass

def main():
    args = parse_args()
    pdb_dir = Path(args.pdb_dir)
    df = pd.read_csv(args.input, sep="\t")
    need = {"pdb_id","chain_id_a","chain_id_b"}
    miss = need - set(df.columns)
    if miss:
        print(f"[ERROR] missing columns: {miss}", file=sys.stderr); sys.exit(1)

    tasks = [(row, pdb_dir, args.assembly) for _, row in df.iterrows()]
    if args.n_proc==1:
        results = [_job(r, pdb_dir, args.assembly) for r in df.itertuples(index=False)]
    else:
        with Pool(processes=args.n_proc) as pool:
            results = pool.starmap(_job, tasks, chunksize=32)

    res_df = pd.DataFrame([r for r in results if r is not None])
    out = df.merge(res_df, on=["pdb_id","chain_id_a","chain_id_b"], how="left")
    out.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] written: {args.output} rows={len(out)}",
          f" with_bsa={res_df.shape[0]}", sep="\n")

if __name__ == "__main__":
    main()
