#!/usr/bin/env python
# scripts/add_bsa_with_freesasa.py

"""
TopK/tightファイルに BSA (Buried Surface Area) と ΔSASA 由来の界面残基集合を付与する。
改良点:
  - 進捗ログ: --progress-every（既定100）
  - FreeSASAのstderr抑制: --suppress-stderr（既定ON）
  - 並列: imap_unorderedで逐次集計、失敗はスキップして続行

依存: gemmi, freesasa, pandas
"""

import argparse, sys, os, math, gzip
from pathlib import Path
from multiprocessing import Pool, cpu_count
from contextlib import contextmanager
import pandas as pd
import gemmi, freesasa

# ---------------- CLI ----------------
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--pdb-dir", default="data/raw/pdb_structures")
    ap.add_argument("--n-proc", type=int, default=max(1, cpu_count()-1))
    ap.add_argument("--assembly", action="store_true",
                    help="生物学的会合体（最初のassembly）で評価")
    ap.add_argument("--progress-every", type=int, default=100,
                    help="この件数ごとに進捗を表示")
    ap.add_argument("--suppress-stderr", dest="suppress_stderr", action="store_true", default=True,
                    help="FreeSASAなどのstderrを/dev/nullへ（既定ON）")
    ap.add_argument("--no-suppress-stderr", dest="suppress_stderr", action="store_false",
                    help="stderrを抑制しない")
    return ap.parse_args()

# ---------- stderr suppress (C拡張対応) ----------
@contextmanager
def suppress_stderr_fd(active: bool = True):
    if not active:
        yield
        return
    try:
        stderr_fd = sys.stderr.fileno()
    except Exception:
        # 環境によってはfilenoが無い
        yield
        return
    with open(os.devnull, "w") as devnull:
        old_stderr_fd = os.dup(stderr_fd)
        try:
            os.dup2(devnull.fileno(), stderr_fd)
            yield
        finally:
            try:
                os.dup2(old_stderr_fd, stderr_fd)
            finally:
                os.close(old_stderr_fd)

# -------------- IO helpers --------------
def _read_structure_gz(path: Path) -> gemmi.Structure:
    return gemmi.read_structure(str(path))

def _write_chains_to_pdb(model: gemmi.Model, chains:set, out_pdb: Path):
    st = gemmi.Structure(); st.name = "sel"
    m = gemmi.Model("1")
    names = {c.name for c in model}
    for ch in list(chains):
        if ch in names:
            m.add_chain(model[ch])
    st.add_model(m)
    st.write_minimal_pdb(str(out_pdb))

def _find_structure_path(pdb_id: str, pdb_dir: Path) -> Path | None:
    for ext in (".bcif.gz",".cif.gz",".cif",".bcif"):
        p = pdb_dir / f"{pdb_id}{ext}"
        if p.exists(): return p
    return None

# -------------- FreeSASA --------------
def _sasa_total(pdb_path: Path, suppress: bool) -> float:
    with suppress_stderr_fd(suppress):
        s = freesasa.Structure(str(pdb_path))
        r = freesasa.calc(s)
    return float(r.totalArea())

def _residue_sasa(pdb_path: Path, suppress: bool):
    """
    residue単位のSASA辞書を返す: { "123" : area, ... } （resSeq+icode）
    """
    with suppress_stderr_fd(suppress):
        s = freesasa.Structure(str(pdb_path))
        r = freesasa.calc(s)
    res_sasa = {}
    for i in range(s.nAtoms()):
        res  = s.residueNumber(i)  # e.g. "123" or "123A"
        res_sasa[res] = res_sasa.get(res, 0.0) + r.atomArea(i)
    return res_sasa

def _resid_set_from_delta(before:dict, after:dict):
    ids = set(before) | set(after)
    iface = set()
    for k in ids:
        if before.get(k,0.0) - after.get(k,0.0) > 1e-6:
            iface.add(k)
    return iface

# -------------- Worker --------------
def _job(task):
    """
    task: (row_dict, pdb_dir:str, use_assembly:bool, suppress_stderr:bool)
    """
    row, pdb_dir_s, use_assembly, suppress = task
    pdb_dir = Path(pdb_dir_s)

    pdb_id = str(row["pdb_id"]).lower()
    ca = str(row["chain_id_a"]); cb = str(row["chain_id_b"])

    gz = _find_structure_path(pdb_id, pdb_dir)
    if gz is None:
        return None

    try:
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

        outdir = Path(".tmp_bsa"); outdir.mkdir(exist_ok=True)
        complex_pdb = outdir / f"{pdb_id}_{ca}{cb}_AB.pdb"
        a_pdb       = outdir / f"{pdb_id}_{ca}_A.pdb"
        b_pdb       = outdir / f"{pdb_id}_{cb}_B.pdb"
        _write_chains_to_pdb(model, {ca,cb}, complex_pdb)
        _write_chains_to_pdb(model, {ca}, a_pdb)
        _write_chains_to_pdb(model, {cb}, b_pdb)

        try:
            sasa_a  = _sasa_total(a_pdb, suppress)
            sasa_b  = _sasa_total(b_pdb, suppress)
            sasa_ab = _sasa_total(complex_pdb, suppress)
            bsa     = sasa_a + sasa_b - sasa_ab

            sa_a = _residue_sasa(a_pdb, suppress)
            sa_b = _residue_sasa(b_pdb, suppress)
            sab  = _residue_sasa(complex_pdb, suppress)

            iface_a = _resid_set_from_delta(sa_a, sab)
            iface_b = _resid_set_from_delta(sa_b, sab)

            return {
                "pdb_id": pdb_id, "chain_id_a": ca, "chain_id_b": cb,
                "sasa_a": sasa_a, "sasa_b": sasa_b, "sasa_ab": sasa_ab,
                "bsa_total": bsa,
                "iface_sasa_a": ",".join(sorted(iface_a)),
                "iface_sasa_b": ",".join(sorted(iface_b)),
                "n_iface_sasa_a": len(iface_a),
                "n_iface_sasa_b": len(iface_b),
            }
        finally:
            # 中間PDBは掃除（必要ならコメントアウト）
            for p in (complex_pdb, a_pdb, b_pdb):
                try: p.unlink()
                except Exception: pass

    except Exception:
        # 失敗はNoneで返す（stderrは既に抑制済み）
        return None

# -------------- Main --------------
def main():
    args = parse_args()
    pdb_dir = Path(args.pdb_dir)

    df = pd.read_csv(args.input, sep="\t")
    need = {"pdb_id","chain_id_a","chain_id_b"}
    miss = need - set(df.columns)
    if miss:
        print(f"[ERROR] missing columns: {miss}", file=sys.stderr); sys.exit(1)

    # タスク化
    tasks = []
    for row in df.itertuples(index=False):
        tasks.append((
            row._asdict(),
            str(pdb_dir),
            bool(args.assembly),
            bool(args.suppress_stderr),
        ))
    total = len(tasks)
    print(f"[INFO] tasks: {total}  n_proc={args.n_proc}  assembly={args.assembly}  suppress_stderr={args.suppress_stderr}")

    results = []
    processed = ok = fail = 0

    if args.n_proc == 1:
        it = map(_job, tasks)
    else:
        pool = Pool(processes=args.n_proc)
        it = pool.imap_unordered(_job, tasks, chunksize=32)

    try:
        for res in it:
            processed += 1
            if res is None:
                fail += 1
            else:
                ok += 1
                results.append(res)
            if processed % max(1, args.progress_every) == 0:
                print(f"[INFO] progress: {processed}/{total}  ok={ok}  fail={fail}")
    finally:
        if args.n_proc != 1:
            pool.close(); pool.join()

    res_df = pd.DataFrame(results)
    out_df = df.merge(res_df, on=["pdb_id","chain_id_a","chain_id_b"], how="left")
    out_df.to_csv(args.output, sep="\t", index=False)

    print("[INFO] done")
    print(f"[INFO] processed={processed} ok={ok} fail={fail}")
    print(f"[INFO] written: {args.output} rows={len(out_df)} with_bsa={res_df.shape[0]}")

if __name__ == "__main__":
    main()
