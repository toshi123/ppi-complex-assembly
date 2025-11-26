#!/usr/bin/env python
# scripts/add_bsa_with_freesasa.py
"""
TopK/tight などの TSV に BSA(= SASA(A)+SASA(B)-SASA(AB)) と
ΔSASA 由来の界面残基集合を付与する決定版。

改良点:
  - タスクを (pdb_id, chain_id_a, chain_id_b) でユニーク化 → 計算重複を排除
  - マージは m:1 を厳格検証 (validate="m:1") → 行数爆増を根本防止
  - ポリマー残基のみを書き出し (水/リガンド/HETATMを除外) → BSAの物理整合性向上
  - AB の原子数 ≤ A の原子数 + B の原子数 を強制チェック
  - 小さな負の BSA は丸め誤差として 0.0 にクランプ、大きな負は無効化
  - 並列 (imap_unordered) + 進捗ログ (--progress-every)
  - FreeSASA の stderr を/dev/null へ送る (--suppress-stderr; 既定ON)

依存: gemmi, freesasa, pandas
"""

import argparse
import os
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
from contextlib import contextmanager
import pandas as pd
import gemmi
import freesasa

# ---------------- CLI ----------------
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="入力 TSV (tight/TopK など)")
    ap.add_argument("--output", required=True, help="出力 TSV (BSA 列を付与)")
    ap.add_argument("--pdb-dir", default="data/raw/pdb_structures",
                    help="PDB 構造 (.cif/.bcif[.gz]) の格納ディレクトリ")
    ap.add_argument("--n-proc", type=int, default=max(1, cpu_count()-1))
    ap.add_argument("--assembly", action="store_true",
                    help="生物学的会合体(assemblies[0])で評価を試みる")
    ap.add_argument("--progress-every", type=int, default=200,
                    help="この件数ごとに進捗ログを表示")
    ap.add_argument("--suppress-stderr", dest="suppress_stderr",
                    action="store_true", default=True,
                    help="FreeSASA等のstderrを/dev/nullへ (既定ON)")
    ap.add_argument("--no-suppress-stderr", dest="suppress_stderr",
                    action="store_false", help="stderrを抑制しない")
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
        # 環境によっては fileno が無い
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
def _find_structure_path(pdb_id: str, pdb_dir: Path) -> Path | None:
    pdb_id = pdb_id.lower()
    for ext in (".bcif.gz", ".cif.gz", ".cif", ".bcif"):
        p = pdb_dir / f"{pdb_id}{ext}"
        if p.exists():
            return p
    return None

def _read_structure(path: Path) -> gemmi.Structure:
    return gemmi.read_structure(str(path))

def _make_model(struct: gemmi.Structure, use_assembly: bool) -> gemmi.Model:
    if use_assembly and len(struct.assemblies) > 0:
        try:
            st_asm = struct.make_assembly(struct.assemblies[0].id)
            return st_asm[0]
        except Exception:
            pass
    return struct[0]

def _write_polymer_chains_to_pdb(model: gemmi.Model, chains: set, out_pdb: Path):
    """
    指定チェーンのポリマー残基のみを書き出す（HETATM/水は除外）。
    """
    st = gemmi.Structure()
    st.name = "sel"
    m = gemmi.Model("1")
    names = {c.name for c in model}
    for ch in list(chains):
        if ch in names:
            new_c = gemmi.Chain(ch)
            for res in model[ch]:
                # タンパク/核酸などのポリマー残基のみ
                if res.is_polymer():
                    new_c.add_residue(res)
            if len(new_c):
                m.add_chain(new_c)
    st.add_model(m)
    st.write_minimal_pdb(str(out_pdb))

# -------------- FreeSASA wrappers --------------
def _sasa_total(pdb_path: Path, suppress_stderr: bool) -> float:
    with suppress_stderr_fd(suppress_stderr):
        s = freesasa.Structure(str(pdb_path))
        r = freesasa.calc(s)
    return float(r.totalArea())

def _residue_sasa(pdb_path: Path, suppress_stderr: bool):
    """
    residue 単位の SASA を返す: { "123" : area, ... } （resSeq+icode）
    """
    with suppress_stderr_fd(suppress_stderr):
        s = freesasa.Structure(str(pdb_path))
        r = freesasa.calc(s)
    res_sasa = {}
    for i in range(s.nAtoms()):
        res = s.residueNumber(i)  # 例: "123" or "123A"
        res_sasa[res] = res_sasa.get(res, 0.0) + r.atomArea(i)
    return res_sasa

def _resid_set_from_delta(before: dict, after: dict):
    # ΔSASA = SASA(単独) - SASA(複合体) > 0 を界面残基とする
    ids = set(before) | set(after)
    iface = set()
    for k in ids:
        if before.get(k, 0.0) - after.get(k, 0.0) > 1e-6:
            iface.add(k)
    return iface

def _atom_count(pdb_path: Path) -> int:
    with suppress_stderr_fd(True):
        s = freesasa.Structure(str(pdb_path))
    return s.nAtoms()

# -------------- Worker --------------
def _job(task):
    """
    task: (pdb_id, chain_id_a, chain_id_b, pdb_dir:str, use_assembly:bool, suppress_stderr:bool)
    戻り dict はキー (pdb_id, chain_id_a, chain_id_b) で一意。
    """
    pdb_id, ca, cb, pdb_dir_s, use_assembly, suppress = task
    pdb_dir = Path(pdb_dir_s)
    gz = _find_structure_path(pdb_id, pdb_dir)
    if gz is None:
        return None

    try:
        st = _read_structure(gz)
        st.remove_waters()
        model = _make_model(st, use_assembly)

        names = {c.name for c in model}
        if ca not in names or cb not in names:
            return None

        outdir = Path(".tmp_bsa"); outdir.mkdir(exist_ok=True)
        complex_pdb = outdir / f"{pdb_id}_{ca}{cb}_AB.pdb"
        a_pdb       = outdir / f"{pdb_id}_{ca}_A.pdb"
        b_pdb       = outdir / f"{pdb_id}_{cb}_B.pdb"

        # ポリマー限定で書き出し（A/B/ABすべて同じルール）
        _write_polymer_chains_to_pdb(model, {ca, cb}, complex_pdb)
        _write_polymer_chains_to_pdb(model, {ca}, a_pdb)
        _write_polymer_chains_to_pdb(model, {cb}, b_pdb)

        # AB に余計な原子が混ざっていないか（原子数で検証）
        atoms_a  = _atom_count(a_pdb)
        atoms_b  = _atom_count(b_pdb)
        atoms_ab = _atom_count(complex_pdb)
        if atoms_ab > atoms_a + atoms_b:
            # assembly のコピー混入や別鎖混入
            for p in (complex_pdb, a_pdb, b_pdb):
                try: p.unlink()
                except Exception: pass
            return None

        # SASA → BSA
        sasa_a  = _sasa_total(a_pdb, suppress)
        sasa_b  = _sasa_total(b_pdb, suppress)
        sasa_ab = _sasa_total(complex_pdb, suppress)
        bsa     = sasa_a + sasa_b - sasa_ab

        # 浮動小数誤差のクランプ
        if bsa < 0 and abs(bsa) < 1e-3:
            bsa = 0.0
        elif bsa < 0:
            # 明らかに不正
            for p in (complex_pdb, a_pdb, b_pdb):
                try: p.unlink()
                except Exception: pass
            return None

        # ΔSASA で界面残基集合
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

    except Exception:
        return None
    finally:
        # 中間PDBは掃除（必要ならコメントアウト）
        for p in ("AB", "A", "B"):
            try:
                Path(f".tmp_bsa/{pdb_id}_{ca}{cb}_{p}.pdb").unlink()
            except Exception:
                pass
        # 旧名（安全清掃）
        for p in (f".tmp_bsa/{pdb_id}_{ca}{cb}_AB.pdb",
                  f".tmp_bsa/{pdb_id}_{ca}_A.pdb",
                  f".tmp_bsa/{pdb_id}_{cb}_B.pdb"):
            try: Path(p).unlink()
            except Exception:
                pass

# -------------- Main --------------
def main():
    args = parse_args()
    pdb_dir = Path(args.pdb_dir)

    df = pd.read_csv(args.input, sep="\t")
    need = {"pdb_id", "chain_id_a", "chain_id_b"}
    miss = need - set(df.columns)
    if miss:
        print(f"[ERROR] missing columns: {sorted(miss)}", file=sys.stderr)
        sys.exit(1)

    # ---- タスクをユニークキーで生成（直積・重複計算を防ぐ） ----
    key_cols = ["pdb_id", "chain_id_a", "chain_id_b"]
    uniq = (df[key_cols].drop_duplicates()
            .astype({"pdb_id": str, "chain_id_a": str, "chain_id_b": str}))
    tasks = [(r.pdb_id, r.chain_id_a, r.chain_id_b,
              str(pdb_dir), bool(args.assembly), bool(args.suppress_stderr))
             for r in uniq.itertuples(index=False)]
    total = len(tasks)
    print(f"[INFO] tasks (unique keys): {total}  n_proc={args.n_proc}  assembly={args.assembly}  suppress_stderr={args.suppress_stderr}")

    # ---- 並列実行（逐次進捗） ----
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
    # 右表をキーでユニーク化（保険）
    if not res_df.empty:
        res_df = res_df.drop_duplicates(subset=key_cols)

    # ---- m:1 を厳格検証して安全にマージ ----
    out_df = df.merge(res_df, on=key_cols, how="left", validate="m:1")
    if len(out_df) != len(df):
        raise SystemExit(f"[ERROR] row count changed: input={len(df)} -> output={len(out_df)}")

    out_df.to_csv(args.output, sep="\t", index=False)
    print("[INFO] done")
    print(f"[INFO] input rows={len(df)}  output rows={len(out_df)}")
    print(f"[INFO] unique keys={total}  ok={ok}  fail={fail}")
    print(f"[INFO] written: {args.output}")

if __name__ == "__main__":
    main()
