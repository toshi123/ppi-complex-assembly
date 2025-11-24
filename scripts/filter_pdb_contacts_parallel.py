#!/usr/bin/env python
"""
filter_pdb_contacts_parallel.py  (parallel + interface residue lists)

IntActのペアとPDBテンプレートTSV (例: data/interim/intact_pairs_with_pdb_templates.tsv) を読み込み、
PDB構造（BinaryCIF優先）を取得し、A/Bチェーン間に実接触がある候補のみを出力。
出力には接触メトリクスに加え、界面残基IDリスト (iface_res_a, iface_res_b) を付与。

強化ポイント:
  - auth/labelずれ解消（mmCIFから対応表を作成）
  - Biological assembly対応 (--use-assembly)
  - 接触メトリクス出力（最短距離/接触原子数/接触残基数）＋ 残基IDリスト
  - BinaryCIF優先、サイズ門番、レジューム（<output>.done）
  - PDB単位の並列処理 (--n-proc)

依存:
  - pandas, requests, numpy, gemmi, biopython (Bio.PDB.PDBList)
"""

import argparse
import gzip
import logging
import math
import shutil
import sys
import warnings
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Set

import numpy as np
import pandas as pd
import requests
import gemmi
from Bio.PDB import PDBList
from multiprocessing import Pool, cpu_count

warnings.simplefilter("ignore")  # Bio.PDBの冗長警告抑制

# -------------------- ログ --------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# -------------------- 引数 --------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Filter IntAct-PDB mapped pairs by actual 3D contacts "
                    "(assembly/auth-label aware, outputs interface residue lists, parallel)."
    )
    p.add_argument("--input",  type=str, required=True,
                   help="入力TSV（pdb_id, chain_id_a, chain_id_b 必須）")
    p.add_argument("--output", type=str, required=True,
                   help="接触ありのみ追記保存。メトリクス＋iface_res列を付与")
    p.add_argument("--pdb-dir", type=str, default="data/raw/pdb_structures",
                   help="構造ファイル保存先（.bcif/.cifを.gzで保存）")
    p.add_argument("--distance-cutoff", type=float, default=5.0,
                   help="重原子間の接触距離閾値 [Å]")
    p.add_argument("--min-contacts", type=int, default=5,
                   help="接触原子ペア数の下限")
    p.add_argument("--max-file-size", type=float, default=50.0,
                   help="扱う圧縮ファイルの最大サイズ [MB] 超過はスキップ")
    p.add_argument("--use-assembly", action="store_true",
                   help="生物学的会合体を展開（最初のassembly）")
    p.add_argument("--prefer-binarycif", action="store_true",
                   help="BinaryCIF(.bcif)を優先")
    p.add_argument("--n-proc", type=int, default=max(1, cpu_count()-1),
                   help="並列プロセス数（PDB単位）")
    p.add_argument("--max-rows", type=int, default=None,
                   help="テスト用: 先頭N行のみ処理（PDB単位で間引かれる）")
    return p.parse_args()


# -------------------- DL/IOユーティリティ --------------------
def _download(url: str, out_path: Path, timeout: int = 30) -> bool:
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        out_path.write_bytes(r.content)
        return True
    except Exception:
        return False

def _gzip_file(src: Path, dst_gz: Path):
    with open(src, "rb") as f_in, gzip.open(dst_gz, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    try:
        src.unlink()
    except Exception:
        pass

def _gz_size_mb(path: Path) -> float:
    return path.stat().st_size / (1024 * 1024)

def ensure_structure(
    pdb_id: str,
    pdb_dir: Path,
    max_file_size_mb: float = 50.0,
    prefer_binarycif: bool = True
) -> Tuple[Optional[gemmi.Structure], Optional[Path]]:
    """構造をローカルに確保しgemmiで読み込む。戻り値: (Structure, gz_path)"""
    pdb_id = pdb_id.lower()
    pdb_dir.mkdir(parents=True, exist_ok=True)

    # ローカル探索順
    exts = [".bcif", ".cif"] if prefer_binarycif else [".cif", ".bcif"]
    candidates: List[Path] = []
    for ext in exts:
        candidates += [pdb_dir / f"{pdb_id}{ext}.gz", pdb_dir / f"{pdb_id}{ext}"]

    chosen: Optional[Path] = next((c for c in candidates if c.exists()), None)
    if chosen is None:
        # ダウンロード
        if prefer_binarycif:
            bcif = pdb_dir / f"{pdb_id}.bcif"
            if not _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.bcif", bcif):
                cif = pdb_dir / f"{pdb_id}.cif"
                if not _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif):
                    # Bio.PDB fallback
                    try:
                        pdbl = PDBList(verbose=False)
                        path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(pdb_dir), file_format="mmCif", overwrite=False)
                        chosen = Path(path) if path and Path(path).exists() else None
                    except Exception:
                        chosen = None
                else:
                    chosen = cif
            else:
                chosen = bcif
        else:
            cif = pdb_dir / f"{pdb_id}.cif"
            if not _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif):
                bcif = pdb_dir / f"{pdb_id}.bcif"
                if not _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.bcif", bcif):
                    try:
                        pdbl = PDBList(verbose=False)
                        path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(pdb_dir), file_format="mmCif", overwrite=False)
                        chosen = Path(path) if path and Path(path).exists() else None
                    except Exception:
                        chosen = None
                else:
                    chosen = bcif
            else:
                chosen = cif

    if chosen is None:
        return None, None

    gz_path = chosen if chosen.suffix == ".gz" else chosen.with_suffix(chosen.suffix + ".gz")
    if chosen.suffix != ".gz":
        try:
            _gzip_file(chosen, gz_path)
        except Exception as e:
            logger.warning(f"gzip failed for {chosen}: {e}")
            gz_path = chosen  # fallback（非圧縮のまま）

    try:
        if gz_path.exists() and _gz_size_mb(gz_path) > max_file_size_mb:
            logger.warning(f"Skipping {pdb_id}: file size { _gz_size_mb(gz_path):.1f} MB > {max_file_size_mb} MB")
            return None, None
    except Exception:
        pass

    try:
        st = gemmi.read_structure(str(gz_path))
        st.remove_waters()
        return st, gz_path
    except Exception as e:
        logger.warning(f"Failed to read structure {gz_path}: {e}")
        return None, None


# -------------------- auth/label 解決 --------------------
def build_chain_name_resolver(cif_gz_path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """mmCIFから label_asym_id ↔ auth_asym_id の対応表を作る。"""
    try:
        doc = gemmi.cif.read(str(cif_gz_path))
        block = doc.sole_block()
        label = block.find_values('_atom_site.label_asym_id')
        auth  = block.find_values('_atom_site.auth_asym_id')
        lab2auth, auth2lab = {}, {}
        for la, au in zip(label, auth):
            lab2auth.setdefault(la, au)
            auth2lab.setdefault(au, la)
        return lab2auth, auth2lab
    except Exception:
        return {}, {}

def resolve_chain_name(model: gemmi.Model, chain_id: str,
                       lab2auth: Dict[str, str], auth2lab: Dict[str, str]) -> Optional[str]:
    names = {c.name for c in model}
    if chain_id in names:
        return chain_id
    if chain_id in auth2lab and auth2lab[chain_id] in names:
        return auth2lab[chain_id]
    if chain_id in lab2auth and lab2auth[chain_id] in names:
        return lab2auth[chain_id]
    return None


# -------------------- Assembly 展開 --------------------
def get_model_for_contact(struct: gemmi.Structure, use_assembly: bool) -> gemmi.Model:
    """use_assembly=Trueなら最初のBiol. Assemblyを展開（あれば）。なければASU。"""
    if use_assembly:
        try:
            asm_ids = [a.id for a in struct.assemblies]
        except Exception:
            asm_ids = []
        if asm_ids:
            try:
                asm_struct = struct.make_assembly(asm_ids[0])
                return asm_struct[0]
            except Exception:
                pass
    return struct[0]


# -------------------- 接触メトリクス＋残基集合 --------------------
def _res_uid(res: gemmi.Residue) -> str:
    """残基IDの簡易表現（auth/label混在対策まではしない）。"""
    icode = res.seqid.icode if res.seqid.icode not in (None, '') else ''
    return f"{res.seqid.num}{icode}"

def contact_metrics_and_sets(model: gemmi.Model,
                             cell: gemmi.UnitCell,
                             chain_id_a: str,
                             chain_id_b: str,
                             cutoff: float) -> Tuple[float, int, int, int, Set[str], Set[str]]:
    """
    A/Bチェーン間の:
      - 最短重原子間距離
      - 接触原子数
      - A/B側の接触残基数
      - A/B側の接触残基ID集合 (iface_res_a/b)
    接触がなければ (inf,0,0,0,∅,∅)。
    """
    names = {c.name for c in model}
    if chain_id_a not in names or chain_id_b not in names:
        return math.inf, 0, 0, 0, set(), set()

    ca = model[chain_id_a]
    ns = gemmi.NeighborSearch(model, cell, cutoff).populate()

    min_d = math.inf
    atoms_a = set()
    atoms_b = set()
    res_a: Set[str] = set()
    res_b: Set[str] = set()

    for res in ca:
        if not isinstance(res, gemmi.Residue):
            continue
        for atom in res:
            if atom.is_hydrogen():
                continue
            for mark in ns.find_atoms(atom.pos, radius=cutoff):
                cra = mark.to_cra(model)
                if cra.chain.name != chain_id_b or cra.atom.is_hydrogen():
                    continue
                d = atom.pos.dist(cra.atom.pos)
                if d < min_d:
                    min_d = d
                atoms_a.add((res.seqid, atom.name))
                atoms_b.add((cra.residue.seqid, cra.atom.name))
                res_a.add(_res_uid(res))
                res_b.add(_res_uid(cra.residue))

    n_atoms = max(len(atoms_a), len(atoms_b))
    return (float(min_d) if math.isfinite(min_d) else math.inf,
            int(n_atoms), len(res_a), len(res_b), res_a, res_b)


# -------------------- ワーカー（PDB単位） --------------------
def process_one_pdb(args_tuple) -> Tuple[str, List[Dict]]:
    """
    1つのPDBについて、接触ありの行（辞書のリスト）を返す。
    並列実行用に引数はタプルで受ける。
    """
    (pdb_id, group_df, cfg) = args_tuple
    pdb_id_str = str(pdb_id)
    kept_rows: List[Dict] = []

    st, gz_path = ensure_structure(
        pdb_id_str,
        cfg["pdb_dir"],
        max_file_size_mb=cfg["max_file_size"],
        prefer_binarycif=cfg["prefer_binarycif"],
    )
    if st is None or gz_path is None:
        return pdb_id_str, kept_rows  # 読めない→空で返す

    # auth/label対応、assemblyモデル
    lab2auth, auth2lab = build_chain_name_resolver(gz_path)
    model = get_model_for_contact(st, use_assembly=cfg["use_assembly"])

    for _, row in group_df.iterrows():
        ca_raw = str(row["chain_id_a"])
        cb_raw = str(row["chain_id_b"])
        ca = resolve_chain_name(model, ca_raw, lab2auth, auth2lab)
        cb = resolve_chain_name(model, cb_raw, lab2auth, auth2lab)
        if not ca or not cb:
            continue

        min_d, n_atoms, n_ra, n_rb, set_a, set_b = contact_metrics_and_sets(
            model, st.cell, ca, cb, cfg["distance_cutoff"]
        )
        if (math.isfinite(min_d) and min_d <= cfg["distance_cutoff"]) and (n_atoms >= cfg["min_contacts"]):
            out = row.to_dict()
            out.update({
                "min_atom_distance": float(min_d),
                "n_atom_contacts":   int(n_atoms),
                "n_res_contacts_a":  int(n_ra),
                "n_res_contacts_b":  int(n_rb),
                "iface_res_a":       ",".join(sorted(set_a, key=lambda x: (len(x), x))),
                "iface_res_b":       ",".join(sorted(set_b, key=lambda x: (len(x), x))),
            })
            kept_rows.append(out)

    return pdb_id_str, kept_rows


# -------------------- メイン --------------------
def main():
    args = parse_args()
    input_path  = Path(args.input)
    output_path = Path(args.output)
    pdb_dir     = Path(args.pdb_dir)

    if not input_path.exists():
        logger.error(f"Input not found: {input_path}")
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Reading input: {input_path}")
    df = pd.read_csv(input_path, sep="\t")
    if args.max_rows:
        df = df.head(args.max_rows).copy()

    required = {"pdb_id", "chain_id_a", "chain_id_b"}
    if not required.issubset(df.columns):
        logger.error(f"Missing required columns: {required - set(df.columns)}")
        sys.exit(1)

    # レジューム
    done_log_path = output_path.with_suffix(".done")
    processed_pdbs = set()
    if done_log_path.exists():
        try:
            with open(done_log_path, "r") as f:
                processed_pdbs = {line.strip() for line in f if line.strip()}
            logger.info(f"Resuming: already processed {len(processed_pdbs)} PDBs.")
        except Exception as e:
            logger.warning(f"Could not read progress log: {e}")

    # 書き出しヘッダ制御（メインプロセスのみが書く）
    write_header = not output_path.exists() or len(processed_pdbs) == 0

    grouped = df.groupby("pdb_id")
    items: List[Tuple[str, pd.DataFrame]] = [(pid, g) for pid, g in grouped]
    total_pdbs = len(items)
    logger.info(f"Processing {len(df)} rows across {total_pdbs} unique PDBs with {args.n_proc} proc(s)...")

    # レジューム: 既処理PDBを除外
    tasks = []
    cfg = {
        "pdb_dir": pdb_dir,
        "distance_cutoff": args.distance_cutoff,
        "min_contacts": args.min_contacts,
        "max_file_size": args.max_file_size,
        "use_assembly": args.use_assembly,
        "prefer_binarycif": args.prefer_binarycif,
    }
    for pid, g in items:
        if str(pid) in processed_pdbs:
            continue
        tasks.append((pid, g, cfg))

    # 並列実行（順序は問わない → 書き込みはメインで逐次）
    if args.n_proc == 1 or len(tasks) == 0:
        iterator = map(process_one_pdb, tasks)
    else:
        pool = Pool(processes=args.n_proc)
        iterator = pool.imap_unordered(process_one_pdb, tasks, chunksize=4)

    processed_count = 0
    try:
        for pdb_id_str, kept_rows in iterator:
            processed_count += 1
            # 追記
            if kept_rows:
                pd.DataFrame(kept_rows).to_csv(
                    output_path, sep="\t", index=False, mode="a", header=write_header
                )
                write_header = False
            # done更新
            with open(done_log_path, "a") as f_done:
                f_done.write(f"{pdb_id_str}\n")

            if processed_count % 50 == 0:
                logger.info(f"[{processed_count}/{max(1,len(tasks))}] processed (kept so far appended)")

    finally:
        if args.n_proc > 1 and 'pool' in locals():
            pool.close()
            pool.join()

    logger.info("Filtering complete.")
    logger.info(f"Output: {output_path}")

if __name__ == "__main__":
    main()
