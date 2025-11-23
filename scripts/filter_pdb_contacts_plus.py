#!/usr/bin/env python
"""
filter_pdb_contacts_plus.py

IntActのペアとPDBのマッピング情報 (例: data/interim/intact_pairs_with_pdb_templates.tsv) を読み込み、
PDBの立体構造ファイル（BinaryCIF優先）を取得し、
「物理的な接触 (Heavy atom distance <= 距離閾値)」があるペアだけを抽出する。

強化ポイント（4点）:
  1) チェーン名のauth/labelずれ対策（mmCIFから対応表を作り両面で解決）
  2) Biological assemblyの展開に対応（--use-assembly）
  3) True/Falseではなく接触メトリクスも出力（最短距離・接触原子数・接触残基数）
  4) 実務Tipsを反映（BinaryCIFを優先ダウンロード、サイズ門番、レジューム＆逐次追記）

依存:
  - pandas
  - gemmi
  - requests
  - biopython (Bio.PDB.PDBList)  ※mmCIFダウンロードのフォールバック用

使用例:
  python scripts/filter_pdb_contacts_plus.py \
    --input  data/interim/intact_pairs_with_pdb_templates_distinct_chains.tsv \
    --output data/interim/intact_pairs_with_pdb_contacts.tsv \
    --pdb-dir data/raw/pdb_structures \
    --distance-cutoff 5.0 \
    --min-contacts 5 \
    --use-assembly \
    --prefer-binarycif \
    --max-file-size 50

レジューム:
  出力と同じ場所に <output>.done を作り、PDB単位で進捗を記録。途中停止後の再開が可能。
"""

import argparse
import gzip
import logging
import math
import shutil
import sys
import warnings
from pathlib import Path
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd
import requests
import gemmi
from Bio.PDB import PDBList


# -------------------- ログ設定 --------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)
warnings.simplefilter("ignore")  # Bio.PDBの冗長警告を抑制

CURRENT_PDB_ID = None  # 割り込み時表示用


# -------------------- 引数 --------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Filter IntAct-PDB mapped pairs by actual 3D structural contacts "
                    "(with auth/label resolution, assembly support, metrics output)."
    )
    p.add_argument("--input",  type=str, required=True,
                   help="入力TSV（pdb_id, chain_id_a, chain_id_b 必須）")
    p.add_argument("--output", type=str, required=True,
                   help="接触ありの行だけを逐次追記で出力。メトリクス列を追加。")
    p.add_argument("--pdb-dir", type=str, default="data/raw/pdb_structures",
                   help="構造ファイル保存先（.bcif/.cifを.gzで保存）")
    p.add_argument("--distance-cutoff", type=float, default=5.0,
                   help="重原子間の接触距離閾値 [Å]")
    p.add_argument("--min-contacts", type=int, default=5,
                   help="接触原子ペア数の下限（小さい偶然接触の弾きに）")
    p.add_argument("--max-file-size", type=float, default=50.0,
                   help="扱う圧縮ファイルの最大サイズ [MB] 超過はスキップ")
    p.add_argument("--use-assembly", action="store_true",
                   help="Biological assemblyを展開（対応IDがあれば第一アセンブリを使用）")
    p.add_argument("--prefer-binarycif", action="store_true",
                   help="BinaryCIF(.bcif)を優先して取得/使用")
    p.add_argument("--max-rows", type=int, default=None,
                   help="テスト用: 先頭N行のみ処理")
    return p.parse_args()


# -------------------- ダウンロード/読み込み --------------------
def _gz_size_mb(path: Path) -> float:
    return path.stat().st_size / (1024 * 1024)


def _gzip_file(src: Path, dst_gz: Path):
    with open(src, "rb") as f_in, gzip.open(dst_gz, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    try:
        src.unlink()  # 元を消す
    except Exception:
        pass


def _download(url: str, out_path: Path, timeout: int = 30) -> bool:
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        out_path.write_bytes(r.content)
        return True
    except Exception as e:
        logger.debug(f"Download failed [{url}]: {e}")
        return False


def ensure_structure(
    pdb_id: str,
    pdb_dir: Path,
    max_file_size_mb: float = 50.0,
    prefer_binarycif: bool = True
) -> Tuple[Optional[gemmi.Structure], Optional[Path]]:
    """
    構造をローカルに確保しgemmiで読み込む:
      - .bcif.gz / .bcif / .cif.gz / .cif の順で探索
      - 見つからなければダウンロード（BinaryCIF優先; 失敗時にmmCIF）
      - 圧縮サイズが大きすぎる場合はスキップ
    戻り値: (Structure or None, 実際に読み込んだ圧縮ファイルパス or None)
    """
    pdb_id = pdb_id.lower()
    pdb_dir.mkdir(parents=True, exist_ok=True)

    # 探索順を構築
    exts = [".bcif", ".cif"]
    if not prefer_binarycif:
        exts = [".cif", ".bcif"]

    # 1) 既存ファイルの探索（.gz優先）
    candidates: List[Path] = []
    for ext in exts:
        candidates.append(pdb_dir / f"{pdb_id}{ext}.gz")
        candidates.append(pdb_dir / f"{pdb_id}{ext}")

    chosen: Optional[Path] = None
    for cand in candidates:
        if cand.exists():
            chosen = cand
            break

    if chosen is None:
        # 2) ダウンロードを試みる
        if prefer_binarycif:
            # BinaryCIF → 失敗したら mmCIF
            tried = []
            bcif = pdb_dir / f"{pdb_id}.bcif"
            if _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.bcif", bcif):
                chosen = bcif
            else:
                tried.append("bcif")
                cif = pdb_dir / f"{pdb_id}.cif"
                if not _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif):
                    # 最後にBio.PDBのmmCif fallback
                    tried.append("cif")
                    try:
                        pdbl = PDBList(verbose=False)
                        path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(pdb_dir), file_format="mmCif", overwrite=False)
                        if path and Path(path).exists():
                            chosen = Path(path)
                        else:
                            logger.warning(f"Download failed via requests and PDBList: {pdb_id}")
                            return None, None
                    except Exception as e:
                        logger.warning(f"PDBList fallback failed for {pdb_id}: {e}")
                        return None, None
                else:
                    chosen = cif
        else:
            # mmCIF優先
            cif = pdb_dir / f"{pdb_id}.cif"
            if _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", cif):
                chosen = cif
            else:
                bcif = pdb_dir / f"{pdb_id}.bcif"
                if _download(f"https://files.rcsb.org/download/{pdb_id.upper()}.bcif", bcif):
                    chosen = bcif
                else:
                    try:
                        pdbl = PDBList(verbose=False)
                        path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(pdb_dir), file_format="mmCif", overwrite=False)
                        if path and Path(path).exists():
                            chosen = Path(path)
                        else:
                            logger.warning(f"Download failed via requests and PDBList: {pdb_id}")
                            return None, None
                    except Exception as e:
                        logger.warning(f"PDBList fallback failed for {pdb_id}: {e}")
                        return None, None

    # 3) .gz化してサイズチェック（.gzじゃなければ.gzに）
    gz_path = chosen if chosen.suffix == ".gz" else chosen.with_suffix(chosen.suffix + ".gz")
    if not chosen.suffix == ".gz":
        try:
            _gzip_file(chosen, gz_path)
        except Exception as e:
            logger.warning(f"gzip failed for {chosen}: {e}")
            gz_path = chosen  # 失敗したら非圧縮のままにする

    try:
        if gz_path.exists() and _gz_size_mb(gz_path) > max_file_size_mb:
            logger.warning(f"Skipping {pdb_id}: file size { _gz_size_mb(gz_path):.1f} MB > {max_file_size_mb} MB")
            return None, None
    except Exception:
        pass

    # 4) gemmiで読み込み
    try:
        st = gemmi.read_structure(str(gz_path))
        st.remove_waters()
        return st, gz_path
    except Exception as e:
        logger.warning(f"Failed to read structure {gz_path}: {e}")
        return None, None


# -------------------- auth/label 対応 --------------------
def build_chain_name_resolver(cif_gz_path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """mmCIFから label_asym_id ↔ auth_asym_id のラフな対応表を作る。"""
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
    """与えられたchain_idがlabel/authのどちらでも、model内の名前に解決して返す。"""
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
    """
    use_assembly=True の場合、最初の生物学的アセンブリを展開（存在すれば）。
    無ければASU（model 0）を返す。
    """
    if use_assembly:
        try:
            asm_ids = [a.id for a in struct.assemblies]
        except Exception:
            asm_ids = []
        if asm_ids:
            try:
                asm_struct = struct.make_assembly(asm_ids[0])
                return asm_struct[0]
            except Exception as e:
                logger.debug(f"make_assembly failed, fallback to ASU: {e}")
    return struct[0]


# -------------------- 接触メトリクス --------------------
def contact_metrics(model: gemmi.Model,
                    cell: gemmi.UnitCell,
                    chain_id_a: str,
                    chain_id_b: str,
                    cutoff: float) -> Tuple[float, int, int, int]:
    names = {c.name for c in model}
    if chain_id_a not in names or chain_id_b not in names:
        return math.inf, 0, 0, 0

    ca = model[chain_id_a]
    cb = model[chain_id_b]

    # ← ここを Structure の cell に
    ns = gemmi.NeighborSearch(model, cell, cutoff).populate()
    min_d = math.inf

    atoms_a = set()
    atoms_b = set()
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

    n_atoms = max(len(atoms_a), len(atoms_b))
    n_res_a = len({f"{s.num}{s.icode}" for s, _ in atoms_a})
    n_res_b = len({f"{s.num}{s.icode}" for s, _ in atoms_b})
    return (float(min_d) if math.isfinite(min_d) else math.inf,
            int(n_atoms), int(n_res_a), int(n_res_b))

# -------------------- メイン --------------------
def main():
    global CURRENT_PDB_ID
    args = parse_args()

    input_path  = Path(args.input)
    output_path = Path(args.output)
    pdb_dir     = Path(args.pdb_dir)

    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
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

    # レジューム用ログ
    done_log_path = output_path.with_suffix(".done")
    processed_pdbs = set()
    if done_log_path.exists():
        try:
            with open(done_log_path, "r") as f:
                processed_pdbs = {line.strip() for line in f if line.strip()}
            logger.info(f"Resuming: already processed {len(processed_pdbs)} PDBs.")
        except Exception as e:
            logger.warning(f"Could not read progress log: {e}")

    # 追記ヘッダ制御
    write_header = not output_path.exists() or len(processed_pdbs) == 0

    grouped = df.groupby("pdb_id")
    total_pdbs = len(grouped)
    logger.info(f"Processing {len(df)} pairs across {total_pdbs} unique PDBs...")

    processed_count = 0

    try:
        for pdb_id, group in grouped:
            processed_count += 1
            CURRENT_PDB_ID = str(pdb_id)

            if str(pdb_id) in processed_pdbs:
                if processed_count % 500 == 0:
                    logger.info(f"Skipping processed {processed_count}/{total_pdbs} ...")
                continue

            if processed_count % 100 == 0:
                logger.info(f"[{processed_count}/{total_pdbs}] PDB {pdb_id}")

            # 構造確保
            st, gz_path = ensure_structure(
                str(pdb_id),
                pdb_dir,
                max_file_size_mb=args.max_file_size,
                prefer_binarycif=args.prefer_binarycif,
            )
            if st is None or gz_path is None:
                # 読めなかった/門前で落とした → 既処理扱いにして先へ
                with open(done_log_path, "a") as f_done:
                    f_done.write(f"{pdb_id}\n")
                continue

            # auth/label対応表
            lab2auth, auth2lab = build_chain_name_resolver(gz_path)

            # assembly or ASU のモデルを取得
            model = get_model_for_contact(st, use_assembly=args.use_assembly)

            # このPDBで接触ありの行を集める
            kept_rows = []
            for idx, row in group.iterrows():
                ca_raw = str(row["chain_id_a"])
                cb_raw = str(row["chain_id_b"])

                ca = resolve_chain_name(model, ca_raw, lab2auth, auth2lab)
                cb = resolve_chain_name(model, cb_raw, lab2auth, auth2lab)
                if not ca or not cb:
                    continue

                min_d, n_atoms, n_ra, n_rb = contact_metrics(
                    model, st.cell, ca, cb, args.distance_cutoff
                )
                if (math.isfinite(min_d) and min_d <= args.distance_cutoff) and (n_atoms >= args.min_contacts):
                    out = row.to_dict()
                    out.update({
                        "min_atom_distance": float(min_d),
                        "n_atom_contacts":   int(n_atoms),
                        "n_res_contacts_a":  int(n_ra),
                        "n_res_contacts_b":  int(n_rb),
                    })
                    kept_rows.append(out)

            if kept_rows:
                pd.DataFrame(kept_rows).to_csv(
                    output_path, sep="\t", index=False, mode="a", header=write_header
                )
                write_header = False

            # 処理済みマーク
            with open(done_log_path, "a") as f_done:
                f_done.write(f"{pdb_id}\n")

    except KeyboardInterrupt:
        logger.info("\n" + "="*56)
        logger.info("Interrupted by user.")
        logger.info(f"Currently processing PDB ID: {CURRENT_PDB_ID}")
        logger.info("Re-run the same command to resume.")
        logger.info("="*56)
        sys.exit(1)

    logger.info("Filtering complete.")
    logger.info(f"Output: {output_path}")

if __name__ == "__main__":
    main()
