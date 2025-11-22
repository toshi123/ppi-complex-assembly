#!/usr/bin/env python
"""
filter_pdb_contacts.py

IntActのペアとPDBのマッピング情報 (data/interim/intact_pairs_with_pdb_templates.tsv) を読み込み、
実際にPDBの立体構造ファイル (.cif.gz) をダウンロードして、
物理的な接触 (Heavy atom distance < 5Å) があるペアだけを抽出する。

依存ライブラリ:
- pandas
- gemmi
- biopython (Bio.PDB.PDBList)

変更点:
- 処理済みPDBのリストを管理し、途中停止後の再開（リジューム）を可能にした。
- 結果を逐次追記保存するように変更。
"""

import argparse
import gzip
import logging
import shutil
import sys
from pathlib import Path
import warnings

import gemmi
import pandas as pd
from Bio.PDB import PDBList

import time

# ログ設定
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

# 現在処理中のPDB IDを保持するグローバル変数（割り込み時の表示用）
CURRENT_PDB_ID = None

# BioPythonの警告抑制 (PDBListがverboseなため)
warnings.simplefilter("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter IntAct-PDB mapped pairs by actual 3D structural contact."
    )
    parser.add_argument(
        "--input",
        type=str,
        default="data/interim/intact_pairs_with_pdb_templates.tsv",
        help="Input TSV file path",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="data/interim/intact_pairs_with_pdb_contacts.tsv",
        help="Output TSV file path",
    )
    parser.add_argument(
        "--pdb-dir",
        type=str,
        default="data/raw/pdb_structures",
        help="Directory to store PDB structures (compressed)",
    )
    parser.add_argument(
        "--distance-cutoff",
        type=float,
        default=5.0,
        help="Distance cutoff for contact in Angstrom (Heavy atoms). Default: 5.0",
    )
    parser.add_argument(
        "--max-file-size",
        type=float,
        default=50.0,
        help="Max PDB file size (MB) to process. Larger files are skipped. Default: 50.0",
    )
    return parser.parse_args()


def get_pdb_structure(pdb_id: str, pdb_dir: Path, max_file_size_mb: float = 50.0) -> gemmi.Structure | None:
    """
    PDB構造を取得する。
    1. ローカルに {pdb_id}.cif.gz があれば読み込む。
    2. なければダウンロードして gzip 圧縮して保存する。
    3. gemmiでパースして返す。
    
    max_file_size_mb: 圧縮ファイルのサイズ上限 (MB)。これを超えるとスキップする。
    """
    pdb_id = pdb_id.lower()
    cif_gz_path = pdb_dir / f"{pdb_id}.cif.gz"
    
    # 1. 既存チェック
    if cif_gz_path.exists():
        # サイズチェック
        size_mb = cif_gz_path.stat().st_size / (1024 * 1024)
        if size_mb > max_file_size_mb:
            logger.warning(f"Skipping {pdb_id}: file size {size_mb:.1f} MB > {max_file_size_mb} MB")
            return None

        try:
            return gemmi.read_structure(str(cif_gz_path))
        except Exception as e:
            logger.warning(f"Failed to read existing {cif_gz_path}: {e}")
            return None

    # 2. ダウンロード (Bio.PDB.PDBList)
    pdbl = PDBList(verbose=False)
    try:
        # retrieve_pdb_file はダウンロードしたファイルパスを返す
        # pdir に保存される。file_format='mmCif' -> .cif
        downloaded_file = pdbl.retrieve_pdb_file(
            pdb_id, pdir=str(pdb_dir), file_format="mmCif", overwrite=False
        )
        
        if not downloaded_file or not Path(downloaded_file).exists():
            logger.warning(f"Download failed for {pdb_id}")
            return None

        downloaded_path = Path(downloaded_file)
        
        # ダウンロードされたファイルが .cif なら .cif.gz に圧縮する
        if downloaded_path.suffix == ".cif":
            with open(downloaded_path, "rb") as f_in:
                with gzip.open(cif_gz_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            # 元の .cif は削除
            downloaded_path.unlink()
            final_path = cif_gz_path
        else:
            final_path = downloaded_path

        # 圧縮後のサイズチェック
        size_mb = final_path.stat().st_size / (1024 * 1024)
        if size_mb > max_file_size_mb:
            logger.warning(f"Skipping {pdb_id}: file size {size_mb:.1f} MB > {max_file_size_mb} MB")
            # ファイルは残しておくが、ロードはしない
            return None

        return gemmi.read_structure(str(final_path))

    except Exception as e:
        logger.warning(f"Error retrieving/parsing {pdb_id}: {e}")
        return None


def check_chain_contact(
    structure: gemmi.Structure,
    chain_id_a: str,
    chain_id_b: str,
    cutoff: float
) -> bool:
    """
    structure 内の chain_id_a と chain_id_b の間に、
    cutoff (Å) 以内の距離にある原子ペア(水素除く)が存在するか判定する。
    """
    try:
        model = structure[0]
    except IndexError:
        return False

    # チェーンの存在確認
    present_chains = {c.name for c in model}
    if chain_id_a not in present_chains or chain_id_b not in present_chains:
        return False

    chain_a = model[chain_id_a]
    
    # NeighborSearch を作成 (モデル全体)
    ns = gemmi.NeighborSearch(model, structure.cell, cutoff).populate()

    # Chain A の原子をループし、近傍に Chain B の原子があるか探す
    for res in chain_a:
        for atom in res:
            # 水素は除外
            if atom.element.name == "H":
                continue
            
            # 近傍原子を検索
            marks = ns.find_atoms(atom.pos, radius=cutoff)
            
            for mark in marks:
                cra = mark.to_cra(model)
                other_chain_name = cra.chain.name
                
                if other_chain_name == chain_id_b:
                    if cra.atom.element.name == "H":
                        continue
                    return True

    return False


def main():
    global CURRENT_PDB_ID
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    pdb_dir = Path(args.pdb_dir)
    
    # 完了済みPDBリストのファイルパス (出力ファイルと同じ場所に .done をつける)
    done_log_path = output_path.with_suffix(".done")

    pdb_dir.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)

    logger.info(f"Reading input: {input_path}")
    df = pd.read_csv(input_path, sep="\t")

    required_cols = ["pdb_id", "chain_id_a", "chain_id_b"]
    if not set(required_cols).issubset(df.columns):
        logger.error(f"Missing required columns: {required_cols}")
        sys.exit(1)

    # ----- リジューム処理の準備 -----
    processed_pdbs = set()
    write_header = True

    if done_log_path.exists():
        logger.info(f"Found progress log: {done_log_path}. Resuming...")
        try:
            with open(done_log_path, "r") as f:
                processed_pdbs = {line.strip() for line in f if line.strip()}
            logger.info(f"Already processed {len(processed_pdbs)} PDBs.")
        except Exception as e:
            logger.warning(f"Could not read progress log: {e}")
    
    # 出力ファイルが既に存在する場合は、追記モードにするためヘッダーを書かない
    if output_path.exists() and len(processed_pdbs) > 0:
        write_header = False
    elif output_path.exists() and len(processed_pdbs) == 0:
        # 出力ファイルだけあってログがない場合は、安全のため上書きするか、
        # ユーザーに判断を委ねるべきだが、今回は上書きスタートとする
        logger.warning("Output file exists but no done log found. Starting over.")
        write_header = True
    else:
        write_header = True

    # ---------------------------

    grouped = df.groupby("pdb_id")
    total_pdbs = len(grouped)
    
    logger.info(f"Processing {len(df)} pairs across {total_pdbs} unique PDB structures...")

    processed_count = 0
    
    try:
        for pdb_id, group in grouped:
            processed_count += 1
            CURRENT_PDB_ID = str(pdb_id)
            
            # スキップ判定
            if str(pdb_id) in processed_pdbs:
                if processed_count % 1000 == 0:
                    logger.info(f"Skipping processed {processed_count}/{total_pdbs} ...")
                continue

            if processed_count % 100 == 0:
                logger.info(f"Processing {processed_count}/{total_pdbs} PDBs... (current: {pdb_id})")

            start_time = time.time()
            structure = get_pdb_structure(str(pdb_id), pdb_dir, args.max_file_size)
            
            # 構造が取れなかった場合も「処理済み」としてマークするか？
            # -> 「取れなかった」＝「接触なし」と同義として進めるのが安全（無限ループ防止）
            if structure is None:
                # 処理済みリストに追加
                with open(done_log_path, "a") as f_done:
                    f_done.write(f"{pdb_id}\n")
                continue

            # このPDBでの有効なペア
            valid_indices_in_pdb = []
            for idx, row in group.iterrows():
                chain_a = str(row["chain_id_a"])
                chain_b = str(row["chain_id_b"])

                if check_chain_contact(structure, chain_a, chain_b, args.distance_cutoff):
                    valid_indices_in_pdb.append(idx)
            
            # 接触があった場合のみ出力ファイルに書き込み
            if valid_indices_in_pdb:
                batch_df = df.loc[valid_indices_in_pdb].copy()
                # 追記モードで書き込み
                batch_df.to_csv(output_path, sep="\t", index=False, mode='a', header=write_header)
                write_header = False # 最初の一回書いたらあとはFalse
            
            # 処理完了ログに追記
            with open(done_log_path, "a") as f_done:
                f_done.write(f"{pdb_id}\n")
            
            elapsed = time.time() - start_time
            if elapsed > 10.0:
                logger.info(f"Slow PDB: {pdb_id} took {elapsed:.1f}s ({len(group)} pairs)")

    except KeyboardInterrupt:
        logger.info("\n" + "="*50)
        logger.info(f"Interrupted by user.")
        logger.info(f"Currently processing PDB ID: {CURRENT_PDB_ID}")
        logger.info("You can resume this process later by running the same command.")
        logger.info("="*50)
        sys.exit(1)

    logger.info(f"Filtering complete.")
    logger.info(f"Output saved to: {output_path}")


if __name__ == "__main__":
    main()
