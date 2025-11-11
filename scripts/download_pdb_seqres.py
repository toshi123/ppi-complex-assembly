#!/usr/bin/env python
"""
download_pdb_seqres.py

PDB の SEQRES ベースの配列データ (pdb_seqres.txt) を
data/raw/ 以下にダウンロードするスクリプト。
"""

import argparse
import datetime
from pathlib import Path
import sys
import textwrap
import urllib.request

# 典型的な URL（変わる可能性はあるので適宜更新）
DEFAULT_PDB_SEQRES_URL = (
    "https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download PDB SEQRES (pdb_seqres.txt) into data/raw/",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            python scripts/download_pdb_seqres.py
            python scripts/download_pdb_seqres.py --output data/raw/pdb_seqres_latest.txt
            """
        ),
    )
    parser.add_argument(
        "--url",
        type=str,
        default=DEFAULT_PDB_SEQRES_URL,
        help=f"PDB SEQRES ファイルの URL (default: {DEFAULT_PDB_SEQRES_URL})",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="保存先ファイルパス。指定がない場合は data/raw/pdb_seqres_YYYYMMDD.txt を自動生成",
    )
    return parser.parse_args()


def ensure_raw_dir() -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    raw_dir = repo_root / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    return raw_dir


def default_output_path(raw_dir: Path) -> Path:
    today = datetime.date.today().strftime("%Y%m%d")
    return raw_dir / f"pdb_seqres_{today}.txt"


def download_file(url: str, output_path: Path) -> None:
    print(f"[INFO] Downloading from:\n  {url}")
    print(f"[INFO] Saving to:\n  {output_path}")

    try:
        with urllib.request.urlopen(url) as response:
            with open(output_path, "wb") as f:
                chunk_size = 1024 * 1024
                total = 0
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    total += len(chunk)
                    if total % (50 * 1024 * 1024) < chunk_size:
                        print(f"[INFO] Downloaded ~{total / (1024**2):.1f} MB ...")
    except Exception as e:
        print(f"[ERROR] Failed to download pdb_seqres: {e}", file=sys.stderr)
        if output_path.exists():
            try:
                output_path.unlink()
            except OSError:
                pass
        sys.exit(1)

    print("[INFO] Download completed.")


def main() -> None:
    args = parse_args()
    raw_dir = ensure_raw_dir()

    if args.output is None:
        output_path = default_output_path(raw_dir)
    else:
        output_path = Path(args.output)
        if not output_path.is_absolute():
            repo_root = Path(__file__).resolve().parents[1]
            output_path = (repo_root / output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

    download_file(args.url, output_path)


if __name__ == "__main__":
    main()
