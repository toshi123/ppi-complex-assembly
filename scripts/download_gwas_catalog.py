#!/usr/bin/env python
"""
download_gwas_catalog.py

GWAS Catalog の association summary (TSV) を
data/raw/ 以下にダウンロードするスクリプト。
"""

import argparse
import datetime
from pathlib import Path
import sys
import textwrap
import urllib.request

# 代表的な latest リリースの URL（一例）
DEFAULT_GWAS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv.gz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download GWAS Catalog associations TSV into data/raw/",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            python scripts/download_gwas_catalog.py
            python scripts/download_gwas_catalog.py --output data/raw/gwas_catalog_associations.tsv.gz
            """
        ),
    )
    parser.add_argument(
        "--url",
        type=str,
        default=DEFAULT_GWAS_URL,
        help=f"GWAS Catalog associations TSV の URL (default: {DEFAULT_GWAS_URL})",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help=(
            "保存先ファイルパス。指定がない場合は "
            "data/raw/gwas_catalog_associations_YYYYMMDD.tsv.gz を自動生成"
        ),
    )
    return parser.parse_args()


def ensure_raw_dir() -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    raw_dir = repo_root / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    return raw_dir


def default_output_path(raw_dir: Path) -> Path:
    today = datetime.date.today().strftime("%Y%m%d")
    return raw_dir / f"gwas_catalog_associations_{today}.tsv.gz"


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
        print(f"[ERROR] Failed to download GWAS Catalog: {e}", file=sys.stderr)
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
