#!/usr/bin/env python
"""
download_uniprot_human.py

UniProtKB から Homo sapiens (Human) [9606] の配列を FASTA 形式で取得し、
data/raw/ 以下に保存するスクリプト。

※ URL は UniProt の API 仕様変更で変わる可能性があるので、
   問題が出たら README / このスクリプトの URL を更新すること。
"""

import argparse
import datetime
from pathlib import Path
import sys
import textwrap
import urllib.parse
import urllib.request


# UniProtKB API (stream) を想定した URL テンプレ
# query は URL エンコードして付加する
UNIPROT_STREAM_BASE = "https://rest.uniprot.org/uniprotkb/stream"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download UniProt human proteome (FASTA) into data/raw/",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            python scripts/download_uniprot_human.py
            python scripts/download_uniprot_human.py --reviewed true
            python scripts/download_uniprot_human.py --output data/raw/uniprot_hsapiens.fasta
            """
        ),
    )
    parser.add_argument(
        "--reviewed",
        type=str,
        default="false",
        choices=["true", "false"],
        help="Swiss-Prot (reviewed) のみに限定するか (true/false, default: false)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="保存先ファイルパス。指定がない場合は data/raw/uniprot_human_YYYYMMDD.fasta を自動生成",
    )
    return parser.parse_args()


def ensure_raw_dir() -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    raw_dir = repo_root / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    return raw_dir


def default_output_path(raw_dir: Path, reviewed: bool) -> Path:
    today = datetime.date.today().strftime("%Y%m%d")
    tag = "swissprot" if reviewed else "all"
    return raw_dir / f"uniprot_human_{tag}_{today}.fasta"


def build_uniprot_url(reviewed: bool) -> str:
    query = 'organism_id:9606'
    if reviewed:
        query += ' AND reviewed:true'

    params = {
        "format": "fasta",
        "compressed": "false",
        "query": query,
    }
    return f"{UNIPROT_STREAM_BASE}?{urllib.parse.urlencode(params)}"


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
        print(f"[ERROR] Failed to download UniProt human: {e}", file=sys.stderr)
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

    reviewed = args.reviewed.lower() == "true"
    url = build_uniprot_url(reviewed)

    if args.output is None:
        output_path = default_output_path(raw_dir, reviewed)
    else:
        output_path = Path(args.output)
        if not output_path.is_absolute():
            repo_root = Path(__file__).resolve().parents[1]
            output_path = (repo_root / output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

    download_file(url, output_path)


if __name__ == "__main__":
    main()
