#!/usr/bin/env python
"""
download_intact.py

IntAct から相互作用データ（PSI-MI TAB; .txt または .tsv）をダウンロードして
data/raw/ 以下に保存するスクリプト。

- デフォルトでは「全 IntAct」の psimitab を落とす。
- 将来的に「ヒトだけを落とす」に変えたくなったら、このスクリプトを拡張する。
"""

import argparse
import datetime
import os
from pathlib import Path
import sys
import textwrap

import urllib.request


# IntAct の psimitab（全データ）の公開 URL
DEFAULT_INTACT_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download IntAct PSIMITAB file into data/raw/",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            # デフォルトのURLから data/raw/intact_YYYYMMDD.txt として保存
            python scripts/download_intact.py

            # 保存ファイル名を指定
            python scripts/download_intact.py --output data/raw/intact_all.txt

            # URLを自分で指定（ミラーや別バージョンを使う場合）
            python scripts/download_intact.py \\
                --url https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt
            """
        ),
    )
    parser.add_argument(
        "--url",
        type=str,
        default=DEFAULT_INTACT_URL,
        help=f"IntAct psimitab を取得する URL (default: {DEFAULT_INTACT_URL})",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help=(
            "保存先ファイルパス。指定がない場合は "
            "data/raw/intact_YYYYMMDD.txt を自動生成"
        ),
    )
    return parser.parse_args()


def ensure_raw_dir() -> Path:
    """data/raw ディレクトリがなければ作る。"""
    repo_root = Path(__file__).resolve().parents[1]  # scripts/ の一つ上
    raw_dir = repo_root / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    return raw_dir


def default_output_path(raw_dir: Path) -> Path:
    """日付付きのデフォルトファイル名を返す。"""
    today = datetime.date.today().strftime("%Y%m%d")
    return raw_dir / f"intact_{today}.txt"


def download_file(url: str, output_path: Path) -> None:
    """URL からファイルをダウンロードして output_path に保存する。"""
    print(f"[INFO] Downloading IntAct from:\n  {url}")
    print(f"[INFO] Saving to:\n  {output_path}")

    try:
        with urllib.request.urlopen(url) as response:
            # ストリーミングで書き出し（大きいファイルでもOK）
            with open(output_path, "wb") as f:
                chunk_size = 1024 * 1024  # 1 MB
                total = 0
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    total += len(chunk)
                    # 軽い進捗メッセージ（煩わしければ消してOK）
                    if total % (50 * 1024 * 1024) < chunk_size:
                        print(f"[INFO] Downloaded ~{total / (1024**2):.1f} MB ...")
    except Exception as e:
        print(f"[ERROR] Failed to download IntAct: {e}", file=sys.stderr)
        if output_path.exists():
            try:
                os.remove(output_path)
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

        # 相対パスならリポジトリ root からの相対とみなす
        if not output_path.is_absolute():
            repo_root = Path(__file__).resolve().parents[1]
            output_path = (repo_root / output_path).resolve()

        # 親ディレクトリがない場合は作っておく
        output_path.parent.mkdir(parents=True, exist_ok=True)

    download_file(args.url, output_path)


if __name__ == "__main__":
    main()
