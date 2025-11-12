#!/usr/bin/env python
"""
build_human_db.py

UniProt ヒト配列 (FASTA) から mmseqs2 の DB (human_db) を作成するスクリプト。

- 入力:
    - data/raw/uniprot_human_current.fasta
      （もしくは --input で指定したファイル）
- 出力:
    - data/interim/mmseqs/human_db   （もしくは --mmseqs-db で指定）

pdb_chain_db 側と対になる「Query側」の DB を作るイメージ。
"""

import argparse
from pathlib import Path
import sys
import textwrap
import shutil
import subprocess


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build mmseqs2 DB (human_db) from UniProt human FASTA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            # デフォルト: data/raw/uniprot_human_current.fasta -> data/interim/mmseqs/human_db
            python scripts/build_human_db.py

            # 入力FASTAとDB名を明示
            python scripts/build_human_db.py \\
                --input data/raw/uniprot_human_all_20251110.fasta \\
                --mmseqs-db data/interim/mmseqs/human_db_20251110

            # mmseqs の実行ファイル名を指定（モジュール環境など）
            python scripts/build_human_db.py --mmseqs-bin /path/to/mmseqs
            """
        ),
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help=(
            "UniProt ヒト FASTA のパス。"
            "指定がない場合は data/raw/uniprot_human_current.fasta を使う。"
        ),
    )
    parser.add_argument(
        "--mmseqs-db",
        type=str,
        default=None,
        help=(
            "作成する mmseqs DB のパス。"
            "指定がない場合は data/interim/mmseqs/human_db を使う。"
        ),
    )
    parser.add_argument(
        "--mmseqs-bin",
        type=str,
        default="mmseqs",
        help="mmseqs2 の実行ファイル名またはパス (default: mmseqs)",
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path]:
    """入力・出力パスをリポジトリルートから解決する。"""
    repo_root = Path(__file__).resolve().parents[1]

    # input FASTA
    if args.input is None:
        input_path = repo_root / "data" / "raw" / "uniprot_human_current.fasta"
    else:
        input_path = Path(args.input)
        if not input_path.is_absolute():
            input_path = (repo_root / input_path).resolve()

    # mmseqs DB prefix
    if args.mmseqs_db is None:
        db_path = repo_root / "data" / "interim" / "mmseqs" / "human_db"
    else:
        db_path = Path(args.mmseqs_db)
        if not db_path.is_absolute():
            db_path = (repo_root / db_path).resolve()
    db_path.parent.mkdir(parents=True, exist_ok=True)

    return input_path, db_path


def check_mmseqs_available(mmseqs_bin: str) -> None:
    """mmseqs2 が利用可能かどうかを簡易チェックする。"""
    mmseqs_path = shutil.which(mmseqs_bin)
    if mmseqs_path is None:
        print(
            f"[ERROR] mmseqs2 executable '{mmseqs_bin}' not found in PATH.",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"[INFO] Using mmseqs2 executable: {mmseqs_path}")


def build_mmseqs_db(mmseqs_bin: str, fasta_path: Path, db_path: Path) -> None:
    """mmseqs createdb を呼び出して DB を作成する。"""
    if not fasta_path.exists():
        print(f"[ERROR] Input FASTA not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Building mmseqs DB from FASTA: {fasta_path}")
    print(f"[INFO] mmseqs DB prefix: {db_path}")

    cmd = [
        mmseqs_bin,
        "createdb",
        str(fasta_path),
        str(db_path),
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] mmseqs createdb failed: {e}", file=sys.stderr)
        sys.exit(1)

    print("[INFO] mmseqs DB created successfully.")


def main() -> None:
    args = parse_args()
    input_path, db_path = resolve_paths(args)

    print(f"[INFO] Input FASTA : {input_path}")
    print(f"[INFO] Output DB   : {db_path}")

    check_mmseqs_available(args.mmseqs_bin)
    build_mmseqs_db(args.mmseqs_bin, input_path, db_path)


if __name__ == "__main__":
    main()