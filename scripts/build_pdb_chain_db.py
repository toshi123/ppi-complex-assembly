#!/usr/bin/env python
"""
build_pdb_chain_db.py

pdb_seqres ファイルから「PDBチェーンごとのFASTA」を作成し、
必要に応じて mmseqs2 の DB も作成するスクリプト。

- 入力:
    - data/raw/pdb_seqres_current.txt
      （もしくは --input で指定したファイル）
- 出力:
    - data/interim/pdb_chain_current.fasta  （もしくは --fasta-output で指定）
    - （オプション）mmseqs2 DB:
        data/interim/mmseqs/pdb_chain_db  （もしくは --mmseqs-db で指定）

まずは FASTA だけ作っておいて、
mmseqs2 を入れた後に --build-mmseqs を付けて再実行、という運用もOK。
"""

import argparse
from pathlib import Path
import sys
import textwrap
import gzip
import shutil
import subprocess


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build PDB chain FASTA (and optionally mmseqs2 DB) from pdb_seqres.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            # デフォルト: data/raw/pdb_seqres_current.txt -> data/interim/pdb_chain_current.fasta
            python scripts/build_pdb_chain_db.py

            # 入力ファイルと出力FASTAを明示
            python scripts/build_pdb_chain_db.py \\
                --input data/raw/pdb_seqres_20251111.txt \\
                --fasta-output data/interim/pdb_chain_20251111.fasta

            # FASTA を作った上で mmseqs2 DB も作成
            python scripts/build_pdb_chain_db.py --build-mmseqs

            # mmseqs2 DB のパスを明示
            python scripts/build_pdb_chain_db.py \\
                --build-mmseqs \\
                --mmseqs-db data/interim/mmseqs/pdb_chain_db
            """
        ),
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help=(
            "pdb_seqres ファイルのパス。"
            "指定がない場合は data/raw/pdb_seqres_current.txt を使う。 "
            "拡張子が .gz の場合は gzip とみなす。"
        ),
    )
    parser.add_argument(
        "--fasta-output",
        type=str,
        default=None,
        help=(
            "チェーンごとの FASTA の出力パス。"
            "指定がない場合は data/interim/pdb_chain_current.fasta を使う。"
        ),
    )
    parser.add_argument(
        "--build-mmseqs",
        action="store_true",
        help=(
            "このフラグを付けると、生成した FASTA から mmseqs2 DB "
            "(createdb) を作成する。"
        ),
    )
    parser.add_argument(
        "--mmseqs-db",
        type=str,
        default=None,
        help=(
            "--build-mmseqs が指定されたときに作成する mmseqs DB のパス。"
            "指定がない場合は data/interim/mmseqs/pdb_chain_db を使う。"
        ),
    )
    parser.add_argument(
        "--mmseqs-bin",
        type=str,
        default="mmseqs",
        help="mmseqs2 の実行ファイル名またはパス (default: mmseqs)",
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path, Path | None]:
    """入力・出力パスをリポジトリルートから解決する。"""
    repo_root = Path(__file__).resolve().parents[1]

    # input
    if args.input is None:
        input_path = repo_root / "data" / "raw" / "pdb_seqres_current.txt"
    else:
        input_path = Path(args.input)
        if not input_path.is_absolute():
            input_path = (repo_root / input_path).resolve()

    # fasta output
    if args.fasta_output is None:
        fasta_path = repo_root / "data" / "interim" / "pdb_chain_current.fasta"
    else:
        fasta_path = Path(args.fasta_output)
        if not fasta_path.is_absolute():
            fasta_path = (repo_root / fasta_path).resolve()
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    # mmseqs db path
    mmseqs_db_path: Path | None = None
    if args.build_mmseqs:
        if args.mmseqs_db is None:
            mmseqs_db_path = repo_root / "data" / "interim" / "mmseqs" / "pdb_chain_db"
        else:
            mmseqs_db_path = Path(args.mmseqs_db)
            if not mmseqs_db_path.is_absolute():
                mmseqs_db_path = (repo_root / mmseqs_db_path).resolve()
        mmseqs_db_path.parent.mkdir(parents=True, exist_ok=True)

    return input_path, fasta_path, mmseqs_db_path


def open_maybe_gzip(path: Path):
    """拡張子から gzip / 通常テキストを切り替えて開く。"""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_pdb_seqres_header(raw_header: str) -> str | None:
    """
    pdb_seqres のヘッダー行から「PDBID_Chain」の ID を作る。

    想定される形式の例:
    - '1ABC:A|PDBID|CHAIN|SEQUENCE'
    - '1abc_A|PDBID|CHAIN|SEQUENCE'
    - '1ABC:A'
    - '1ABC_A'
    など。

    戻り値:
        '1ABC_A' のような ID（PDB ID は大文字4文字）を返す。
        パースに失敗した場合は None。
    """
    header = raw_header.strip()

    if not header:
        return None

    # 最初のトークンだけ見る（後ろにコメント等がつくことがあるため）
    token = header.split()[0]

    # '1ABC:A|...' のような形式の場合、'|' の前だけ取る
    if "|" in token:
        token = token.split("|", 1)[0]

    # ここで想定される token は、
    #   '1ABC:A' もしくは '1ABC_A' など
    pdbid = None
    chain = None

    if ":" in token:
        # 例: '1ABC:A'
        pdbid, chain = token.split(":", 1)
    elif "_" in token:
        # 例: '1ABC_A'
        pdbid, chain = token.split("_", 1)
    else:
        # 例: '1ABCX' のように pdbid + chain がくっついている場合を苦し紛れに扱う
        if len(token) >= 5:
            pdbid = token[:4]
            chain = token[4:]
        else:
            pdbid = token[:4]
            chain = ""  # チェーン不明扱い

    if not pdbid:
        return None

    pdbid = pdbid[:4].upper()

    if chain is None or chain == "":
        # チェーンが取れない場合はスキップしたほうが安全
        return None

    # chain 内の区切り文字はそのまま残す（例: 'AA' など複数文字チェーンもありうる）
    return f"{pdbid}_{chain}"


def build_pdb_chain_fasta(input_path: Path, fasta_path: Path) -> None:
    """
    pdb_seqres ファイルからチェーンごとの FASTA を作成する。

    - 入力は FASTA 形式を想定。
    - 各エントリのヘッダを parse して 'PDBID_Chain' 形式に正規化する。
    """
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Reading pdb_seqres from: {input_path}")
    print(f"[INFO] Writing PDB chain FASTA to: {fasta_path}")

    n_records = 0
    n_skipped = 0

    with open_maybe_gzip(input_path) as fin, open(fasta_path, "w", encoding="utf-8") as fout:
        current_header = None
        current_seq_lines: list[str] = []

        def flush_record():
            nonlocal n_records, n_skipped
            if current_header is None:
                return
            if not current_seq_lines:
                return
            parsed_id = parse_pdb_seqres_header(current_header)
            if parsed_id is None:
                n_skipped += 1
                return
            seq = "".join(current_seq_lines).replace(" ", "").replace("\n", "")
            if not seq:
                n_skipped += 1
                return
            fout.write(f">{parsed_id}\n")
            # 60文字ごとに改行（お好みで）
            for i in range(0, len(seq), 60):
                fout.write(seq[i : i + 60] + "\n")
            n_records += 1

        for line in fin:
            if not line:
                continue
            if line.startswith(">"):
                # 直前のレコードを flush
                flush_record()
                current_header = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line.strip())

        # 最後のレコードを flush
        flush_record()

    print(f"[INFO] PDB chain FASTA records written: {n_records}")
    print(f"[INFO] Skipped records                 : {n_skipped}")


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
    print(f"[INFO] Building mmseqs DB from FASTA: {fasta_path}")
    print(f"[INFO] mmseqs DB path: {db_path}")

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
    input_path, fasta_path, mmseqs_db_path = resolve_paths(args)

    print(f"[INFO] Input pdb_seqres : {input_path}")
    print(f"[INFO] Output FASTA     : {fasta_path}")
    if args.build_mmseqs:
        print(f"[INFO] mmseqs DB        : {mmseqs_db_path}")
    else:
        print("[INFO] mmseqs DB build  : disabled (use --build-mmseqs to enable)")

    # FASTA を構築
    build_pdb_chain_fasta(input_path, fasta_path)

    # 必要なら mmseqs2 DB を構築
    if args.build_mmseqs:
        if mmseqs_db_path is None:
            print("[ERROR] mmseqs_db_path is None despite --build-mmseqs.", file=sys.stderr)
            sys.exit(1)
        check_mmseqs_available(args.mmseqs_bin)
        build_mmseqs_db(args.mmseqs_bin, fasta_path, mmseqs_db_path)


if __name__ == "__main__":
    main()