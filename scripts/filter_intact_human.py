#!/usr/bin/env python
"""
filter_intact_human.py

IntAct の PSIMITAB ファイルから「ヒト同士の相互作用ペア（UniProt ID）」だけを抽出し、
TSV 形式で data/interim/ 以下に書き出すスクリプト。

- 入力: data/raw/intact_*.txt （または intact_current.txt）を想定
- 出力: data/interim/intact_human_pairs.tsv （デフォルト）
  - uniprot_a
  - uniprot_b
  （ペアはソートしているので A-B と B-A は同一扱い）

必要であれば後で、
- 別列に PubMed ID や detection method 等を追加する余地あり。
"""

import argparse
from pathlib import Path
import sys
import textwrap
import gzip


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter IntAct PSIMITAB to human–human UniProt pairs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples
            --------
            # デフォルト入力 (data/raw/intact_current.txt) → デフォルト出力
            python scripts/filter_intact_human.py

            # 入力・出力を明示指定
            python scripts/filter_intact_human.py \\
                --input data/raw/intact_20251106.txt \\
                --output data/interim/intact_human_pairs_20251106.tsv

            # 重複ペアを残したい場合（A-B と B-A、同じペアの多重行も残す）
            python scripts/filter_intact_human.py --keep-duplicates
            """
        ),
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help=(
            "IntAct PSIMITAB ファイルのパス。"
            "指定がない場合は data/raw/intact_current.txt を使う。"
        ),
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help=(
            "出力TSVファイルのパス。"
            "指定がない場合は data/interim/intact_human_pairs.tsv を使う。"
        ),
    )
    parser.add_argument(
        "--keep-duplicates",
        action="store_true",
        help=(
            "デフォルトでは (uniprot_a, uniprot_b) ペアをソートしてユニーク化する。"
            "このフラグを付けると、元の行ごとに出力し重複を残す。"
        ),
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path]:
    """入力・出力パスをリポジトリルートからのパスとして解決する。"""
    repo_root = Path(__file__).resolve().parents[1]

    # input
    if args.input is None:
        input_path = repo_root / "data" / "raw" / "intact_current.txt"
    else:
        input_path = Path(args.input)
        if not input_path.is_absolute():
            input_path = (repo_root / input_path).resolve()

    # output
    if args.output is None:
        output_path = repo_root / "data" / "interim" / "intact_human_pairs.tsv"
    else:
        output_path = Path(args.output)
        if not output_path.is_absolute():
            output_path = (repo_root / output_path).resolve()

    # 出力ディレクトリ作成
    output_path.parent.mkdir(parents=True, exist_ok=True)

    return input_path, output_path


def open_maybe_gzip(path: Path):
    """拡張子から gzip / 通常テキストをいい感じに判断して開く。"""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def extract_uniprot_id(id_field: str) -> str | None:
    """
    IntAct PSIMITAB の "id(s) interactor A/B" 列から UniProt ID を1つ取り出す。

    例: "uniprotkb:P12345|ensembl:ENSP0000..." → "P12345"
    UniProt エントリが見つからない場合は None を返す。
    """
    if not id_field:
        return None

    # 'uniprotkb:P12345|chebi:XXXX' のような形式を想定
    candidates = id_field.split("|")
    for token in candidates:
        token = token.strip()
        if token.lower().startswith("uniprotkb:"):
            return token.split(":", 1)[1]
    return None


def is_human_taxid(tax_field: str) -> bool:
    """
    taxidフィールドに 'taxid:9606' が含まれているかどうか。
    IntAct では "taxid:9606(Homo sapiens)" のような形式が多い。
    """
    return "taxid:9606" in tax_field


def filter_intact_human_pairs(
    input_path: Path,
    output_path: Path,
    keep_duplicates: bool = False,
) -> None:
    """
    IntAct PSIMITAB からヒト同士の UniProt ペアを抽出して TSV に出力する。
    """
    print(f"[INFO] Input : {input_path}")
    print(f"[INFO] Output: {output_path}")
    print(f"[INFO] keep_duplicates = {keep_duplicates}")

    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    total_lines = 0
    skipped_comment = 0
    skipped_nonhuman = 0
    skipped_nonuniprot = 0
    kept_pairs = 0

    pair_set: set[tuple[str, str]] = set()

    with open_maybe_gzip(input_path) as fin, open(output_path, "w", encoding="utf-8") as fout:
        # ヘッダ
        fout.write("uniprot_a\tuniprot_b\n")

        for line in fin:
            total_lines += 1
            if not line.strip():
                continue
            if line.startswith("#"):
                skipped_comment += 1
                continue

            cols = line.rstrip("\n").split("\t")
            # IntAct PSIMITAB は最低15列ある前提
            if len(cols) < 11:
                # 壊れた行はとりあえずスキップ
                continue

            id_a = cols[0]  # id(s) interactor A
            id_b = cols[1]  # id(s) interactor B
            tax_a = cols[9]  # taxid interactor A
            tax_b = cols[10]  # taxid interactor B

            # ヒト同士チェック
            if not (is_human_taxid(tax_a) and is_human_taxid(tax_b)):
                skipped_nonhuman += 1
                continue

            uniprot_a = extract_uniprot_id(id_a)
            uniprot_b = extract_uniprot_id(id_b)

            if not (uniprot_a and uniprot_b):
                skipped_nonuniprot += 1
                continue

            # 自己相互作用を含めるかどうかは後で検討。
            # ひとまず含めるが、必要ならここで uniprot_a == uniprot_b を skip してもよい。
            if keep_duplicates:
                fout.write(f"{uniprot_a}\t{uniprot_b}\n")
                kept_pairs += 1
            else:
                # undirected とみなしてソートしてユニーク化
                a, b = sorted((uniprot_a, uniprot_b))
                pair = (a, b)
                if pair not in pair_set:
                    pair_set.add(pair)

        if not keep_duplicates:
            for a, b in sorted(pair_set):
                fout.write(f"{a}\t{b}\n")
            kept_pairs = len(pair_set)

    print(f"[INFO] Total lines           : {total_lines}")
    print(f"[INFO] Skipped comments      : {skipped_comment}")
    print(f"[INFO] Skipped non-human     : {skipped_nonhuman}")
    print(f"[INFO] Skipped non-UniProt   : {skipped_nonuniprot}")
    print(f"[INFO] Kept pairs            : {kept_pairs}")
    if not keep_duplicates:
        print(f"[INFO] (unique undirected pairs)")


def main() -> None:
    args = parse_args()
    input_path, output_path = resolve_paths(args)
    filter_intact_human_pairs(
        input_path=input_path,
        output_path=output_path,
        keep_duplicates=args.keep_duplicates,
    )


if __name__ == "__main__":
    main()
