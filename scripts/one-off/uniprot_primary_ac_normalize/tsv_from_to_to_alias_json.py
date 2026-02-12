#!/usr/bin/env python3
"""
tsv_from_to_to_alias_json.py

UniProtのID mappingで落とした "TSV (from/to only)" を読み、
alias_ac -> primary_ac のJSON辞書を生成する。

使い方:
  python tsv_from_to_to_alias_json.py \
    --in path/to/uniprot_mapping_all.tsv \
    --out data/processed/uniprot_alias_to_primary.json

オプション:
  --keep-identity  : from==to もJSONに残す（デフォは捨てる）
  --report         : 集計レポートをstdoutに出す
"""

import argparse
import csv
import json
from collections import Counter

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="inp", required=True, help="Input TSV (from/to only)")
    p.add_argument("--out", dest="out", required=True, help="Output JSON path")
    p.add_argument("--keep-identity", action="store_true",
                   help="Keep rows where from==to (default: drop)")
    p.add_argument("--report", action="store_true",
                   help="Print summary report")
    return p.parse_args()

def norm(x: str) -> str:
    return (x or "").strip()

def main():
    args = parse_args()

    alias2primary = {}
    seen = 0
    skipped_empty = 0
    skipped_identity = 0
    conflicts = 0
    to_counts = Counter()

    with open(args.inp, "r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # UniProtのTSVはだいたい "From<TAB>To" だが、念のため柔軟に
        if header and len(header) >= 2:
            h0, h1 = header[0].strip().lower(), header[1].strip().lower()
            has_header = (h0 in {"from", "yourlist", "query"} and h1 in {"to", "uniprot"}) or (h0 == "from" and h1 == "to")
        else:
            has_header = False

        if not has_header:
            # headerが無い（or 判定できない）場合は1行目もデータとして扱う
            if header:
                row = header
                if len(row) >= 2:
                    frm, to = norm(row[0]), norm(row[1])
                    if frm and to:
                        if (not args.keep_identity) and (frm == to):
                            skipped_identity += 1
                        else:
                            if frm in alias2primary and alias2primary[frm] != to:
                                conflicts += 1
                            else:
                                alias2primary[frm] = to
                            to_counts[to] += 1
                        seen += 1
                    else:
                        skipped_empty += 1

        for row in reader:
            if len(row) < 2:
                continue
            frm, to = norm(row[0]), norm(row[1])

            if not frm or not to:
                skipped_empty += 1
                continue

            if (not args.keep_identity) and (frm == to):
                skipped_identity += 1
                continue

            # 同じfromが複数toに割れるのは基本おかしいのでカウント
            if frm in alias2primary and alias2primary[frm] != to:
                conflicts += 1
                # ここは「最初に見たものを優先」して進む（要調査）
                continue

            alias2primary[frm] = to
            to_counts[to] += 1
            seen += 1

    # JSON出力（安定のためソート）
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(dict(sorted(alias2primary.items())), f, ensure_ascii=False, indent=2)

    if args.report:
        print("== TSV -> alias2primary JSON ==")
        print(f"Input rows processed: {seen}")
        print(f"Output mappings:      {len(alias2primary)}")
        print(f"Skipped empty:        {skipped_empty}")
        print(f"Skipped identity:     {skipped_identity}  (from==to)")
        print(f"Conflicts:            {conflicts}  (same 'from' mapped to multiple 'to')")
        # 参考：to側が多い上位
        top = to_counts.most_common(10)
        if top:
            print("Top 10 'to' frequencies:")
            for k, v in top:
                print(f"  {k}\t{v}")

if __name__ == "__main__":
    main()
