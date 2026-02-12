#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
check_primary_ac_normalization.py

目的:
1) alias->primary JSON に含まれる alias（キー）が proteins.uniprot_id に残っていないか確認
2) 主要テーブルの参照整合性（孤児行）がないか確認

実行（repoルートで）:
  python scripts/oneoff/check_primary_ac_normalization.py \
    --db data/processed/human_proteome.sqlite \
    --alias-json data/processed/uniprot_alias_to_primary.json

終了コード:
  0: OK
  1: NG（alias残存 or 孤児あり）
"""

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Iterable, List, Tuple


SQL_ORPHAN_CHECKS: List[Tuple[str, str]] = [
    ("protein_go has orphan uniprot_id",
     """
     SELECT COUNT(*) FROM protein_go g
     LEFT JOIN proteins p ON p.uniprot_id = g.uniprot_id
     WHERE p.uniprot_id IS NULL
     """),
    ("protein_locations has orphan uniprot_id",
     """
     SELECT COUNT(*) FROM protein_locations l
     LEFT JOIN proteins p ON p.uniprot_id = l.uniprot_id
     WHERE p.uniprot_id IS NULL
     """),
    ("ppi_partners has orphan uniprot_id",
     """
     SELECT COUNT(*) FROM ppi_partners x
     LEFT JOIN proteins p ON p.uniprot_id = x.uniprot_id
     WHERE p.uniprot_id IS NULL
     """),
    ("ppi_partners has orphan partner_id (partner not in proteins)",
     """
     SELECT COUNT(*) FROM ppi_partners x
     LEFT JOIN proteins p ON p.uniprot_id = x.partner_id
     WHERE p.uniprot_id IS NULL
     """),
    ("ppi_interfaces has orphan uniprot_id",
     """
     SELECT COUNT(*) FROM ppi_interfaces x
     LEFT JOIN proteins p ON p.uniprot_id = x.uniprot_id
     WHERE p.uniprot_id IS NULL
     """),
    ("ppi_interfaces has orphan partner_id (partner not in proteins)",
     """
     SELECT COUNT(*) FROM ppi_interfaces x
     LEFT JOIN proteins p ON p.uniprot_id = x.partner_id
     WHERE p.uniprot_id IS NULL
     """),
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="SQLite db path")
    p.add_argument("--alias-json", required=True, help="alias->primary JSON path")
    p.add_argument("--chunk-size", type=int, default=800,
                   help="Chunk size for IN() queries (<=999 is safe for SQLite default)")
    p.add_argument("--show-sample", type=int, default=20,
                   help="Show up to N sample remaining aliases if any")
    return p.parse_args()


def chunked(xs: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(xs), n):
        yield xs[i:i+n]


def count_in_proteins(conn: sqlite3.Connection, ac_list: List[str], chunk_size: int) -> int:
    """Count how many IDs in ac_list exist in proteins.uniprot_id (chunked)."""
    cur = conn.cursor()
    total = 0
    for chunk in chunked(ac_list, chunk_size):
        q = f"SELECT COUNT(*) FROM proteins WHERE uniprot_id IN ({','.join(['?']*len(chunk))})"
        total += cur.execute(q, chunk).fetchone()[0]
    return total


def sample_remaining_aliases(conn: sqlite3.Connection, ac_list: List[str], chunk_size: int, limit: int) -> List[str]:
    """Return a small sample of aliases still present in proteins.uniprot_id."""
    cur = conn.cursor()
    out: List[str] = []
    for chunk in chunked(ac_list, chunk_size):
        q = f"SELECT uniprot_id FROM proteins WHERE uniprot_id IN ({','.join(['?']*len(chunk))}) LIMIT {limit - len(out)}"
        rows = cur.execute(q, chunk).fetchall()
        out.extend([r[0] for r in rows])
        if len(out) >= limit:
            break
    return out


def main():
    args = parse_args()
    db_path = Path(args.db)
    alias_path = Path(args.alias_json)

    if not db_path.exists():
        raise FileNotFoundError(f"DB not found: {db_path}")
    if not alias_path.exists():
        raise FileNotFoundError(f"Alias JSON not found: {alias_path}")

    alias_map = json.loads(alias_path.read_text(encoding="utf-8"))
    aliases = [k for k, v in alias_map.items() if k and v]  # keys
    primaries = list({v for v in alias_map.values() if v})

    print("== Primary AC normalization checks ==")
    print(f"DB:         {db_path}")
    print(f"Alias JSON:  {alias_path}")
    print(f"Aliases:     {len(aliases):,}")
    print(f"Primaries:   {len(primaries):,}")
    print("")

    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()

    # 0) quick sanity: proteins count
    n_proteins = cur.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
    print(f"[Info] proteins rows: {n_proteins:,}")

    # 1) alias remaining check
    print("\n[Check 1] aliases should NOT remain in proteins.uniprot_id")
    n_alias_remaining = count_in_proteins(conn, aliases, args.chunk_size)
    print(f"  alias remaining in proteins: {n_alias_remaining:,}  (expected: 0)")

    sample = []
    if n_alias_remaining > 0 and args.show_sample > 0:
        sample = sample_remaining_aliases(conn, aliases, args.chunk_size, args.show_sample)
        if sample:
            print("  sample remaining aliases:")
            for s in sample:
                print(f"    - {s}")

    # 1b) (optional but useful) primaries existence sample
    print("\n[Info] primaries existence (sample)")
    sample_primaries = primaries[: min(500, len(primaries))]
    n_primary_found = count_in_proteins(conn, sample_primaries, args.chunk_size)
    print(f"  primaries found in proteins (sample {len(sample_primaries)}): {n_primary_found:,}")

    # 2) orphan checks
    print("\n[Check 2] orphan rows (should be 0)")
    orphan_total = 0
    orphan_details = []
    for title, sql in SQL_ORPHAN_CHECKS:
        n = cur.execute(sql).fetchone()[0]
        orphan_details.append((title, n))
        orphan_total += n

    for title, n in orphan_details:
        print(f"  {title}: {n:,}")

    conn.close()

    # verdict
    ok = (n_alias_remaining == 0) and (orphan_total == 0)
    print("\n== Verdict ==")
    if ok:
        print("OK: primary AC normalization looks consistent (no remaining aliases, no orphans).")
        raise SystemExit(0)
    else:
        print("NG: something is off.")
        if n_alias_remaining != 0:
            print("  - aliases still present in proteins (normalization incomplete or wrong alias JSON).")
        if orphan_total != 0:
            print("  - orphan rows exist (some tables reference IDs not in proteins).")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
