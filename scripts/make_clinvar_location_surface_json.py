#!/usr/bin/env python
"""
variant_location_summary.py

指定の手順に従って、以下を実行するスクリプト:

1. タンパク質の分類を取得
   - human_proteome.sqlite / human_residues.parquet から
     各 UniProt AC + 配列番号 ごとの
     - 局在 (protein_locations.location)
     - surface/buried/interface ラベル
     を取得する。
2. 変異データを取得
   - data/raw/Clinvar.json から
     - タンパク質 AC
     - 配列番号 (Position)
     - アミノ酸変化 (before -> after)
     - Clinical_significance (Disease, Polymorphism など)
     を取得する。
3. 変異を 1. の分類にマッピングして集計
4. JSON で出力
   - 「局在先aかつ表面、アミノ酸A→B、Disease: n個」のような形で取得できる構造

使い方:
    python scripts/variant_location_summary.py \
        --sqlite data/processed/human_proteome.sqlite \
        --parquet data/processed/human_residues.parquet \
        --clinvar data/raw/Clinvar.json \
        --output results/clinvar_location_surface_summary.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Any

# プロジェクトの `src` ディレクトリを import path に追加して `ppi_complex` を解決できるようにする
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from ppi_complex.db import HumanProteomeDB


def load_clinvar_variants(path: Path) -> Dict[str, Any]:
    """
    Clinvar.json を読み込み、 {uniprot_id: [variant_dict, ...]} を返す。

    期待している JSON 構造 (例):
    {
        "Q92610": [
            {
                "Variant": {"before": "Ser", "after": "Pro"},
                "Position": "98",
                "Clinical_significance": "Polymorphism",
                ...
            },
            ...
        ],
        ...
    }
    """
    with path.open("r") as f:
        data = json.load(f)
    return data


def classify_clinical_significance(label: str) -> str:
    """
    Clinical_significance を大まかに "Disease" / "Polymorphism" / "Other" に分類。
    必要に応じてここを調整してください。
    """
    if label is None:
        return "Other"

    l = label.lower()
    if "disease" in l or "pathogenic" in l or "likely pathogenic" in l:
        return "Disease"
    if "benign" in l or "likely benign" in l or "polymorphism" in l:
        return "Polymorphism"
    return "Other"


def build_variant_summary(
    db: HumanProteomeDB,
    clinvar: Dict[str, Any],
) -> Dict[str, Any]:
    """
    指定された手順 1〜4 を行い、集計結果の辞書を返す。

    返り値の構造 (例):
    {
      "<location>": {
        "surface": {
          "Ser->Pro": {"Disease": 10, "Polymorphism": 5, "Other": 1},
          ...
        },
        "buried": {
          ...
        },
        "interface": {
          ...
        }
      },
      ...
    }
    """

    # location -> surface_category -> aa_change -> significance -> count
    counts: Dict[str, Dict[str, Dict[str, Dict[str, int]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    )

    # Clinvar.json は {uniprot_id: [variants...]} 形式
    for uniprot_id, variants in clinvar.items():
        # 1-a: そのタンパク質の局在 (複数ある場合は全て使う)
        locations = db.get_protein_locations(uniprot_id)
        if not locations:
            # 局在情報が無いタンパク質はスキップ
            continue

        # 1-b: AC+position に対する surface/buried/interface 等のラベルは
        # Parquet (residues ビュー) から取得する。
        # ここでは各変異ごとに該当位置の1残基情報を取る。
        for var in variants:
            try:
                pos = int(var.get("Position"))
            except (TypeError, ValueError):
                continue

            residue = db.get_residue(uniprot_id, pos)
            if residue is None:
                # 残基情報が無ければスキップ
                continue

            aa_before = var.get("Variant", {}).get("before")
            aa_after = var.get("Variant", {}).get("after")
            if not aa_before or not aa_after:
                continue
            aa_change = f"{aa_before}->{aa_after}"

            # surface/buried/interface の情報
            # Parquet 側のカラム名を仮定: "surface" カラムに
            # 'surface' / 'buried' / 'unknown' などが入っている構造を想定。
            # interface は interface_partners != [] で判定する。
            surface_label = residue.get("surface")
            # interface_partners が空でなければ interface とみなす
            interface_partners = residue.get("interface_partners")

            # カテゴリは優先順位付きで決める:
            #   1. interface_partners が非空なら "interface"
            #   2. そうでなければ surface カラムの値を使用
            if interface_partners and interface_partners != "[]":
                region_category = "interface"
            else:
                if surface_label in ("surface", "buried"):
                    region_category = surface_label
                else:
                    region_category = "unknown"

            clin_sig_raw = var.get("Clinical_significance")
            clin_sig_cat = classify_clinical_significance(clin_sig_raw)

            # 2: 局在 × surface/buried/interface × AA変化 × Disease/Polymorphism
            for loc in locations:
                counts[loc][region_category][aa_change][clin_sig_cat] += 1

    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize Clinvar variants by localization and surface/buried/interface."
    )
    parser.add_argument(
        "--sqlite",
        type=Path,
        default=Path("data/processed/human_proteome.sqlite"),
        help="Path to human_proteome.sqlite",
    )
    parser.add_argument(
        "--parquet",
        type=Path,
        default=Path("data/processed/human_residues.parquet"),
        help="Path to human_residues.parquet",
    )
    parser.add_argument(
        "--clinvar",
        type=Path,
        default=Path("data/raw/Clinvar.json"),
        help="Path to Clinvar.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/clinvar_location_surface_summary.json"),
        help="Output JSON file path",
    )

    args = parser.parse_args()

    # 出力先ディレクトリを作成
    if args.output.parent:
        args.output.parent.mkdir(parents=True, exist_ok=True)

    clinvar_data = load_clinvar_variants(args.clinvar)

    with HumanProteomeDB(
        sqlite_path=str(args.sqlite), parquet_path=str(args.parquet)
    ) as db:
        summary = build_variant_summary(db, clinvar_data)

    # JSON で保存
    with args.output.open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print(f"Wrote summary to {args.output}")


if __name__ == "__main__":
    main()


