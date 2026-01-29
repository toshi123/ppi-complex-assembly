#!/usr/bin/env python
"""
merge_residue_annotations.py

複数のデータソースを統合し、タンパク質ごとに残基レベルのアノテーションを
まとめたJSONLを出力する。

データソース:
1. human_monomer_surface.jsonl - 単量体PDBからのsurface/buriedラベル
2. interfaces.jsonl - PPIインターフェース情報
3. IDEAL.xml.gz - 天然変性領域（ディスオーダー）情報

出力フォーマット:
{
  "uniprot_id": "P12345",
  "sequence": "MKLLV...",
  "length": 500,
  "residues": {
    "1": {
      "aa": "M",
      "surface": "surface",      // surface/buried/unknown
      "rsa": 0.35,
      "sasa": 45.2,
      "disorder": false,         // from IDEAL
      "interface_partners": []   // UniProt IDs of interaction partners
    },
    ...
  },
  "summary": {
    "n_surface": 200,
    "n_buried": 150,
    "n_unknown": 150,
    "n_disorder": 50,
    "n_interface": 30,
    "ppi_partners": ["Q67890", "P11111"]
  }
}

使い方:
python scripts/merge_residue_annotations.py \
  --surface    data/processed/human_monomer_surface.jsonl \
  --interfaces data/processed/interfaces/interfaces.jsonl \
  --ideal      data/raw/IDEAL.xml.gz \
  --output     data/processed/residue_annotations.jsonl
"""

import argparse
import gzip
import json
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge surface labels, interface info, and disorder regions into unified JSONL"
    )
    p.add_argument("--surface", required=True,
                   help="Surface/buried labels JSONL (human_monomer_surface.jsonl)")
    p.add_argument("--interfaces", required=True,
                   help="Interface JSONL (interfaces.jsonl)")
    p.add_argument("--ideal", required=True,
                   help="IDEAL disorder XML (IDEAL.xml.gz)")
    p.add_argument("--output", default="data/processed/residue_annotations.jsonl",
                   help="Output merged JSONL")
    p.add_argument("--human-only", action="store_true", default=True,
                   help="Only include human proteins from IDEAL (default: True)")
    p.add_argument("--log-every", type=int, default=5000)
    return p.parse_args()


def parse_ideal_xml(ideal_path: str, human_only: bool = True) -> dict:
    """
    IDEALのXMLをパースし、UniProt ID → ディスオーダー領域のマッピングを返す
    
    Returns:
        dict: {uniprot_id: set of (start, end) tuples for disorder regions}
    """
    logging.info(f"Parsing IDEAL XML: {ideal_path}")
    
    disorder_map = defaultdict(set)
    
    def _open(path):
        if str(path).endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, "r", encoding="utf-8")
    
    with _open(ideal_path) as f:
        # XMLを逐次パース（メモリ効率）
        context = ET.iterparse(f, events=("end",))
        
        current_uniprots = []
        current_organism = None
        current_regions = []
        
        for event, elem in context:
            if elem.tag == "uniprot":
                current_uniprots.append(elem.text)
            elif elem.tag == "source_organism":
                current_organism = elem.text
            elif elem.tag == "Region":
                order_disorder = elem.findtext("order_disorder")
                if order_disorder == "disorder":
                    start = elem.findtext("region_start")
                    end = elem.findtext("region_end")
                    if start and end:
                        current_regions.append((int(start), int(end)))
            elif elem.tag == "IDEAL_entry":
                # エントリ終了時に処理
                if human_only and current_organism != "Homo sapiens":
                    pass  # skip non-human
                else:
                    for up in current_uniprots:
                        for region in current_regions:
                            disorder_map[up].add(region)
                
                # リセット
                current_uniprots = []
                current_organism = None
                current_regions = []
                
                # メモリ解放
                elem.clear()
    
    logging.info(f"IDEAL: loaded disorder regions for {len(disorder_map)} UniProt IDs")
    return dict(disorder_map)


def build_disorder_positions(disorder_regions: set) -> set:
    """
    ディスオーダー領域（start, end）のセットから、ディスオーダー位置のセットを作成
    """
    positions = set()
    for start, end in disorder_regions:
        positions.update(range(start, end + 1))
    return positions


def load_interfaces(interfaces_path: str) -> dict:
    """
    インターフェースJSONLを読み込み、UniProt ID → インターフェース情報のマッピングを返す
    
    Returns:
        dict: {uniprot_id: {position: [partner_uniprot_ids]}}
    """
    logging.info(f"Loading interfaces: {interfaces_path}")
    
    interface_map = defaultdict(lambda: defaultdict(set))
    
    def _open(path):
        if str(path).endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, "r", encoding="utf-8")
    
    n_entries = 0
    with _open(interfaces_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            n_entries += 1
            
            ua = rec["uniprot"]["a"]
            ub = rec["uniprot"]["b"]
            
            # Side A の残基
            for res in rec["interface"].get("residues_a", []):
                sp_pos = res.get("uniprot_pos")
                if sp_pos:
                    interface_map[ua][sp_pos].add(ub)
            
            # Side B の残基
            for res in rec["interface"].get("residues_b", []):
                sp_pos = res.get("uniprot_pos")
                if sp_pos:
                    interface_map[ub][sp_pos].add(ua)
    
    logging.info(f"Interfaces: loaded {n_entries} entries, {len(interface_map)} unique proteins")
    return dict(interface_map)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    
    # 1. IDEALからディスオーダー情報を読み込み
    disorder_map = parse_ideal_xml(args.ideal, human_only=args.human_only)
    
    # 2. インターフェース情報を読み込み
    interface_map = load_interfaces(args.interfaces)
    
    # 3. Surface JSONLを読み込みながら、統合して出力
    logging.info(f"Processing surface JSONL: {args.surface}")
    
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    
    def _open_read(path):
        if str(path).endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, "r", encoding="utf-8")
    
    n_processed = 0
    n_with_disorder = 0
    n_with_interface = 0
    
    with _open_read(args.surface) as fin, open(out, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            
            rec = json.loads(line)
            uniprot_id = rec["uniprot_id"]
            sequence = rec.get("sequence", "")
            length = rec.get("length", len(sequence))
            residue_labels = rec.get("residue_labels", {})
            
            # ディスオーダー位置を取得
            disorder_regions = disorder_map.get(uniprot_id, set())
            disorder_positions = build_disorder_positions(disorder_regions)
            has_disorder = len(disorder_positions) > 0
            
            # インターフェース情報を取得
            interfaces = interface_map.get(uniprot_id, {})
            has_interface = len(interfaces) > 0
            
            # 統合した残基アノテーションを構築
            merged_residues = {}
            all_partners = set()
            
            n_surface = 0
            n_buried = 0
            n_unknown = 0
            n_disorder_res = 0
            n_interface_res = 0
            
            for pos_str, info in residue_labels.items():
                pos = int(pos_str)
                aa = info.get("aa")
                label = info.get("label", "unknown")
                rsa = info.get("rsa")
                sasa = info.get("sasa")
                
                # ディスオーダー判定
                is_disorder = pos in disorder_positions
                
                # インターフェースパートナー
                partners = list(interfaces.get(pos, []))
                all_partners.update(partners)
                
                merged_residues[pos_str] = {
                    "aa": aa,
                    "surface": label,
                    "disorder": is_disorder,
                    "interface_partners": partners
                }
                
                # RSA/SASAがあれば追加
                if rsa is not None:
                    merged_residues[pos_str]["rsa"] = rsa
                if sasa is not None:
                    merged_residues[pos_str]["sasa"] = sasa
                
                # カウント
                if label == "surface":
                    n_surface += 1
                elif label == "buried":
                    n_buried += 1
                else:
                    n_unknown += 1
                
                if is_disorder:
                    n_disorder_res += 1
                if partners:
                    n_interface_res += 1
            
            # 出力レコード
            out_rec = {
                "uniprot_id": uniprot_id,
                "sequence": sequence,
                "length": length,
                "residues": merged_residues,
                "summary": {
                    "n_surface": n_surface,
                    "n_buried": n_buried,
                    "n_unknown": n_unknown,
                    "n_disorder": n_disorder_res,
                    "n_interface": n_interface_res,
                    "ppi_partners": sorted(all_partners)
                }
            }
            
            fout.write(json.dumps(out_rec, ensure_ascii=False) + "\n")
            
            n_processed += 1
            if has_disorder:
                n_with_disorder += 1
            if has_interface:
                n_with_interface += 1
            
            if n_processed % args.log_every == 0:
                logging.info(f"Processed {n_processed} proteins...")
    
    logging.info(f"Done! Output: {out}")
    logging.info(f"Total proteins: {n_processed}")
    logging.info(f"  with disorder info: {n_with_disorder}")
    logging.info(f"  with interface info: {n_with_interface}")


if __name__ == "__main__":
    main()
