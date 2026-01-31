#!/usr/bin/env python
"""
parse_uniprot_xml.py

UniProt XMLをパースし、タンパク質ごとの特徴情報をJSONLで出力する。

抽出する情報:
1. 残基レベル（位置指定）:
   - active_site: 活性部位
   - modified_residue: PTM修飾部位
   - disulfide_bond: ジスルフィド結合（残基ペア）
   - binding_site: 金属/イオン結合部位

2. 領域レベル（start-end）:
   - signal_peptide: シグナルペプチド
   - transmembrane: 膜貫通領域
   - topological_domain: トポロジー（cytoplasmic/extracellular など）
   - dna_binding: DNA/RNA結合領域
   - coiled_coil: コイルドコイル
   - motif: 配列モチーフ
   - domain: ドメイン

3. タンパク質レベル:
   - subcellular_locations: 詳細な細胞内局在
   - tissue_expression: 組織発現パターン（Bgeeより）

使い方:
python scripts/parse_uniprot_xml.py \
  --input data/raw/UP000005640_9606.xml.gz \
  --output data/processed/uniprot_features.jsonl
"""

import argparse
import gzip
import json
import logging
import re
from collections import defaultdict
from xml.etree import ElementTree as ET


# UniProt XML namespace
NS = {"up": "https://uniprot.org/uniprot"}


def parse_args():
    p = argparse.ArgumentParser(description="Parse UniProt XML to extract protein features")
    p.add_argument("--input", required=True, help="UniProt XML file (gzipped)")
    p.add_argument("--output", default="data/processed/uniprot_features.jsonl",
                   help="Output JSONL file")
    p.add_argument("--log-every", type=int, default=5000)
    return p.parse_args()


def open_file(path: str):
    """gzip対応でファイルを開く"""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def get_position(loc_elem):
    """<location>要素から位置情報を取得"""
    begin = loc_elem.find("up:begin", NS) if loc_elem is not None else None
    end = loc_elem.find("up:end", NS) if loc_elem is not None else None
    position = loc_elem.find("up:position", NS) if loc_elem is not None else None
    
    if position is not None:
        pos = position.get("position")
        return int(pos) if pos and pos.isdigit() else None, None
    
    start = None
    stop = None
    if begin is not None:
        pos = begin.get("position")
        start = int(pos) if pos and pos.isdigit() else None
    if end is not None:
        pos = end.get("position")
        stop = int(pos) if pos and pos.isdigit() else None
    
    return start, stop


def parse_features(entry):
    """エントリからfeature情報を抽出"""
    # 残基レベル
    active_sites = []
    modified_residues = []
    disulfide_bonds = []
    binding_sites = []
    
    # 領域レベル
    signal_peptides = []
    transmembrane_regions = []
    topological_domains = []
    dna_binding_regions = []
    coiled_coils = []
    motifs = []
    domains = []
    
    for feature in entry.findall("up:feature", NS):
        ftype = feature.get("type")
        desc = feature.get("description", "")
        loc = feature.find("up:location", NS)
        start, end = get_position(loc)
        
        if ftype == "active site":
            if start:
                active_sites.append({
                    "position": start,
                    "description": desc
                })
        
        elif ftype == "modified residue":
            if start:
                modified_residues.append({
                    "position": start,
                    "modification": desc
                })
        
        elif ftype == "disulfide bond":
            if start and end:
                disulfide_bonds.append({
                    "position1": start,
                    "position2": end
                })
            elif start:
                # 単一位置の場合（interchain など）
                disulfide_bonds.append({
                    "position1": start,
                    "position2": None,
                    "note": desc
                })
        
        elif ftype == "binding site":
            if start:
                binding_sites.append({
                    "position": start if end is None else None,
                    "start": start if end else None,
                    "end": end,
                    "ligand": desc
                })
        
        elif ftype == "signal peptide":
            if start and end:
                signal_peptides.append({
                    "start": start,
                    "end": end
                })
        
        elif ftype == "transmembrane region":
            if start and end:
                transmembrane_regions.append({
                    "start": start,
                    "end": end,
                    "description": desc
                })
        
        elif ftype == "topological domain":
            if start and end:
                # description から IN/OUT を判定
                side = None
                desc_lower = desc.lower()
                if "cytoplasmic" in desc_lower or "lumenal" not in desc_lower and "extracellular" not in desc_lower:
                    if "cytoplasmic" in desc_lower:
                        side = "IN"
                if "extracellular" in desc_lower or "lumenal" in desc_lower:
                    side = "OUT"
                
                topological_domains.append({
                    "start": start,
                    "end": end,
                    "description": desc,
                    "side": side
                })
        
        elif ftype == "DNA-binding region":
            if start and end:
                dna_binding_regions.append({
                    "start": start,
                    "end": end,
                    "description": desc
                })
        
        elif ftype == "coiled-coil region":
            if start and end:
                coiled_coils.append({
                    "start": start,
                    "end": end
                })
        
        elif ftype == "short sequence motif":
            if start and end:
                motifs.append({
                    "start": start,
                    "end": end,
                    "description": desc
                })
        
        elif ftype == "domain":
            if start and end:
                domains.append({
                    "start": start,
                    "end": end,
                    "name": desc
                })
    
    return {
        "residue_features": {
            "active_sites": active_sites,
            "modified_residues": modified_residues,
            "disulfide_bonds": disulfide_bonds,
            "binding_sites": binding_sites
        },
        "region_features": {
            "signal_peptides": signal_peptides,
            "transmembrane": transmembrane_regions,
            "topological_domains": topological_domains,
            "dna_binding": dna_binding_regions,
            "coiled_coils": coiled_coils,
            "motifs": motifs,
            "domains": domains
        }
    }


def parse_subcellular_locations(entry):
    """細胞内局在を抽出"""
    locations = []
    
    for comment in entry.findall("up:comment[@type='subcellular location']", NS):
        for subloc in comment.findall("up:subcellularLocation", NS):
            loc_chain = []
            for loc in subloc.findall("up:location", NS):
                if loc.text:
                    loc_chain.append(loc.text)
            if loc_chain:
                locations.append(loc_chain)
    
    return locations


def parse_tissue_expression(entry):
    """組織発現パターンを抽出（Bgeeより）"""
    for dbref in entry.findall("up:dbReference[@type='Bgee']", NS):
        for prop in dbref.findall("up:property[@type='expression patterns']", NS):
            value = prop.get("value", "")
            if value:
                # "Expressed in lateral globus pallidus and 104 other cell types or tissues"
                # から組織名を抽出
                match = re.match(r"Expressed in (.+?) and \d+ other", value)
                if match:
                    return {
                        "primary_tissue": match.group(1),
                        "full_description": value
                    }
                # "and N other" がない場合
                match2 = re.match(r"Expressed in (.+)$", value)
                if match2:
                    return {
                        "primary_tissue": match2.group(1),
                        "full_description": value
                    }
    return None


def parse_entry(entry):
    """1つのエントリを解析"""
    # アクセッション番号
    accessions = [acc.text for acc in entry.findall("up:accession", NS)]
    primary_ac = accessions[0] if accessions else None
    
    if not primary_ac:
        return None
    
    # 遺伝子名
    gene_names = []
    for gene in entry.findall("up:gene/up:name", NS):
        gene_names.append({
            "name": gene.text,
            "type": gene.get("type")
        })
    
    # 特徴情報
    features = parse_features(entry)
    
    # 細胞内局在
    subcellular_locations = parse_subcellular_locations(entry)
    
    # 組織発現
    tissue_expression = parse_tissue_expression(entry)
    
    return {
        "uniprot_id": primary_ac,
        "accessions": accessions,
        "gene_names": gene_names,
        **features,
        "subcellular_locations": subcellular_locations,
        "tissue_expression": tissue_expression
    }


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    
    logging.info(f"Parsing UniProt XML: {args.input}")
    
    n_entries = 0
    n_with_features = 0
    
    with open_file(args.input) as fin, open(args.output, "w") as fout:
        # iterparseで逐次処理（メモリ効率）
        context = ET.iterparse(fin, events=("end",))
        
        for event, elem in context:
            # entry要素の終了時に処理
            if elem.tag == "{https://uniprot.org/uniprot}entry":
                result = parse_entry(elem)
                
                if result:
                    # 何かしらの特徴があるかチェック
                    has_features = (
                        any(result["residue_features"].values()) or
                        any(result["region_features"].values()) or
                        result["subcellular_locations"] or
                        result["tissue_expression"]
                    )
                    
                    if has_features:
                        n_with_features += 1
                    
                    fout.write(json.dumps(result, ensure_ascii=False) + "\n")
                    n_entries += 1
                    
                    if n_entries % args.log_every == 0:
                        logging.info(f"Processed {n_entries} entries...")
                
                # メモリ解放
                elem.clear()
    
    logging.info(f"Done! Processed {n_entries} entries")
    logging.info(f"  with features: {n_with_features}")
    logging.info(f"Output: {args.output}")


if __name__ == "__main__":
    main()
