#!/usr/bin/env python
"""
convert_to_sqlite_parquet.py

JSONLおよびJSONファイルをSQLite + Parquet形式に変換する。

入力:
- residue_annotations.jsonl: 残基レベルのアノテーション
- ac2de.json: タンパク質名（description）
- ac2go.json: GO terms
- ac2seq.json: 配列（residue_annotationsにない場合の補完用）
- ac2subcellularlocations_up_classified4.json: 細胞内局在

出力:
- human_proteome.sqlite: タンパク質レベルのメタデータ
- human_residues.parquet: 残基レベルのアノテーション

使い方:
python scripts/convert_to_sqlite_parquet.py \
    --annotations data/processed/residue_annotations.jsonl \
    --ac2de /path/to/ac2json/ac2de.json \
    --ac2go /path/to/ac2json/ac2go.json \
    --ac2seq /path/to/ac2json/ac2seq.json \
    --ac2loc /path/to/ac2json/ac2subcellularlocations_up_classified4.json \
    --out-sqlite data/processed/human_proteome.sqlite \
    --out-parquet data/processed/human_residues.parquet
"""

import argparse
import gzip
import json
import logging
import sqlite3
from pathlib import Path
from collections import defaultdict

import pyarrow as pa
import pyarrow.parquet as pq


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert JSONL/JSON to SQLite + Parquet"
    )
    p.add_argument("--annotations", required=True,
                   help="residue_annotations.jsonl")
    p.add_argument("--ac2de", required=True,
                   help="ac2de.json (protein descriptions)")
    p.add_argument("--ac2go", required=True,
                   help="ac2go.json (GO terms)")
    p.add_argument("--ac2seq", required=True,
                   help="ac2seq.json (sequences, for proteins not in annotations)")
    p.add_argument("--ac2loc", required=True,
                   help="ac2subcellularlocations_up_classified4.json")
    p.add_argument("--uniprot-features",
                   help="uniprot_features.jsonl (optional, from parse_uniprot_xml.py)")
    p.add_argument("--out-sqlite", default="data/processed/human_proteome.sqlite",
                   help="Output SQLite database")
    p.add_argument("--out-parquet", default="data/processed/human_residues.parquet",
                   help="Output Parquet file")
    p.add_argument("--batch-size", type=int, default=10000,
                   help="Batch size for Parquet writes")
    p.add_argument("--log-every", type=int, default=5000)
    return p.parse_args()


def load_json(path: str) -> dict:
    """JSONファイルを読み込む"""
    logging.info(f"Loading JSON: {path}")
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    logging.info(f"  Loaded {len(data)} entries")
    return data


def open_file(path: str):
    """gzip対応でファイルを開く"""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def create_sqlite_schema(conn: sqlite3.Connection):
    """SQLiteスキーマを作成"""
    cursor = conn.cursor()
    
    # タンパク質基本情報
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS proteins (
            uniprot_id TEXT PRIMARY KEY,
            sequence TEXT,
            length INTEGER,
            description TEXT,
            n_surface INTEGER,
            n_buried INTEGER,
            n_unknown INTEGER,
            n_disorder INTEGER,
            n_interface INTEGER,
            surface_pdb_sources TEXT
        )
    """)
    
    # GO terms
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS protein_go (
            uniprot_id TEXT,
            category TEXT,
            go_term TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_protein_go_uniprot 
        ON protein_go(uniprot_id)
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_protein_go_category 
        ON protein_go(category)
    """)
    
    # 細胞内局在
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS protein_locations (
            uniprot_id TEXT,
            location TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_protein_locations_uniprot 
        ON protein_locations(uniprot_id)
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_protein_locations_location 
        ON protein_locations(location)
    """)
    
    # PPIパートナー
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ppi_partners (
            uniprot_id TEXT,
            partner_id TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_ppi_partners_uniprot 
        ON ppi_partners(uniprot_id)
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_ppi_partners_partner 
        ON ppi_partners(partner_id)
    """)
    
    # PPIインターフェイス詳細（パートナーごと）
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ppi_interfaces (
            uniprot_id TEXT,
            partner_id TEXT,
            pdb_id TEXT,
            chains TEXT,
            residues TEXT,
            n_residues INTEGER,
            bsa_total REAL,
            n_atom_contacts INTEGER,
            min_atom_distance REAL,
            stats TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id),
            PRIMARY KEY (uniprot_id, partner_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_ppi_interfaces_uniprot 
        ON ppi_interfaces(uniprot_id)
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_ppi_interfaces_partner 
        ON ppi_interfaces(partner_id)
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_ppi_interfaces_pdb 
        ON ppi_interfaces(pdb_id)
    """)
    
    # === UniProt features テーブル ===
    
    # 残基レベル特徴（活性部位、PTM、結合部位）
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS uniprot_residue_features (
            uniprot_id TEXT,
            position INTEGER,
            feature_type TEXT,
            description TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_uniprot_residue_features_uniprot 
        ON uniprot_residue_features(uniprot_id)
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_uniprot_residue_features_type 
        ON uniprot_residue_features(feature_type)
    """)
    
    # ジスルフィド結合（2残基のペア）
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS uniprot_disulfide_bonds (
            uniprot_id TEXT,
            position1 INTEGER,
            position2 INTEGER,
            note TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_uniprot_disulfide_bonds_uniprot 
        ON uniprot_disulfide_bonds(uniprot_id)
    """)
    
    # 領域レベル特徴（シグナルペプチド、膜貫通、ドメインなど）
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS uniprot_regions (
            uniprot_id TEXT,
            region_type TEXT,
            start_pos INTEGER,
            end_pos INTEGER,
            description TEXT,
            side TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_uniprot_regions_uniprot 
        ON uniprot_regions(uniprot_id)
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_uniprot_regions_type 
        ON uniprot_regions(region_type)
    """)
    
    # 詳細な細胞内局在（UniProtから）
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS uniprot_subcellular_locations (
            uniprot_id TEXT,
            location_chain TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_uniprot_subcellular_locations_uniprot 
        ON uniprot_subcellular_locations(uniprot_id)
    """)
    
    # 組織発現パターン
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS uniprot_expression (
            uniprot_id TEXT PRIMARY KEY,
            primary_tissue TEXT,
            full_description TEXT,
            FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id)
        )
    """)
    
    conn.commit()


def get_parquet_schema():
    """Parquetスキーマを定義"""
    return pa.schema([
        ("uniprot_id", pa.string()),
        ("position", pa.uint16()),
        ("aa", pa.string()),
        ("surface", pa.string()),
        ("rsa", pa.float32()),
        ("sasa", pa.float32()),
        ("disorder", pa.bool_()),
        ("interface_partners", pa.string()),  # JSON配列として格納
    ])


class ParquetBatchWriter:
    """バッチ処理でParquetに書き込むクラス"""
    
    def __init__(self, path: str, schema: pa.Schema, batch_size: int = 10000):
        self.path = path
        self.schema = schema
        self.batch_size = batch_size
        self.buffer = []
        self.writer = None
        self.total_rows = 0
    
    def _init_writer(self):
        if self.writer is None:
            self.writer = pq.ParquetWriter(
                self.path, 
                self.schema,
                compression='zstd',
                compression_level=3
            )
    
    def add_rows(self, rows: list):
        """行を追加（バッファがいっぱいになったらフラッシュ）"""
        self.buffer.extend(rows)
        if len(self.buffer) >= self.batch_size:
            self._flush()
    
    def _flush(self):
        """バッファをParquetに書き込む"""
        if not self.buffer:
            return
        
        self._init_writer()
        
        # 列ごとにデータを整理
        columns = {field.name: [] for field in self.schema}
        for row in self.buffer:
            for field in self.schema:
                columns[field.name].append(row.get(field.name))
        
        # PyArrow配列に変換
        arrays = []
        for field in self.schema:
            if field.type == pa.float32():
                # NoneをNaNとして扱う
                arr = pa.array(columns[field.name], type=pa.float32())
            elif field.type == pa.bool_():
                arr = pa.array(columns[field.name], type=pa.bool_())
            elif field.type == pa.uint16():
                arr = pa.array(columns[field.name], type=pa.uint16())
            else:
                arr = pa.array(columns[field.name], type=pa.string())
            arrays.append(arr)
        
        table = pa.Table.from_arrays(arrays, schema=self.schema)
        self.writer.write_table(table)
        
        self.total_rows += len(self.buffer)
        self.buffer = []
    
    def close(self):
        """残りのバッファをフラッシュしてクローズ"""
        self._flush()
        if self.writer:
            self.writer.close()
        return self.total_rows


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    
    # 出力ディレクトリの作成
    Path(args.out_sqlite).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_parquet).parent.mkdir(parents=True, exist_ok=True)
    
    # 既存ファイルがあれば削除
    if Path(args.out_sqlite).exists():
        Path(args.out_sqlite).unlink()
    if Path(args.out_parquet).exists():
        Path(args.out_parquet).unlink()
    
    # 補助データの読み込み
    ac2de = load_json(args.ac2de)
    ac2go = load_json(args.ac2go)
    ac2loc = load_json(args.ac2loc)
    ac2seq = load_json(args.ac2seq)
    
    # SQLite接続
    logging.info(f"Creating SQLite database: {args.out_sqlite}")
    conn = sqlite3.connect(args.out_sqlite)
    create_sqlite_schema(conn)
    cursor = conn.cursor()
    
    # Parquetライター
    logging.info(f"Creating Parquet file: {args.out_parquet}")
    parquet_writer = ParquetBatchWriter(
        args.out_parquet,
        get_parquet_schema(),
        batch_size=args.batch_size
    )
    
    # メイン処理: residue_annotations.jsonl を読み込みながら変換
    logging.info(f"Processing annotations: {args.annotations}")
    
    n_proteins = 0
    n_residues = 0
    processed_acs = set()
    
    # GO terms と locations のバッチ挿入用
    go_batch = []
    loc_batch = []
    ppi_batch = []
    ppi_iface_batch = []
    
    with open_file(args.annotations) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            rec = json.loads(line)
            uniprot_id = rec["uniprot_id"]
            sequence = rec.get("sequence", "")
            length = rec.get("length", len(sequence))
            residues = rec.get("residues", {})
            summary = rec.get("summary", {})
            surface_pdb_sources = rec.get("surface_pdb_sources", [])
            ppi_interfaces = rec.get("ppi_interfaces", {})
            
            processed_acs.add(uniprot_id)
            
            # タンパク質情報をSQLiteに挿入
            description = ac2de.get(uniprot_id, "")
            cursor.execute("""
                INSERT OR REPLACE INTO proteins 
                (uniprot_id, sequence, length, description, 
                 n_surface, n_buried, n_unknown, n_disorder, n_interface,
                 surface_pdb_sources)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                uniprot_id,
                sequence,
                length,
                description,
                summary.get("n_surface", 0),
                summary.get("n_buried", 0),
                summary.get("n_unknown", 0),
                summary.get("n_disorder", 0),
                summary.get("n_interface", 0),
                json.dumps(surface_pdb_sources) if surface_pdb_sources else "[]"
            ))
            
            # GO terms
            if uniprot_id in ac2go:
                go_data = ac2go[uniprot_id]
                for category in ["F", "C", "P"]:
                    for term in go_data.get(category, []):
                        go_batch.append((uniprot_id, category, term))
            
            # Locations
            if uniprot_id in ac2loc:
                for loc in ac2loc[uniprot_id]:
                    loc_batch.append((uniprot_id, loc))
            
            # PPI partners
            for partner in summary.get("ppi_partners", []):
                ppi_batch.append((uniprot_id, partner))
            
            # PPI interfaces（パートナーごとの詳細）
            for partner_id, iface_info in ppi_interfaces.items():
                ppi_iface_batch.append((
                    uniprot_id,
                    partner_id,
                    iface_info.get("pdb_id"),
                    json.dumps(iface_info.get("chains", [])),
                    json.dumps(iface_info.get("residues", [])),
                    iface_info.get("n_residues", 0),
                    iface_info.get("bsa_total"),
                    iface_info.get("n_atom_contacts"),
                    iface_info.get("min_atom_distance"),
                    json.dumps(iface_info.get("stats")) if iface_info.get("stats") else None
                ))
            
            # 残基データをParquetに追加
            residue_rows = []
            for pos_str, info in residues.items():
                pos = int(pos_str)
                partners = info.get("interface_partners", [])
                partners_json = json.dumps(partners) if partners else "[]"
                
                residue_rows.append({
                    "uniprot_id": uniprot_id,
                    "position": pos,
                    "aa": info.get("aa", ""),
                    "surface": info.get("surface", "unknown"),
                    "rsa": info.get("rsa"),
                    "sasa": info.get("sasa"),
                    "disorder": info.get("disorder", False),
                    "interface_partners": partners_json
                })
                n_residues += 1
            
            parquet_writer.add_rows(residue_rows)
            
            n_proteins += 1
            
            # バッチ挿入（定期的に）
            if n_proteins % 1000 == 0:
                if go_batch:
                    cursor.executemany(
                        "INSERT INTO protein_go (uniprot_id, category, go_term) VALUES (?, ?, ?)",
                        go_batch
                    )
                    go_batch = []
                if loc_batch:
                    cursor.executemany(
                        "INSERT INTO protein_locations (uniprot_id, location) VALUES (?, ?)",
                        loc_batch
                    )
                    loc_batch = []
                if ppi_batch:
                    cursor.executemany(
                        "INSERT INTO ppi_partners (uniprot_id, partner_id) VALUES (?, ?)",
                        ppi_batch
                    )
                    ppi_batch = []
                if ppi_iface_batch:
                    cursor.executemany("""
                        INSERT OR REPLACE INTO ppi_interfaces 
                        (uniprot_id, partner_id, pdb_id, chains, residues, n_residues,
                         bsa_total, n_atom_contacts, min_atom_distance, stats)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, ppi_iface_batch)
                    ppi_iface_batch = []
                conn.commit()
            
            if n_proteins % args.log_every == 0:
                logging.info(f"  Processed {n_proteins} proteins, {n_residues:,} residues...")
    
    # 残りのバッチを挿入
    if go_batch:
        cursor.executemany(
            "INSERT INTO protein_go (uniprot_id, category, go_term) VALUES (?, ?, ?)",
            go_batch
        )
    if loc_batch:
        cursor.executemany(
            "INSERT INTO protein_locations (uniprot_id, location) VALUES (?, ?)",
            loc_batch
        )
    if ppi_batch:
        cursor.executemany(
            "INSERT INTO ppi_partners (uniprot_id, partner_id) VALUES (?, ?)",
            ppi_batch
        )
    if ppi_iface_batch:
        cursor.executemany("""
            INSERT OR REPLACE INTO ppi_interfaces 
            (uniprot_id, partner_id, pdb_id, chains, residues, n_residues,
             bsa_total, n_atom_contacts, min_atom_distance, stats)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, ppi_iface_batch)
    
    # ac2seq にあるが residue_annotations にないタンパク質を追加
    # （残基データはないが、メタデータとして保持）
    logging.info("Adding proteins from ac2seq that are not in annotations...")
    n_added = 0
    for uniprot_id, seq in ac2seq.items():
        if uniprot_id not in processed_acs:
            description = ac2de.get(uniprot_id, "")
            cursor.execute("""
                INSERT OR IGNORE INTO proteins 
                (uniprot_id, sequence, length, description,
                 n_surface, n_buried, n_unknown, n_disorder, n_interface)
                VALUES (?, ?, ?, ?, 0, 0, 0, 0, 0)
            """, (uniprot_id, seq, len(seq), description))
            
            # GO terms
            if uniprot_id in ac2go:
                go_data = ac2go[uniprot_id]
                for category in ["F", "C", "P"]:
                    for term in go_data.get(category, []):
                        cursor.execute(
                            "INSERT INTO protein_go (uniprot_id, category, go_term) VALUES (?, ?, ?)",
                            (uniprot_id, category, term)
                        )
            
            # Locations
            if uniprot_id in ac2loc:
                for loc in ac2loc[uniprot_id]:
                    cursor.execute(
                        "INSERT INTO protein_locations (uniprot_id, location) VALUES (?, ?)",
                        (uniprot_id, loc)
                    )
            
            n_added += 1
            if n_added % 10000 == 0:
                conn.commit()
                logging.info(f"  Added {n_added} additional proteins...")
    
    conn.commit()
    logging.info(f"  Added {n_added} proteins from ac2seq")
    
    # UniProt features の読み込みと挿入
    if args.uniprot_features:
        logging.info(f"Loading UniProt features: {args.uniprot_features}")
        n_uniprot = 0
        
        residue_feature_batch = []
        disulfide_batch = []
        region_batch = []
        subloc_batch = []
        expression_batch = []
        
        with open_file(args.uniprot_features) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                rec = json.loads(line)
                uniprot_id = rec["uniprot_id"]
                
                # 残基レベル特徴
                residue_features = rec.get("residue_features", {})
                
                # Active sites
                for feat in residue_features.get("active_sites", []):
                    residue_feature_batch.append((
                        uniprot_id,
                        feat.get("position"),
                        "active_site",
                        feat.get("description", "")
                    ))
                
                # Modified residues (PTM)
                for feat in residue_features.get("modified_residues", []):
                    residue_feature_batch.append((
                        uniprot_id,
                        feat.get("position"),
                        "ptm",
                        feat.get("modification", "")
                    ))
                
                # Binding sites
                for feat in residue_features.get("binding_sites", []):
                    pos = feat.get("position")
                    if pos:
                        residue_feature_batch.append((
                            uniprot_id,
                            pos,
                            "binding_site",
                            feat.get("ligand", "")
                        ))
                    else:
                        # 領域として扱う場合
                        start = feat.get("start")
                        end = feat.get("end")
                        if start and end:
                            region_batch.append((
                                uniprot_id,
                                "binding_region",
                                start,
                                end,
                                feat.get("ligand", ""),
                                None
                            ))
                
                # Disulfide bonds
                for feat in residue_features.get("disulfide_bonds", []):
                    disulfide_batch.append((
                        uniprot_id,
                        feat.get("position1"),
                        feat.get("position2"),
                        feat.get("note", "")
                    ))
                
                # 領域レベル特徴
                region_features = rec.get("region_features", {})
                
                # Signal peptides
                for feat in region_features.get("signal_peptides", []):
                    region_batch.append((
                        uniprot_id,
                        "signal_peptide",
                        feat.get("start"),
                        feat.get("end"),
                        "",
                        None
                    ))
                
                # Transit peptides (mitochondria, chloroplast targeting)
                for feat in region_features.get("transit_peptides", []):
                    region_batch.append((
                        uniprot_id,
                        "transit_peptide",
                        feat.get("start"),
                        feat.get("end"),
                        feat.get("description", ""),
                        None
                    ))
                
                # Transmembrane
                for feat in region_features.get("transmembrane", []):
                    region_batch.append((
                        uniprot_id,
                        "transmembrane",
                        feat.get("start"),
                        feat.get("end"),
                        feat.get("description", ""),
                        None
                    ))
                
                # Topological domains (IN/OUT)
                for feat in region_features.get("topological_domains", []):
                    region_batch.append((
                        uniprot_id,
                        "topological_domain",
                        feat.get("start"),
                        feat.get("end"),
                        feat.get("description", ""),
                        feat.get("side")
                    ))
                
                # DNA/RNA binding
                for feat in region_features.get("dna_binding", []):
                    region_batch.append((
                        uniprot_id,
                        "dna_binding",
                        feat.get("start"),
                        feat.get("end"),
                        feat.get("description", ""),
                        None
                    ))
                
                # Coiled-coils
                for feat in region_features.get("coiled_coils", []):
                    region_batch.append((
                        uniprot_id,
                        "coiled_coil",
                        feat.get("start"),
                        feat.get("end"),
                        "",
                        None
                    ))
                
                # Motifs
                for feat in region_features.get("motifs", []):
                    region_batch.append((
                        uniprot_id,
                        "motif",
                        feat.get("start"),
                        feat.get("end"),
                        feat.get("description", ""),
                        None
                    ))
                
                # Domains
                for feat in region_features.get("domains", []):
                    region_batch.append((
                        uniprot_id,
                        "domain",
                        feat.get("start"),
                        feat.get("end"),
                        feat.get("name", ""),
                        None
                    ))
                
                # 細胞内局在
                for loc_chain in rec.get("subcellular_locations", []):
                    subloc_batch.append((
                        uniprot_id,
                        " > ".join(loc_chain)  # ["Cytoplasm", "Cell cortex"] → "Cytoplasm > Cell cortex"
                    ))
                
                # 発現パターン
                expr = rec.get("tissue_expression")
                if expr:
                    expression_batch.append((
                        uniprot_id,
                        expr.get("primary_tissue", ""),
                        expr.get("full_description", "")
                    ))
                
                n_uniprot += 1
                
                # バッチ挿入
                if n_uniprot % 5000 == 0:
                    if residue_feature_batch:
                        cursor.executemany("""
                            INSERT INTO uniprot_residue_features 
                            (uniprot_id, position, feature_type, description)
                            VALUES (?, ?, ?, ?)
                        """, residue_feature_batch)
                        residue_feature_batch = []
                    if disulfide_batch:
                        cursor.executemany("""
                            INSERT INTO uniprot_disulfide_bonds 
                            (uniprot_id, position1, position2, note)
                            VALUES (?, ?, ?, ?)
                        """, disulfide_batch)
                        disulfide_batch = []
                    if region_batch:
                        cursor.executemany("""
                            INSERT INTO uniprot_regions 
                            (uniprot_id, region_type, start_pos, end_pos, description, side)
                            VALUES (?, ?, ?, ?, ?, ?)
                        """, region_batch)
                        region_batch = []
                    if subloc_batch:
                        cursor.executemany("""
                            INSERT INTO uniprot_subcellular_locations 
                            (uniprot_id, location_chain)
                            VALUES (?, ?)
                        """, subloc_batch)
                        subloc_batch = []
                    if expression_batch:
                        cursor.executemany("""
                            INSERT OR REPLACE INTO uniprot_expression 
                            (uniprot_id, primary_tissue, full_description)
                            VALUES (?, ?, ?)
                        """, expression_batch)
                        expression_batch = []
                    conn.commit()
                    logging.info(f"  Processed {n_uniprot} UniProt entries...")
        
        # 残りのバッチを挿入
        if residue_feature_batch:
            cursor.executemany("""
                INSERT INTO uniprot_residue_features 
                (uniprot_id, position, feature_type, description)
                VALUES (?, ?, ?, ?)
            """, residue_feature_batch)
        if disulfide_batch:
            cursor.executemany("""
                INSERT INTO uniprot_disulfide_bonds 
                (uniprot_id, position1, position2, note)
                VALUES (?, ?, ?, ?)
            """, disulfide_batch)
        if region_batch:
            cursor.executemany("""
                INSERT INTO uniprot_regions 
                (uniprot_id, region_type, start_pos, end_pos, description, side)
                VALUES (?, ?, ?, ?, ?, ?)
            """, region_batch)
        if subloc_batch:
            cursor.executemany("""
                INSERT INTO uniprot_subcellular_locations 
                (uniprot_id, location_chain)
                VALUES (?, ?)
            """, subloc_batch)
        if expression_batch:
            cursor.executemany("""
                INSERT OR REPLACE INTO uniprot_expression 
                (uniprot_id, primary_tissue, full_description)
                VALUES (?, ?, ?)
            """, expression_batch)
        
        conn.commit()
        logging.info(f"  Loaded {n_uniprot} UniProt feature entries")
    
    # ファイナライズ
    conn.close()
    total_residues = parquet_writer.close()
    
    # 結果サマリー
    logging.info("=" * 50)
    logging.info("Conversion complete!")
    logging.info(f"SQLite: {args.out_sqlite}")
    logging.info(f"  - {n_proteins + n_added} proteins")
    logging.info(f"Parquet: {args.out_parquet}")
    logging.info(f"  - {total_residues:,} residue records")
    
    # ファイルサイズ
    sqlite_size = Path(args.out_sqlite).stat().st_size / (1024 * 1024)
    parquet_size = Path(args.out_parquet).stat().st_size / (1024 * 1024)
    logging.info(f"File sizes:")
    logging.info(f"  SQLite:  {sqlite_size:.1f} MB")
    logging.info(f"  Parquet: {parquet_size:.1f} MB")


if __name__ == "__main__":
    main()
