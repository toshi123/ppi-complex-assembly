#!/usr/bin/env python
"""
make_residue_labels_json.py

各タンパク質（UniProt ID）について、残基ごとに以下のラベルを付与したJSONLを出力する：
  - surface: 表面残基（RSA > 閾値）
  - buried: 埋没残基（RSA ≤ 閾値）
  - interface: インターフェース残基（PPIで接触）

処理:
1. TSVからインターフェース残基を収集（全相互作用相手を統合）
2. PDB構造から単独チェーンのSASAを計算（FreeSASA）
3. RSA（相対SASA）で surface/buried を判定
4. interface はTSVの iface_res_*/iface_sasa_* から取得
5. ラベルの優先順位: interface > surface > buried

使い方:
python scripts/make_residue_labels_json.py \
  --input  data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
  --upstream data/interim/intact_pairs_with_pdb_contacts.tight_bsa.band.tsv \
  --sifts   data/raw/pdb_chain_uniprot_current.csv.gz \
  --fasta   data/raw/uniprot_human_all_20251110.fasta \
  --pdb-dir data/raw/pdb_structures \
  --output  data/processed/residue_labels.jsonl
"""

import argparse
import gzip
import json
import logging
import os
import sys
from pathlib import Path
from collections import defaultdict
from contextlib import contextmanager

import pandas as pd
import numpy as np
import gemmi
import freesasa

# -------------------- 定数 --------------------

# Gly-X-Gly トリペプチド中での最大SASA（Miller et al. 1987 の値をベースに調整）
# RSA = SASA / MAX_SASA[aa] で相対SASAを計算
MAX_SASA = {
    'A': 113.0, 'R': 241.0, 'N': 158.0, 'D': 151.0, 'C': 140.0,
    'Q': 189.0, 'E': 183.0, 'G': 85.0,  'H': 194.0, 'I': 182.0,
    'L': 180.0, 'K': 211.0, 'M': 204.0, 'F': 218.0, 'P': 143.0,
    'S': 122.0, 'T': 146.0, 'W': 259.0, 'Y': 229.0, 'V': 160.0,
}

CANDIDATE_IFACE_A = ["iface_sasa_a", "iface_res_a"]
CANDIDATE_IFACE_B = ["iface_sasa_b", "iface_res_b"]

# -------------------- ファイル読み込みユーティリティ --------------------

def resolve_symlink(path: Path) -> Path:
    """
    シンボリックリンクを解決する。
    macOSで作成されたシンボリックリンクがLinuxで壊れている場合の検出も行う。
    """
    path = Path(path)
    
    # シンボリックリンクの場合は解決を試みる
    if path.is_symlink():
        try:
            resolved = path.resolve()
            if resolved.exists():
                return resolved
        except Exception:
            pass
    
    # ファイルの先頭を確認（macOSシンボリックリンクのメタデータ検出）
    if path.exists():
        try:
            with open(path, 'rb') as f:
                header = f.read(10)
            if header.startswith(b'XSym'):
                # macOSシンボリックリンクのメタファイル
                # 実際のリンク先を読み取る
                with open(path, 'r') as f:
                    lines = f.readlines()
                if len(lines) >= 4:
                    target = lines[3].strip()
                    target_path = path.parent / target
                    if target_path.exists():
                        logging.warning(f"Resolved macOS symlink: {path} -> {target_path}")
                        return target_path
                raise ValueError(
                    f"File appears to be a macOS symlink metadata file: {path}\n"
                    f"Please specify the actual file path instead of the symlink.\n"
                    f"Try: ls -la {path.parent}/ to find the actual file."
                )
        except ValueError:
            raise
        except Exception:
            pass
    
    return path


def read_csv_auto(path, **kwargs):
    """
    gzip圧縮かどうかをマジックバイトで自動判定してCSVを読み込む。
    拡張子が.gzでも実際に非圧縮の場合に対応。
    シンボリックリンクも解決する。
    """
    path = resolve_symlink(Path(path))
    
    # マジックバイトでgzipかどうか判定
    with open(path, 'rb') as f:
        magic = f.read(2)
    
    if magic == b'\x1f\x8b':  # gzip magic number
        return pd.read_csv(path, compression='gzip', **kwargs)
    else:
        return pd.read_csv(path, compression=None, **kwargs)

# -------------------- CLI --------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Build per-protein residue labels (surface/buried/interface) JSONL."
    )
    p.add_argument("--input", required=True,
                   help="Main TSV (e.g., topk_diverse.tsv or tight_bsa.band.tsv)")
    p.add_argument("--upstream", nargs="*", default=[],
                   help="Additional TSVs to merge interface residue info")
    p.add_argument("--sifts", required=True,
                   help="SIFTS chain mapping CSV (pdb_chain_uniprot_current.csv.gz)")
    p.add_argument("--fasta", required=True,
                   help="UniProt FASTA (human)")
    p.add_argument("--pdb-dir", default="data/raw/pdb_structures",
                   help="Directory containing PDB structure files (.cif.gz)")
    p.add_argument("--output", default="data/processed/residue_labels.jsonl",
                   help="Output JSONL")
    p.add_argument("--rsa-threshold", type=float, default=0.25,
                   help="RSA threshold for surface/buried (default: 0.25)")
    p.add_argument("--n-proc", type=int, default=1,
                   help="Number of parallel processes (not yet implemented)")
    p.add_argument("--log-every", type=int, default=100)
    p.add_argument("--suppress-stderr", action="store_true", default=True,
                   help="Suppress FreeSASA stderr")
    p.add_argument("--debug", action="store_true",
                   help="Enable debug logging for SASA computation")
    return p.parse_args()

# -------------------- stderr抑制 --------------------

@contextmanager
def suppress_stderr_fd(active: bool = True):
    if not active:
        yield
        return
    try:
        stderr_fd = sys.stderr.fileno()
    except Exception:
        yield
        return
    with open(os.devnull, "w") as devnull:
        old_stderr_fd = os.dup(stderr_fd)
        try:
            os.dup2(devnull.fileno(), stderr_fd)
            yield
        finally:
            try:
                os.dup2(old_stderr_fd, stderr_fd)
            finally:
                os.close(old_stderr_fd)

# -------------------- IO --------------------

def extract_uniprot_ac(s: str) -> str:
    s = str(s).strip()
    if '|' in s:
        parts = s.split('|')
        if len(parts) >= 2:
            return parts[1]
    return s


def load_uniprot_fasta(fasta_path: str) -> dict:
    seqs = {}
    def _open(path):
        return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")
    with _open(fasta_path) as f:
        ac = None
        buf = []
        for line in f:
            if line.startswith(">"):
                if ac and buf:
                    seqs[ac] = "".join(buf).replace("\n", "").strip()
                header = line[1:].strip()
                raw_ac = header.split()[0]
                ac = extract_uniprot_ac(raw_ac)
                buf = []
            else:
                buf.append(line.strip())
        if ac and buf:
            seqs[ac] = "".join(buf).replace("\n", "").strip()
    return seqs


def normalize_res_list(x):
    """カンマ区切りの残基番号リストをintリストに変換"""
    if pd.isna(x):
        return []
    s = str(x).replace(";", ",").replace(" ", ",")
    toks = [t for t in s.split(",") if t.strip() != ""]
    out = []
    for t in toks:
        t = t.strip()
        num = ""
        for ch in t:
            if ch.isdigit():
                num += ch
            else:
                break
        if num:
            out.append(int(num))
    return sorted(set(out))


def pick_first_present(cols, df):
    for c in cols:
        if c in df.columns:
            return c
    return None


def find_structure_path(pdb_id: str, pdb_dir: Path) -> Path | None:
    pdb_id = pdb_id.lower()
    for ext in (".bcif.gz", ".cif.gz", ".cif", ".bcif"):
        p = pdb_dir / f"{pdb_id}{ext}"
        if p.exists():
            return p
    return None

# -------------------- SIFTS --------------------

def build_sifts_index(sifts_df: pd.DataFrame):
    """SIFTS DataFrame から (pdb, chain, uniprot) -> [(pdb_beg, pdb_end, sp_beg, sp_end), ...] のindexを構築"""
    idx = defaultdict(list)
    needed = {"PDB", "CHAIN", "SP_PRIMARY", "SP_BEG", "SP_END"}
    if not needed.issubset(set(sifts_df.columns)):
        raise ValueError(f"SIFTS columns missing. Need {needed}, got {list(sifts_df.columns)}")
    
    for col in ["RES_BEG", "RES_END", "PDB_BEG", "PDB_END", "SP_BEG", "SP_END"]:
        if col in sifts_df.columns:
            sifts_df[col] = pd.to_numeric(sifts_df[col], errors="coerce")
    
    # PDB_BEG/END のフォールバック
    if "RES_BEG" in sifts_df.columns and "PDB_BEG" in sifts_df.columns:
        sifts_df["PDB_BEG"] = sifts_df["PDB_BEG"].fillna(sifts_df["RES_BEG"])
    if "RES_END" in sifts_df.columns and "PDB_END" in sifts_df.columns:
        sifts_df["PDB_END"] = sifts_df["PDB_END"].fillna(sifts_df["RES_END"])
    
    required_cols = ["PDB", "CHAIN", "SP_PRIMARY", "SP_BEG", "SP_END"]
    if "PDB_BEG" in sifts_df.columns and "PDB_END" in sifts_df.columns:
        required_cols += ["PDB_BEG", "PDB_END"]
    
    sifts_df = sifts_df.dropna(subset=required_cols)
    
    for _, r in sifts_df.iterrows():
        key = (str(r["PDB"]).lower(), str(r["CHAIN"]), str(r["SP_PRIMARY"]))
        pdb_beg = int(r.get("PDB_BEG", r.get("SP_BEG", 1)))
        pdb_end = int(r.get("PDB_END", r.get("SP_END", 1)))
        sp_beg = int(r["SP_BEG"])
        sp_end = int(r["SP_END"])
        idx[key].append((pdb_beg, pdb_end, sp_beg, sp_end))
    
    for k in idx:
        idx[k].sort(key=lambda t: t[0])
    return idx


def pdb_to_sp_pos(pdb_pos: int, ranges) -> int | None:
    """PDB残基番号 → UniProt残基番号"""
    for pdb_beg, pdb_end, sp_beg, sp_end in ranges:
        if pdb_beg <= pdb_pos <= pdb_end:
            sp = sp_beg + (pdb_pos - pdb_beg)
            if sp <= sp_end:
                return sp
    return None

# -------------------- SASA計算 --------------------

def is_polymer_residue(res: gemmi.Residue) -> bool:
    """
    残基がポリマー（タンパク質/核酸）かどうかを判定。
    gemmiのバージョン互換性を考慮。
    """
    # 方法1: is_polymer() メソッドがあれば使う
    if hasattr(res, 'is_polymer'):
        try:
            return res.is_polymer()
        except Exception:
            pass
    
    # 方法2: het_flag で判定（'A' = ATOM = ポリマー、'H' = HETATM = 非ポリマー）
    if hasattr(res, 'het_flag'):
        return res.het_flag == 'A'
    
    # 方法3: 残基名で判定（標準アミノ酸 + 核酸）
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT', 'DU',
        'MSE', 'SEC', 'PYL',  # 修飾アミノ酸
    }
    return res.name in standard_residues


def write_chain_to_pdb(model: gemmi.Model, chain_id: str, out_pdb: Path):
    """指定チェーンのポリマー残基のみをPDBファイルに書き出す"""
    st = gemmi.Structure()
    st.name = "sel"
    m = gemmi.Model("1")
    names = {c.name for c in model}
    if chain_id in names:
        new_c = gemmi.Chain(chain_id)
        for res in model[chain_id]:
            if is_polymer_residue(res):
                new_c.add_residue(res)
        if len(new_c):
            m.add_chain(new_c)
    st.add_model(m)
    st.write_minimal_pdb(str(out_pdb))


def compute_residue_sasa(pdb_path: Path, suppress_stderr: bool) -> dict:
    """
    FreeSASAで残基ごとのSASAを計算。
    戻り値: {残基番号(str): {"sasa": float, "aa": str}, ...}
    """
    with suppress_stderr_fd(suppress_stderr):
        try:
            s = freesasa.Structure(str(pdb_path))
            r = freesasa.calc(s)
        except Exception:
            return {}
    
    res_sasa = {}
    for i in range(s.nAtoms()):
        res_num = s.residueNumber(i)  # "123" or "123A"
        res_name = s.residueName(i)   # "ALA", "GLY", etc.
        
        # 3文字コード→1文字コード
        aa1 = gemmi.find_tabulated_residue(res_name).one_letter_code
        if aa1 == 'X':
            aa1 = None
        
        if res_num not in res_sasa:
            res_sasa[res_num] = {"sasa": 0.0, "aa": aa1}
        res_sasa[res_num]["sasa"] += r.atomArea(i)
    
    return res_sasa


def compute_rsa(sasa: float, aa: str) -> float | None:
    """相対SASA (RSA) を計算"""
    if aa is None or aa not in MAX_SASA:
        return None
    max_sasa = MAX_SASA[aa]
    if max_sasa <= 0:
        return None
    return sasa / max_sasa

# -------------------- メイン処理 --------------------

def collect_interface_residues(df: pd.DataFrame, sifts_idx: dict) -> dict:
    """
    TSVからインターフェース残基を収集し、UniProtごとに統合。
    戻り値: {uniprot_id: set(uniprot_pos), ...}
    """
    iface_a_col = pick_first_present(CANDIDATE_IFACE_A, df)
    iface_b_col = pick_first_present(CANDIDATE_IFACE_B, df)
    
    if not iface_a_col or not iface_b_col:
        logging.warning("No interface residue columns found in input")
        return {}
    
    uniprot_iface = defaultdict(set)
    
    for _, r in df.iterrows():
        ua = str(r["uniprot_a"])
        ub = str(r["uniprot_b"])
        pdb = str(r["pdb_id"]).lower()
        ca = str(r["chain_id_a"])
        cb = str(r["chain_id_b"])
        
        # A側
        ra_pdb = normalize_res_list(r.get(iface_a_col))
        ranges_a = sifts_idx.get((pdb, ca, ua), [])
        for pdb_pos in ra_pdb:
            sp_pos = pdb_to_sp_pos(pdb_pos, ranges_a)
            if sp_pos:
                uniprot_iface[ua].add(sp_pos)
        
        # B側
        rb_pdb = normalize_res_list(r.get(iface_b_col))
        ranges_b = sifts_idx.get((pdb, cb, ub), [])
        for pdb_pos in rb_pdb:
            sp_pos = pdb_to_sp_pos(pdb_pos, ranges_b)
            if sp_pos:
                uniprot_iface[ub].add(sp_pos)
    
    return uniprot_iface


def collect_pdb_sources(df: pd.DataFrame) -> dict:
    """
    UniProtごとに使用可能なPDBソースを収集。
    戻り値: {uniprot_id: [(pdb_id, chain_id), ...], ...}
    """
    sources = defaultdict(list)
    
    for _, r in df.iterrows():
        ua = str(r["uniprot_a"])
        ub = str(r["uniprot_b"])
        pdb = str(r["pdb_id"]).lower()
        ca = str(r["chain_id_a"])
        cb = str(r["chain_id_b"])
        
        sources[ua].append((pdb, ca))
        sources[ub].append((pdb, cb))
    
    # 重複除去
    for k in sources:
        sources[k] = list(set(sources[k]))
    
    return sources


def compute_sasa_for_uniprot(
    uniprot_id: str,
    pdb_sources: list,
    sifts_idx: dict,
    pdb_dir: Path,
    suppress_stderr: bool,
    debug: bool = False
) -> dict:
    """
    UniProtの各残基位置についてSASAを計算（複数PDBソースから最良値を採用）。
    戻り値: {uniprot_pos: {"sasa": float, "rsa": float, "aa": str}, ...}
    """
    residue_sasa = {}
    tmpdir = Path(".tmp_sasa")
    tmpdir.mkdir(exist_ok=True)
    
    if not pdb_sources:
        if debug:
            logging.debug(f"  {uniprot_id}: no PDB sources")
        return residue_sasa
    
    for pdb_id, chain_id in pdb_sources:
        struct_path = find_structure_path(pdb_id, pdb_dir)
        if struct_path is None:
            if debug:
                logging.debug(f"  {uniprot_id}: structure not found for {pdb_id}")
            continue
        
        try:
            st = gemmi.read_structure(str(struct_path))
            st.remove_waters()
            model = st[0]
        except Exception as e:
            if debug:
                logging.debug(f"  {uniprot_id}: failed to read {pdb_id}: {e}")
            continue
        
        # チェーンをPDBに書き出し
        tmp_pdb = tmpdir / f"{pdb_id}_{chain_id}.pdb"
        try:
            write_chain_to_pdb(model, chain_id, tmp_pdb)
            pdb_res_sasa = compute_residue_sasa(tmp_pdb, suppress_stderr)
            if debug and not pdb_res_sasa:
                logging.debug(f"  {uniprot_id}: empty SASA result for {pdb_id}_{chain_id}")
        except Exception as e:
            if debug:
                logging.debug(f"  {uniprot_id}: SASA calc failed for {pdb_id}_{chain_id}: {e}")
            continue
        finally:
            try:
                tmp_pdb.unlink()
            except Exception:
                pass
        
        # SIFTSでUniProt位置にマッピング
        ranges = sifts_idx.get((pdb_id, chain_id, uniprot_id), [])
        if debug and not ranges:
            logging.debug(f"  {uniprot_id}: no SIFTS ranges for ({pdb_id}, {chain_id}, {uniprot_id})")
        if not ranges:
            continue
        
        for pdb_res_num, data in pdb_res_sasa.items():
            # pdb_res_num は "123" のような文字列
            try:
                pdb_pos = int(''.join(c for c in pdb_res_num if c.isdigit()))
            except ValueError:
                continue
            
            sp_pos = pdb_to_sp_pos(pdb_pos, ranges)
            if sp_pos is None:
                continue
            
            aa = data["aa"]
            sasa = data["sasa"]
            rsa = compute_rsa(sasa, aa)
            
            # 既存より大きいSASAを採用（より露出した状態を優先）
            if sp_pos not in residue_sasa or (rsa and rsa > residue_sasa[sp_pos].get("rsa", 0)):
                residue_sasa[sp_pos] = {
                    "sasa": sasa,
                    "rsa": rsa,
                    "aa": aa
                }
    
    return residue_sasa


def main():
    args = parse_args()
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=log_level, format="%(levelname)s: %(message)s")
    
    # 入力TSV読み込み
    logging.info(f"Reading input: {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    
    # 上流TSVをマージ（インターフェース残基情報を補完）
    for up_path in args.upstream:
        if not Path(up_path).exists():
            logging.warning(f"Upstream not found: {up_path}")
            continue
        up_df = pd.read_csv(up_path, sep="\t")
        merge_keys = ["uniprot_a", "uniprot_b", "pdb_id", "chain_id_a", "chain_id_b"]
        if not set(merge_keys).issubset(up_df.columns):
            continue
        
        # 不足列を補完
        for col in CANDIDATE_IFACE_A + CANDIDATE_IFACE_B:
            if col not in df.columns and col in up_df.columns:
                up_sub = up_df[merge_keys + [col]].drop_duplicates(subset=merge_keys)
                df = df.merge(up_sub, on=merge_keys, how="left", suffixes=("", "_up"))
                if f"{col}_up" in df.columns:
                    df[col] = df[col].where(df[col].notna(), df[f"{col}_up"])
                    df = df.drop(columns=[f"{col}_up"])
    
    # SIFTS読み込み（gzip自動判定）
    logging.info(f"Reading SIFTS: {args.sifts}")
    sifts_df = read_csv_auto(args.sifts, comment='#', low_memory=False)
    sifts_idx = build_sifts_index(sifts_df)
    
    # UniProt配列読み込み
    logging.info(f"Reading UniProt FASTA: {args.fasta}")
    up_seqs = load_uniprot_fasta(args.fasta)
    
    # インターフェース残基を収集
    logging.info("Collecting interface residues from TSV...")
    uniprot_iface = collect_interface_residues(df, sifts_idx)
    logging.info(f"  Found interface info for {len(uniprot_iface)} UniProt IDs")
    
    # PDBソースを収集
    pdb_sources = collect_pdb_sources(df)
    
    # 対象UniProtリスト
    all_uniprots = set(uniprot_iface.keys()) | set(pdb_sources.keys())
    logging.info(f"Processing {len(all_uniprots)} UniProt IDs...")
    
    # 出力
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_dir = Path(args.pdb_dir)
    
    n_written = 0
    n_no_sasa = 0
    
    with open(out_path, "w") as fw:
        for idx, uniprot_id in enumerate(sorted(all_uniprots)):
            if (idx + 1) % args.log_every == 0:
                logging.info(f"  Processed {idx + 1}/{len(all_uniprots)}")
            
            seq = up_seqs.get(uniprot_id)
            if seq is None:
                continue
            
            # SASA計算
            sources = pdb_sources.get(uniprot_id, [])
            residue_sasa = compute_sasa_for_uniprot(
                uniprot_id, sources, sifts_idx, pdb_dir, args.suppress_stderr,
                debug=args.debug
            )
            
            if not residue_sasa:
                n_no_sasa += 1
            
            # インターフェース残基
            iface_positions = uniprot_iface.get(uniprot_id, set())
            
            # 各残基にラベル付け
            residue_labels = {}
            n_surface = 0
            n_buried = 0
            n_interface = 0
            
            for pos in range(1, len(seq) + 1):
                aa = seq[pos - 1]
                sasa_info = residue_sasa.get(pos, {})
                rsa = sasa_info.get("rsa")
                sasa = sasa_info.get("sasa")
                
                # ラベル決定（優先順位: interface > surface > buried）
                if pos in iface_positions:
                    label = "interface"
                    n_interface += 1
                elif rsa is not None:
                    if rsa > args.rsa_threshold:
                        label = "surface"
                        n_surface += 1
                    else:
                        label = "buried"
                        n_buried += 1
                else:
                    label = "unknown"
                
                residue_labels[str(pos)] = {
                    "aa": aa,
                    "label": label,
                }
                if rsa is not None:
                    residue_labels[str(pos)]["rsa"] = round(rsa, 3)
                if sasa is not None:
                    residue_labels[str(pos)]["sasa"] = round(sasa, 2)
            
            # サマリ
            n_total = len(seq)
            n_unknown = n_total - n_surface - n_buried - n_interface
            
            rec = {
                "uniprot_id": uniprot_id,
                "sequence": seq,
                "length": n_total,
                "residue_labels": residue_labels,
                "summary": {
                    "n_total": n_total,
                    "n_surface": n_surface,
                    "n_buried": n_buried,
                    "n_interface": n_interface,
                    "n_unknown": n_unknown,
                    "surface_frac": round(n_surface / n_total, 3) if n_total > 0 else 0,
                    "buried_frac": round(n_buried / n_total, 3) if n_total > 0 else 0,
                    "interface_frac": round(n_interface / n_total, 3) if n_total > 0 else 0,
                },
                "pdb_sources": [f"{p}_{c}" for p, c in sources],
            }
            
            fw.write(json.dumps(rec, ensure_ascii=False) + "\n")
            n_written += 1
    
    logging.info(f"Written: {out_path}")
    logging.info(f"  Total proteins: {n_written}")
    logging.info(f"  No SASA computed: {n_no_sasa}")


if __name__ == "__main__":
    main()
