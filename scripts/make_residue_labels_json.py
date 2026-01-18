#!/usr/bin/env python
"""
make_residue_labels_json.py

各タンパク質（UniProt ID）について、残基ごとに以下のラベルを付与したJSONLを出力する：
  - surface: 表面残基（単体でも複合体でもRSA > 閾値）
  - interface: 界面残基（単体ではsurface、複合体ではburied）
  - buried: 埋没残基（単体でも複合体でもRSA ≤ 閾値）
  - unknown: SASA計算不可（PDB構造なし等）

処理:
1. TSVからPPIに関与するUniProt IDを収集
2. TSVに複合体構造（pdb_id, chain_id_a, chain_id_b）がある場合は、
   単体と複合体の両方でSASA計算し、差分から界面残基を検出
3. TSVに複合体構造がない場合は、SIFTSから単体チェーンを探してSASA計算
4. RSA（相対SASA）で surface/interface/buried を判定

界面判定ロジック（複合体構造がある場合）:
- 単体RSA > 閾値 かつ 複合体RSA ≤ 閾値 → interface
- 単体RSA > 閾値 かつ 複合体RSA > 閾値 → surface
- 単体RSA ≤ 閾値 → buried（複合体の値に関わらず）

使い方:
python scripts/make_residue_labels_json.py \
  --input  data/interim/intact_pairs_with_pdb_contacts_topk_diverse.tsv \
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

# インターフェース列（参考情報として保持、ラベル付けには使わない）
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


def build_uniprot_to_pdb_index(sifts_df: pd.DataFrame) -> dict:
    """
    SIFTS DataFrame から UniProt ID → [(pdb_id, chain_id, ranges), ...] の逆引きindexを構築。
    UniProtに直接対応するPDB構造を探すために使用。
    """
    idx = defaultdict(list)
    needed = {"PDB", "CHAIN", "SP_PRIMARY", "SP_BEG", "SP_END"}
    if not needed.issubset(set(sifts_df.columns)):
        return idx
    
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
    
    sifts_clean = sifts_df.dropna(subset=required_cols).copy()
    
    # UniProt IDごとにPDB/チェーン情報を集約
    for uniprot_id, group in sifts_clean.groupby("SP_PRIMARY"):
        pdb_chains = {}
        for _, r in group.iterrows():
            pdb_id = str(r["PDB"]).lower()
            chain_id = str(r["CHAIN"])
            key = (pdb_id, chain_id)
            
            pdb_beg = int(r.get("PDB_BEG", r.get("SP_BEG", 1)))
            pdb_end = int(r.get("PDB_END", r.get("SP_END", 1)))
            sp_beg = int(r["SP_BEG"])
            sp_end = int(r["SP_END"])
            
            if key not in pdb_chains:
                pdb_chains[key] = []
            pdb_chains[key].append((pdb_beg, pdb_end, sp_beg, sp_end))
        
        # カバレッジの広いものを優先（sp_end - sp_beg の合計）
        ranked = []
        for (pdb_id, chain_id), ranges in pdb_chains.items():
            coverage = sum(sp_end - sp_beg + 1 for _, _, sp_beg, sp_end in ranges)
            ranges_sorted = sorted(ranges, key=lambda t: t[0])
            ranked.append((coverage, pdb_id, chain_id, ranges_sorted))
        
        # カバレッジ降順でソート
        ranked.sort(key=lambda x: -x[0])
        idx[str(uniprot_id)] = [(pdb_id, chain_id, ranges) for _, pdb_id, chain_id, ranges in ranked]
    
    return idx

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


def write_complex_to_pdb(model: gemmi.Model, chain_ids: list, out_pdb: Path):
    """複数チェーン（複合体）のポリマー残基をPDBファイルに書き出す"""
    st = gemmi.Structure()
    st.name = "complex"
    m = gemmi.Model("1")
    names = {c.name for c in model}
    
    for chain_id in chain_ids:
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


def compute_residue_sasa_by_chain(pdb_path: Path, suppress_stderr: bool) -> dict:
    """
    FreeSASAで残基ごとのSASAを計算（チェーン別）。
    複合体構造で特定チェーンのSASAを取得するために使用。
    
    戻り値: {chain_id: {残基番号(str): {"sasa": float, "aa": str}, ...}, ...}
    """
    with suppress_stderr_fd(suppress_stderr):
        try:
            s = freesasa.Structure(str(pdb_path))
            r = freesasa.calc(s)
        except Exception:
            return {}
    
    res_sasa_by_chain = defaultdict(dict)
    for i in range(s.nAtoms()):
        chain_id = s.chainLabel(i)
        res_num = s.residueNumber(i)  # "123" or "123A"
        res_name = s.residueName(i)   # "ALA", "GLY", etc.
        
        # 3文字コード→1文字コード
        aa1 = gemmi.find_tabulated_residue(res_name).one_letter_code
        if aa1 == 'X':
            aa1 = None
        
        if res_num not in res_sasa_by_chain[chain_id]:
            res_sasa_by_chain[chain_id][res_num] = {"sasa": 0.0, "aa": aa1}
        res_sasa_by_chain[chain_id][res_num]["sasa"] += r.atomArea(i)
    
    return dict(res_sasa_by_chain)


def compute_rsa(sasa: float, aa: str) -> float | None:
    """相対SASA (RSA) を計算"""
    if aa is None or aa not in MAX_SASA:
        return None
    max_sasa = MAX_SASA[aa]
    if max_sasa <= 0:
        return None
    return sasa / max_sasa

# -------------------- メイン処理 --------------------

def collect_ppi_partners(df: pd.DataFrame) -> dict:
    """
    TSVからPPIパートナー情報を収集。
    戻り値: {uniprot_id: set(partner_uniprot_ids), ...}
    """
    partners = defaultdict(set)
    
    # カラム名の正規化（_x, _y サフィックス対応）
    ua_col = None
    ub_col = None
    for col in df.columns:
        if col in ["uniprot_a", "uniprot_a_x", "uniprot_a_y"]:
            ua_col = col
        if col in ["uniprot_b", "uniprot_b_x", "uniprot_b_y"]:
            ub_col = col
    
    if not ua_col or not ub_col:
        logging.warning("No uniprot_a/uniprot_b columns found")
        return partners
    
    for _, r in df.iterrows():
        ua = str(r[ua_col])
        ub = str(r[ub_col])
        partners[ua].add(ub)
        partners[ub].add(ua)
    
    return partners


def build_complex_structure_index(df: pd.DataFrame, sifts_idx: dict) -> dict:
    """
    TSVから各UniProtに対応する複合体構造情報を抽出。
    
    戻り値: {
        uniprot_id: [
            {
                "pdb_id": str,
                "my_chain": str,        # このUniProtのチェーン
                "partner_chain": str,   # パートナーのチェーン
                "partner_uniprot": str,
                "my_ranges": [(pdb_beg, pdb_end, sp_beg, sp_end), ...],
                "partner_ranges": [(pdb_beg, pdb_end, sp_beg, sp_end), ...],
            },
            ...
        ],
        ...
    }
    """
    complex_idx = defaultdict(list)
    
    # カラム名の正規化
    ua_col = ub_col = pdb_col = chain_a_col = chain_b_col = None
    for col in df.columns:
        if col in ["uniprot_a", "uniprot_a_x", "uniprot_a_y"]:
            ua_col = col
        if col in ["uniprot_b", "uniprot_b_x", "uniprot_b_y"]:
            ub_col = col
        if col == "pdb_id":
            pdb_col = col
        if col == "chain_id_a":
            chain_a_col = col
        if col == "chain_id_b":
            chain_b_col = col
    
    if not all([ua_col, ub_col, pdb_col, chain_a_col, chain_b_col]):
        logging.warning("Missing required columns for complex structure index")
        return complex_idx
    
    for _, r in df.iterrows():
        ua = str(r[ua_col])
        ub = str(r[ub_col])
        pdb_id = str(r[pdb_col]).lower()
        chain_a = str(r[chain_a_col])
        chain_b = str(r[chain_b_col])
        
        # SIFTSからマッピング情報を取得
        key_a = (pdb_id, chain_a, ua)
        key_b = (pdb_id, chain_b, ub)
        
        ranges_a = sifts_idx.get(key_a, [])
        ranges_b = sifts_idx.get(key_b, [])
        
        # uniprot_a 側のエントリ
        if ranges_a:  # マッピングがある場合のみ
            complex_idx[ua].append({
                "pdb_id": pdb_id,
                "my_chain": chain_a,
                "partner_chain": chain_b,
                "partner_uniprot": ub,
                "my_ranges": ranges_a,
                "partner_ranges": ranges_b,
            })
        
        # uniprot_b 側のエントリ
        if ranges_b:  # マッピングがある場合のみ
            complex_idx[ub].append({
                "pdb_id": pdb_id,
                "my_chain": chain_b,
                "partner_chain": chain_a,
                "partner_uniprot": ua,
                "my_ranges": ranges_b,
                "partner_ranges": ranges_a,
            })
    
    # 重複除去（同じPDB/チェーンの組み合わせ）
    for uniprot_id in complex_idx:
        seen = set()
        unique = []
        for entry in complex_idx[uniprot_id]:
            key = (entry["pdb_id"], entry["my_chain"], entry["partner_chain"])
            if key not in seen:
                seen.add(key)
                unique.append(entry)
        complex_idx[uniprot_id] = unique
    
    return complex_idx


def collect_target_uniprots(df: pd.DataFrame) -> set:
    """
    TSVからラベル付け対象のUniProt IDを収集。
    """
    uniprots = set()
    
    # カラム名の正規化
    for col in df.columns:
        if col in ["uniprot_a", "uniprot_a_x", "uniprot_a_y"]:
            uniprots.update(df[col].dropna().astype(str).unique())
        if col in ["uniprot_b", "uniprot_b_x", "uniprot_b_y"]:
            uniprots.update(df[col].dropna().astype(str).unique())
    
    return uniprots


def compute_sasa_for_uniprot_monomer(
    uniprot_id: str,
    uniprot_pdb_idx: dict,
    pdb_dir: Path,
    suppress_stderr: bool,
    debug: bool = False,
    max_structures: int = 5
) -> tuple[dict, list]:
    """
    UniProtに**直接対応する**PDB構造からSASAを計算（単体チェーン）。
    SIFTSの逆引きインデックスを使用。
    
    戻り値: (
        {uniprot_pos: {"sasa": float, "rsa": float, "aa": str}, ...},
        [(pdb_id, chain_id), ...]  # 使用したPDBソース
    )
    """
    residue_sasa = {}
    used_sources = []
    tmpdir = Path(".tmp_sasa")
    tmpdir.mkdir(exist_ok=True)
    
    # SIFTSから直接対応するPDB構造を取得
    pdb_sources = uniprot_pdb_idx.get(uniprot_id, [])
    
    if not pdb_sources:
        if debug:
            logging.debug(f"  {uniprot_id}: no direct PDB structures in SIFTS")
        return residue_sasa, used_sources
    
    # カバレッジの高いものから順に処理（最大max_structures個）
    for pdb_id, chain_id, ranges in pdb_sources[:max_structures]:
        struct_path = find_structure_path(pdb_id, pdb_dir)
        if struct_path is None:
            if debug:
                logging.debug(f"  {uniprot_id}: structure file not found for {pdb_id}")
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
        
        if not pdb_res_sasa:
            continue
        
        used_sources.append((pdb_id, chain_id))
        
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
            
            # 既存より大きいRSAを採用（より露出した状態を優先）
            existing_rsa = residue_sasa.get(sp_pos, {}).get("rsa")
            existing_rsa = existing_rsa if existing_rsa is not None else 0
            
            if sp_pos not in residue_sasa or (rsa is not None and rsa > existing_rsa):
                residue_sasa[sp_pos] = {
                    "sasa": sasa,
                    "rsa": rsa,
                    "aa": aa
                }
    
    return residue_sasa, used_sources


def compute_sasa_for_uniprot_complex(
    uniprot_id: str,
    complex_idx: dict,
    pdb_dir: Path,
    suppress_stderr: bool,
    debug: bool = False,
    max_structures: int = 5
) -> tuple[dict, dict, list, str]:
    """
    TSVに記載された複合体構造から、単体と複合体の両方でSASAを計算。
    単体でsurface、複合体でburiedになった残基を界面として検出するために使用。
    
    戻り値: (
        {uniprot_pos: {"sasa": float, "rsa": float, "aa": str}, ...},  # 単体SASA
        {uniprot_pos: {"sasa": float, "rsa": float, "aa": str}, ...},  # 複合体SASA
        [(pdb_id, chain_id, partner_chain), ...]  # 使用したPDBソース
        "complex" or "none"  # 計算モード
    )
    """
    residue_sasa_monomer = {}
    residue_sasa_complex = {}
    used_sources = []
    tmpdir = Path(".tmp_sasa")
    tmpdir.mkdir(exist_ok=True)
    
    # TSVから複合体構造情報を取得
    complex_sources = complex_idx.get(uniprot_id, [])
    
    if not complex_sources:
        if debug:
            logging.debug(f"  {uniprot_id}: no complex structures in TSV")
        return residue_sasa_monomer, residue_sasa_complex, used_sources, "none"
    
    # 複合体構造で計算
    for entry in complex_sources[:max_structures]:
        pdb_id = entry["pdb_id"]
        my_chain = entry["my_chain"]
        partner_chain = entry["partner_chain"]
        my_ranges = entry["my_ranges"]
        
        struct_path = find_structure_path(pdb_id, pdb_dir)
        if struct_path is None:
            if debug:
                logging.debug(f"  {uniprot_id}: structure file not found for {pdb_id}")
            continue
        
        try:
            st = gemmi.read_structure(str(struct_path))
            st.remove_waters()
            model = st[0]
        except Exception as e:
            if debug:
                logging.debug(f"  {uniprot_id}: failed to read {pdb_id}: {e}")
            continue
        
        # 1. 単体（自分のチェーンのみ）でSASA計算
        tmp_pdb_mono = tmpdir / f"{pdb_id}_{my_chain}_mono.pdb"
        sasa_mono_by_chain = {}
        try:
            write_chain_to_pdb(model, my_chain, tmp_pdb_mono)
            sasa_mono_by_chain = compute_residue_sasa_by_chain(tmp_pdb_mono, suppress_stderr)
        except Exception as e:
            if debug:
                logging.debug(f"  {uniprot_id}: monomer SASA calc failed for {pdb_id}: {e}")
        finally:
            try:
                tmp_pdb_mono.unlink()
            except Exception:
                pass
        
        # 2. 複合体（両チェーン）でSASA計算
        chain_ids = [my_chain, partner_chain]
        tmp_pdb_cplx = tmpdir / f"{pdb_id}_{my_chain}_{partner_chain}_complex.pdb"
        sasa_cplx_by_chain = {}
        try:
            write_complex_to_pdb(model, chain_ids, tmp_pdb_cplx)
            sasa_cplx_by_chain = compute_residue_sasa_by_chain(tmp_pdb_cplx, suppress_stderr)
        except Exception as e:
            if debug:
                logging.debug(f"  {uniprot_id}: complex SASA calc failed for {pdb_id}: {e}")
        finally:
            try:
                tmp_pdb_cplx.unlink()
            except Exception:
                pass
        
        # 両方の結果が必要
        if not sasa_mono_by_chain or not sasa_cplx_by_chain:
            if debug:
                logging.debug(f"  {uniprot_id}: incomplete SASA result for {pdb_id}")
            continue
        
        if my_chain not in sasa_mono_by_chain or my_chain not in sasa_cplx_by_chain:
            if debug:
                logging.debug(f"  {uniprot_id}: chain {my_chain} not in SASA result")
            continue
        
        used_sources.append((pdb_id, my_chain, partner_chain))
        
        # 単体SASAの処理
        for pdb_res_num, data in sasa_mono_by_chain[my_chain].items():
            try:
                pdb_pos = int(''.join(c for c in str(pdb_res_num) if c.isdigit()))
            except ValueError:
                continue
            
            sp_pos = pdb_to_sp_pos(pdb_pos, my_ranges)
            if sp_pos is None:
                continue
            
            aa = data["aa"]
            sasa = data["sasa"]
            rsa = compute_rsa(sasa, aa)
            
            # 単体では大きいRSAを採用（より露出した状態を優先）
            existing_rsa = residue_sasa_monomer.get(sp_pos, {}).get("rsa")
            existing_rsa = existing_rsa if existing_rsa is not None else 0
            
            if sp_pos not in residue_sasa_monomer or (rsa is not None and rsa > existing_rsa):
                residue_sasa_monomer[sp_pos] = {
                    "sasa": sasa,
                    "rsa": rsa,
                    "aa": aa
                }
        
        # 複合体SASAの処理
        for pdb_res_num, data in sasa_cplx_by_chain[my_chain].items():
            try:
                pdb_pos = int(''.join(c for c in str(pdb_res_num) if c.isdigit()))
            except ValueError:
                continue
            
            sp_pos = pdb_to_sp_pos(pdb_pos, my_ranges)
            if sp_pos is None:
                continue
            
            aa = data["aa"]
            sasa = data["sasa"]
            rsa = compute_rsa(sasa, aa)
            
            # 複合体では小さいRSAを採用（界面での埋没状態を優先）
            existing_rsa = residue_sasa_complex.get(sp_pos, {}).get("rsa")
            
            if sp_pos not in residue_sasa_complex:
                residue_sasa_complex[sp_pos] = {
                    "sasa": sasa,
                    "rsa": rsa,
                    "aa": aa
                }
            elif rsa is not None and (existing_rsa is None or rsa < existing_rsa):
                residue_sasa_complex[sp_pos] = {
                    "sasa": sasa,
                    "rsa": rsa,
                    "aa": aa
                }
    
    return residue_sasa_monomer, residue_sasa_complex, used_sources, "complex"


def compute_sasa_for_uniprot(
    uniprot_id: str,
    complex_idx: dict,
    uniprot_pdb_idx: dict,
    pdb_dir: Path,
    suppress_stderr: bool,
    debug: bool = False,
    max_structures: int = 5
) -> tuple[dict, dict, list, str]:
    """
    UniProtのSASAを計算。
    1. TSVに複合体構造がある場合は単体と複合体の両方で計算（界面検出用）
    2. TSVに複合体構造がない場合はSIFTSから単体チェーンのみで計算
    
    戻り値: (
        {uniprot_pos: {"sasa": float, "rsa": float, "aa": str}, ...},  # 単体SASA
        {uniprot_pos: {"sasa": float, "rsa": float, "aa": str}, ...},  # 複合体SASA（なければ空dict）
        [使用したPDBソース],
        "complex" or "monomer" or "none"  # 計算モード
    )
    """
    # まず複合体構造で計算を試みる
    sasa_mono, sasa_cplx, used_sources, mode = compute_sasa_for_uniprot_complex(
        uniprot_id, complex_idx, pdb_dir, suppress_stderr, debug, max_structures
    )
    
    if sasa_mono and sasa_cplx:
        return sasa_mono, sasa_cplx, used_sources, "complex"
    
    # 複合体構造がない場合は単体チェーンのみで計算
    residue_sasa, used_sources = compute_sasa_for_uniprot_monomer(
        uniprot_id, uniprot_pdb_idx, pdb_dir, suppress_stderr, debug, max_structures
    )
    
    if residue_sasa:
        return residue_sasa, {}, used_sources, "monomer"
    
    return {}, {}, [], "none"


def main():
    args = parse_args()
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=log_level, format="%(levelname)s: %(message)s")
    
    # 入力TSV読み込み
    logging.info(f"Reading input: {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    
    # SIFTS読み込み（gzip自動判定）
    logging.info(f"Reading SIFTS: {args.sifts}")
    sifts_df = read_csv_auto(args.sifts, comment='#', low_memory=False)
    
    # (pdb, chain, uniprot) -> ranges のインデックスを構築
    logging.info("Building SIFTS index...")
    sifts_idx = build_sifts_index(sifts_df.copy())
    logging.info(f"  Found {len(sifts_idx)} (pdb, chain, uniprot) entries in SIFTS")
    
    # UniProt → PDB の逆引きインデックスを構築（単体用フォールバック）
    logging.info("Building UniProt -> PDB index from SIFTS...")
    uniprot_pdb_idx = build_uniprot_to_pdb_index(sifts_df.copy())
    logging.info(f"  Found direct PDB structures for {len(uniprot_pdb_idx)} UniProt IDs in SIFTS")
    
    # TSVから複合体構造インデックスを構築
    logging.info("Building complex structure index from TSV...")
    complex_idx = build_complex_structure_index(df, sifts_idx)
    logging.info(f"  Found complex structures for {len(complex_idx)} UniProt IDs in TSV")
    
    # UniProt配列読み込み
    logging.info(f"Reading UniProt FASTA: {args.fasta}")
    up_seqs = load_uniprot_fasta(args.fasta)
    
    # PPIパートナー情報を収集（メタデータ用）
    logging.info("Collecting PPI partner info from TSV...")
    ppi_partners = collect_ppi_partners(df)
    logging.info(f"  Found PPI info for {len(ppi_partners)} UniProt IDs")
    
    # 対象UniProtリスト
    target_uniprots = collect_target_uniprots(df)
    logging.info(f"Processing {len(target_uniprots)} UniProt IDs...")
    
    # 出力
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_dir = Path(args.pdb_dir)
    
    n_written = 0
    n_complex = 0
    n_monomer = 0
    n_no_sasa = 0
    
    with open(out_path, "w") as fw:
        for idx, uniprot_id in enumerate(sorted(target_uniprots)):
            if (idx + 1) % args.log_every == 0:
                logging.info(f"  Processed {idx + 1}/{len(target_uniprots)}")
            
            seq = up_seqs.get(uniprot_id)
            if seq is None:
                continue
            
            # SASA計算（複合体優先、なければ単体）
            sasa_mono, sasa_cplx, used_sources, calc_mode = compute_sasa_for_uniprot(
                uniprot_id, complex_idx, uniprot_pdb_idx, pdb_dir, args.suppress_stderr,
                debug=args.debug
            )
            
            if calc_mode == "complex":
                n_complex += 1
            elif calc_mode == "monomer":
                n_monomer += 1
            else:
                n_no_sasa += 1
            
            # PPIパートナー（メタデータ）
            partners = list(ppi_partners.get(uniprot_id, set()))
            
            # 各残基にラベル付け（surface/interface/buried/unknown）
            residue_labels = {}
            n_surface = 0
            n_interface = 0
            n_buried = 0
            
            for pos in range(1, len(seq) + 1):
                aa = seq[pos - 1]
                
                # 単体SASAを基準に使用
                mono_info = sasa_mono.get(pos, {})
                rsa_mono = mono_info.get("rsa")
                sasa_mono_val = mono_info.get("sasa")
                
                # 複合体SASA（ある場合）
                cplx_info = sasa_cplx.get(pos, {})
                rsa_cplx = cplx_info.get("rsa")
                sasa_cplx_val = cplx_info.get("sasa")
                
                # ラベル決定（surface/interface/buried/unknown）
                if rsa_mono is not None:
                    if calc_mode == "complex" and rsa_cplx is not None:
                        # 複合体計算モード：単体と複合体の差分から界面を検出
                        if rsa_mono > args.rsa_threshold:
                            # 単体でsurface
                            if rsa_cplx <= args.rsa_threshold:
                                # 複合体でburied → interface
                                label = "interface"
                                n_interface += 1
                            else:
                                # 複合体でもsurface → surface
                                label = "surface"
                                n_surface += 1
                        else:
                            # 単体でburied → buried（複合体の値に関わらず）
                            label = "buried"
                            n_buried += 1
                    else:
                        # 単体計算モード：単体RSAのみで判定
                        if rsa_mono > args.rsa_threshold:
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
                # 単体RSAを記録（基準値として）
                if rsa_mono is not None:
                    residue_labels[str(pos)]["rsa_mono"] = round(rsa_mono, 3)
                if sasa_mono_val is not None:
                    residue_labels[str(pos)]["sasa_mono"] = round(sasa_mono_val, 2)
                # 複合体RSAも記録（ある場合）
                if rsa_cplx is not None:
                    residue_labels[str(pos)]["rsa_complex"] = round(rsa_cplx, 3)
                if sasa_cplx_val is not None:
                    residue_labels[str(pos)]["sasa_complex"] = round(sasa_cplx_val, 2)
            
            # サマリ
            n_total = len(seq)
            n_unknown = n_total - n_surface - n_interface - n_buried
            
            # used_sources のフォーマット（complex vs monomer で異なる）
            if calc_mode == "complex":
                pdb_sources_str = [f"{p}_{c}+{pc}" for p, c, pc in used_sources]
            else:
                pdb_sources_str = [f"{p}_{c}" for p, c in used_sources]
            
            rec = {
                "uniprot_id": uniprot_id,
                "sequence": seq,
                "length": n_total,
                "residue_labels": residue_labels,
                "summary": {
                    "n_total": n_total,
                    "n_surface": n_surface,
                    "n_interface": n_interface,
                    "n_buried": n_buried,
                    "n_unknown": n_unknown,
                    "surface_frac": round(n_surface / n_total, 3) if n_total > 0 else 0,
                    "interface_frac": round(n_interface / n_total, 3) if n_total > 0 else 0,
                    "buried_frac": round(n_buried / n_total, 3) if n_total > 0 else 0,
                },
                "sasa_calc_mode": calc_mode,
                "pdb_sources": pdb_sources_str,
                "ppi_info": {
                    "has_ppi": len(partners) > 0,
                    "n_partners": len(partners),
                    "partners": sorted(partners),
                    "interface_note": (
                        "Interface residues detected by comparing monomer vs complex SASA."
                        if calc_mode == "complex" else
                        "Monomer SASA only; interface positions cannot be determined."
                        if calc_mode == "monomer" else None
                    ) if partners else None
                },
            }
            
            fw.write(json.dumps(rec, ensure_ascii=False) + "\n")
            n_written += 1
    
    logging.info(f"Written: {out_path}")
    logging.info(f"  Total proteins: {n_written}")
    logging.info(f"  With SASA (complex): {n_complex}")
    logging.info(f"  With SASA (monomer): {n_monomer}")
    logging.info(f"  No SASA: {n_no_sasa}")


if __name__ == "__main__":
    main()
