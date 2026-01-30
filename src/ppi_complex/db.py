"""
db.py

SQLite + Parquet のヒトプロテオームデータベースにアクセスするユーティリティ。

使用例:
    from ppi_complex.db import HumanProteomeDB
    
    db = HumanProteomeDB(
        sqlite_path="data/processed/human_proteome.sqlite",
        parquet_path="data/processed/human_residues.parquet"
    )
    
    # タンパク質情報を取得
    protein = db.get_protein("P31946")
    print(protein["description"])  # "14-3-3 protein zeta/delta"
    
    # 残基情報を取得
    residue = db.get_residue("P31946", 150)
    print(residue["surface"])  # "surface"
    
    # 複数残基を取得
    residues = db.get_residues("P31946", start=100, end=200)
    
    # 表面残基だけ取得
    surface_residues = db.get_surface_residues("P31946")
    
    # 界面残基を取得
    interface_residues = db.get_interface_residues("P31946")
    
    # 任意のSQLクエリ
    df = db.query("SELECT * FROM proteins WHERE length > 1000 LIMIT 10")
"""

import json
from pathlib import Path
from typing import Optional, List, Dict, Any, Union

try:
    import duckdb
except ImportError:
    raise ImportError("duckdb is required. Install with: pip install duckdb")

try:
    import pandas as pd
except ImportError:
    pd = None


class HumanProteomeDB:
    """SQLite + Parquet のヒトプロテオームDBへのインターフェース"""
    
    def __init__(
        self, 
        sqlite_path: str = "data/processed/human_proteome.sqlite",
        parquet_path: str = "data/processed/human_residues.parquet"
    ):
        """
        Parameters
        ----------
        sqlite_path : str
            SQLiteデータベースのパス
        parquet_path : str
            Parquetファイルのパス
        """
        self.sqlite_path = Path(sqlite_path)
        self.parquet_path = Path(parquet_path)
        
        if not self.sqlite_path.exists():
            raise FileNotFoundError(f"SQLite database not found: {sqlite_path}")
        if not self.parquet_path.exists():
            raise FileNotFoundError(f"Parquet file not found: {parquet_path}")
        
        # DuckDB接続
        self.conn = duckdb.connect()
        
        # SQLiteをアタッチ
        self.conn.execute(f"ATTACH '{self.sqlite_path}' AS meta (TYPE sqlite)")
        
        # Parquetファイルをビューとして登録
        self.conn.execute(f"""
            CREATE OR REPLACE VIEW residues AS 
            SELECT * FROM read_parquet('{self.parquet_path}')
        """)
    
    def close(self):
        """接続を閉じる"""
        self.conn.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    # ========== タンパク質レベルのクエリ ==========
    
    def get_protein(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        タンパク質情報を取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        
        Returns
        -------
        dict or None
            タンパク質情報の辞書、見つからない場合はNone
        """
        result = self.conn.execute("""
            SELECT * FROM meta.proteins WHERE uniprot_id = ?
        """, [uniprot_id]).fetchone()
        
        if result is None:
            return None
        
        columns = [desc[0] for desc in self.conn.description]
        return dict(zip(columns, result))
    
    def get_protein_go(self, uniprot_id: str, category: str = None) -> List[Dict[str, str]]:
        """
        タンパク質のGO termsを取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        category : str, optional
            'F' (Function), 'C' (Component), 'P' (Process) でフィルタ
        
        Returns
        -------
        list of dict
            [{"category": "F", "go_term": "GO:0005515"}, ...]
        """
        if category:
            result = self.conn.execute("""
                SELECT category, go_term FROM meta.protein_go 
                WHERE uniprot_id = ? AND category = ?
            """, [uniprot_id, category]).fetchall()
        else:
            result = self.conn.execute("""
                SELECT category, go_term FROM meta.protein_go 
                WHERE uniprot_id = ?
            """, [uniprot_id]).fetchall()
        
        return [{"category": row[0], "go_term": row[1]} for row in result]
    
    def get_protein_locations(self, uniprot_id: str) -> List[str]:
        """
        タンパク質の細胞内局在を取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        
        Returns
        -------
        list of str
            局在のリスト
        """
        result = self.conn.execute("""
            SELECT location FROM meta.protein_locations 
            WHERE uniprot_id = ?
        """, [uniprot_id]).fetchall()
        
        return [row[0] for row in result]
    
    def get_ppi_partners(self, uniprot_id: str) -> List[str]:
        """
        タンパク質のPPIパートナーを取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        
        Returns
        -------
        list of str
            パートナーのUniProt IDリスト
        """
        result = self.conn.execute("""
            SELECT DISTINCT partner_id FROM meta.ppi_partners 
            WHERE uniprot_id = ?
        """, [uniprot_id]).fetchall()
        
        return [row[0] for row in result]
    
    # ========== 残基レベルのクエリ ==========
    
    def get_residue(self, uniprot_id: str, position: int) -> Optional[Dict[str, Any]]:
        """
        特定位置の残基情報を取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        position : int
            残基位置（1-indexed）
        
        Returns
        -------
        dict or None
            残基情報の辞書
        """
        result = self.conn.execute("""
            SELECT * FROM residues 
            WHERE uniprot_id = ? AND position = ?
        """, [uniprot_id, position]).fetchone()
        
        if result is None:
            return None
        
        columns = [desc[0] for desc in self.conn.description]
        row_dict = dict(zip(columns, result))
        
        # interface_partners をパース
        if "interface_partners" in row_dict:
            row_dict["interface_partners"] = json.loads(row_dict["interface_partners"])
        
        return row_dict
    
    def get_residues(
        self, 
        uniprot_id: str, 
        start: int = None, 
        end: int = None
    ) -> List[Dict[str, Any]]:
        """
        残基情報を範囲で取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        start : int, optional
            開始位置（含む）
        end : int, optional
            終了位置（含む）
        
        Returns
        -------
        list of dict
            残基情報のリスト
        """
        query = "SELECT * FROM residues WHERE uniprot_id = ?"
        params = [uniprot_id]
        
        if start is not None:
            query += " AND position >= ?"
            params.append(start)
        if end is not None:
            query += " AND position <= ?"
            params.append(end)
        
        query += " ORDER BY position"
        
        result = self.conn.execute(query, params).fetchall()
        columns = [desc[0] for desc in self.conn.description]
        
        residues = []
        for row in result:
            row_dict = dict(zip(columns, row))
            if "interface_partners" in row_dict:
                row_dict["interface_partners"] = json.loads(row_dict["interface_partners"])
            residues.append(row_dict)
        
        return residues
    
    def get_surface_residues(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        表面残基のみを取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        
        Returns
        -------
        list of dict
            表面残基のリスト
        """
        result = self.conn.execute("""
            SELECT * FROM residues 
            WHERE uniprot_id = ? AND surface = 'surface'
            ORDER BY position
        """, [uniprot_id]).fetchall()
        
        columns = [desc[0] for desc in self.conn.description]
        
        residues = []
        for row in result:
            row_dict = dict(zip(columns, row))
            if "interface_partners" in row_dict:
                row_dict["interface_partners"] = json.loads(row_dict["interface_partners"])
            residues.append(row_dict)
        
        return residues
    
    def get_interface_residues(
        self, 
        uniprot_id: str, 
        partner_id: str = None
    ) -> List[Dict[str, Any]]:
        """
        界面残基を取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        partner_id : str, optional
            特定パートナーとの界面のみを取得
        
        Returns
        -------
        list of dict
            界面残基のリスト
        """
        result = self.conn.execute("""
            SELECT * FROM residues 
            WHERE uniprot_id = ? AND interface_partners != '[]'
            ORDER BY position
        """, [uniprot_id]).fetchall()
        
        columns = [desc[0] for desc in self.conn.description]
        
        residues = []
        for row in result:
            row_dict = dict(zip(columns, row))
            partners = json.loads(row_dict["interface_partners"])
            row_dict["interface_partners"] = partners
            
            # 特定パートナーでフィルタ
            if partner_id is not None:
                if partner_id not in partners:
                    continue
            
            residues.append(row_dict)
        
        return residues
    
    def get_disordered_residues(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        ディスオーダー残基を取得
        
        Parameters
        ----------
        uniprot_id : str
            UniProt アクセッション番号
        
        Returns
        -------
        list of dict
            ディスオーダー残基のリスト
        """
        result = self.conn.execute("""
            SELECT * FROM residues 
            WHERE uniprot_id = ? AND disorder = true
            ORDER BY position
        """, [uniprot_id]).fetchall()
        
        columns = [desc[0] for desc in self.conn.description]
        
        residues = []
        for row in result:
            row_dict = dict(zip(columns, row))
            if "interface_partners" in row_dict:
                row_dict["interface_partners"] = json.loads(row_dict["interface_partners"])
            residues.append(row_dict)
        
        return residues
    
    # ========== 汎用クエリ ==========
    
    def query(self, sql: str, params: list = None) -> Union["pd.DataFrame", List[tuple]]:
        """
        任意のSQLクエリを実行
        
        Parameters
        ----------
        sql : str
            SQLクエリ（meta.proteins, meta.protein_go, residues 等を使用可能）
        params : list, optional
            パラメータ
        
        Returns
        -------
        pd.DataFrame or list of tuple
            pandasがインストールされていればDataFrame、なければタプルのリスト
        """
        if params:
            result = self.conn.execute(sql, params)
        else:
            result = self.conn.execute(sql)
        
        if pd is not None:
            return result.df()
        else:
            return result.fetchall()
    
    # ========== 統計 ==========
    
    def stats(self) -> Dict[str, int]:
        """
        データベース統計を取得
        
        Returns
        -------
        dict
            統計情報
        """
        n_proteins = self.conn.execute(
            "SELECT COUNT(*) FROM meta.proteins"
        ).fetchone()[0]
        
        n_residues = self.conn.execute(
            "SELECT COUNT(*) FROM residues"
        ).fetchone()[0]
        
        n_with_residues = self.conn.execute("""
            SELECT COUNT(DISTINCT uniprot_id) FROM residues
        """).fetchone()[0]
        
        n_interface_residues = self.conn.execute("""
            SELECT COUNT(*) FROM residues WHERE interface_partners != '[]'
        """).fetchone()[0]
        
        n_disordered = self.conn.execute("""
            SELECT COUNT(*) FROM residues WHERE disorder = true
        """).fetchone()[0]
        
        return {
            "n_proteins": n_proteins,
            "n_proteins_with_residue_data": n_with_residues,
            "n_residues": n_residues,
            "n_interface_residues": n_interface_residues,
            "n_disordered_residues": n_disordered
        }


# 便利関数
def open_db(
    sqlite_path: str = "data/processed/human_proteome.sqlite",
    parquet_path: str = "data/processed/human_residues.parquet"
) -> HumanProteomeDB:
    """
    データベースを開く（コンテキストマネージャとして使用可能）
    
    Example
    -------
    >>> with open_db() as db:
    ...     protein = db.get_protein("P31946")
    """
    return HumanProteomeDB(sqlite_path, parquet_path)
