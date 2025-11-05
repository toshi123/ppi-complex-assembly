"""
ppi_complex

Protein–protein interaction (PPI) network–based complex assembly.

- IntAct の PPI からグラフを構築
- PDB / 相同性検索にもとづく複合体テンプレートを同定
- 共通サブユニットの構造重ね合わせにより複合体を拡張
- インターフェースの抽出と疾患変異・物性解析を行う
"""

__version__ = "0.1.0"

# よく使うサブモジュールをここで import しておくと、
# `import ppi_complex as pc` から `pc.homology.run_mmseqs(...)` みたいに触れる
from . import config
from . import download
from . import homology
from . import assembly
from . import interface
from . import analysis

__all__ = [
    "config",
    "download",
    "homology",
    "assembly",
    "interface",
    "analysis",
]
