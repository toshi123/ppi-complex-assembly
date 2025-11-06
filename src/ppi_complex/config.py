"""
設定ファイルの読み込みと管理
"""
import yaml
from pathlib import Path
from typing import Dict, Any

# リポジトリのルートディレクトリを取得
REPO_ROOT = Path(__file__).resolve().parents[2]


def load_config(config_name: str = "project") -> Dict[str, Any]:
    """
    設定ファイルを読み込む
    
    Parameters
    ----------
    config_name : str
        設定ファイル名（拡張子なし）。デフォルトは "project"
    
    Returns
    -------
    Dict[str, Any]
        設定の辞書
    """
    config_path = REPO_ROOT / "config" / f"{config_name}.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"設定ファイルが見つかりません: {config_path}")
    
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    
    return config


def load_paths() -> Dict[str, Path]:
    """
    paths.yaml を読み込んで、パスを Path オブジェクトとして返す
    
    Returns
    -------
    Dict[str, Path]
        パスの辞書（値は Path オブジェクト、REPO_ROOT からの相対パス）
    """
    paths_config = load_config("paths")
    
    # 相対パスを REPO_ROOT からの絶対パスに変換
    paths = {}
    for key, value in paths_config.items():
        if isinstance(value, str):
            # 相対パスの場合は REPO_ROOT からの相対とみなす
            if Path(value).is_absolute():
                paths[key] = Path(value)
            else:
                paths[key] = REPO_ROOT / value
        else:
            paths[key] = value
    
    return paths

