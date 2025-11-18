#!/usr/bin/env bash
set -euo pipefail

# リポジトリルートを決定
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MMSEQS_BIN="${MMSEQS_BIN:-mmseqs}"

HUMAN_DB="${REPO_ROOT}/data/interim/mmseqs/human_db"
PDB_CHAIN_DB="${REPO_ROOT}/data/interim/mmseqs/pdb_chain_db"
RESULT_DB="${REPO_ROOT}/data/interim/mmseqs/human_vs_pdb"
RESULT_TSV="${REPO_ROOT}/data/interim/mmseqs/human_vs_pdb.m8"

echo "[INFO] REPO_ROOT     : ${REPO_ROOT}"
echo "[INFO] HUMAN_DB      : ${HUMAN_DB}"
echo "[INFO] PDB_CHAIN_DB  : ${PDB_CHAIN_DB}"
echo "[INFO] RESULT_DB     : ${RESULT_DB}"
echo "[INFO] RESULT_TSV    : ${RESULT_TSV}"
echo "[INFO] MMSEQS_BIN    : ${MMSEQS_BIN}"

# 入力の存在チェック
if [[ ! -f "${HUMAN_DB}.dbtype" ]]; then
  echo "[ERROR] human_db mmseqs DB not found: ${HUMAN_DB}.dbtype" >&2
  exit 1
fi

if [[ ! -f "${PDB_CHAIN_DB}.dbtype" ]]; then
  echo "[ERROR] pdb_chain_db mmseqs DB not found: ${PDB_CHAIN_DB}.dbtype" >&2
  exit 1
fi

if [[ ! -f "${RESULT_DB}.dbtype" ]]; then
  echo "[ERROR] result DB (human_vs_pdb) not found: ${RESULT_DB}.dbtype" >&2
  echo "[ERROR] 先に run_mmseqs_search.sh を実行してください。" >&2
  exit 1
fi

# すでに TSV があるなら再計算しない
if [[ -f "${RESULT_TSV}" ]]; then
  echo "[INFO] Result TSV already exists. Skip convertalis."
  echo "[INFO] If you want to re-generate it, delete:"
  echo "       ${RESULT_TSV}"
  exit 0
fi

# TSV 変換
echo "[INFO] Running mmseqs convertalis ..."
"${MMSEQS_BIN}" convertalis \
  "${HUMAN_DB}" \
  "${PDB_CHAIN_DB}" \
  "${RESULT_DB}" \
  "${RESULT_TSV}" \
  --format-output "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits"

echo "[INFO] Done. TSV written to: ${RESULT_TSV}"
