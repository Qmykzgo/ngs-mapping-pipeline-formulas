#!/usr/bin/env bash
# =============================================================================
# 01_qc_raw.sh — FastQC для сырых ридов
#
# Что смотрим в FastQC-отчёте:
#   Per base sequence quality → Q30 на большинстве позиций?
#   Per sequence GC content   → соответствует организму (E.coli ~51%)?
#   Adapter Content           → есть адаптеры → нужен тримминг
#   Overrepresented sequences → контаминация / rRNA?
#   Sequence Duplication      → % дубликатов (предварительно)
# =============================================================================
set -euo pipefail
conda activate QC_fastq

DATA="data"; QC="results/qc"; THREADS=4
mkdir -p "$QC"

echo "FastQC на сырых ридах..."
fastqc --threads "$THREADS" --outdir "$QC" \
  "$DATA/SRR292770_1.fastq.gz" "$DATA/SRR292770_2.fastq.gz"

echo "Готово: $QC/"
echo "Открой HTML-отчёт в браузере."
