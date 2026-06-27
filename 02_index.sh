#!/usr/bin/env bash
# =============================================================================
# 02_index.sh — Построить FM-индекс референса
#
# Что происходит при индексировании:
#   1. Все суффиксы референса сортируются → Suffix Array
#   2. BWT(G) = последние буквы отсортированных ротаций
#   3. FM-index = BWT + таблицы Count[] и Occ[]
#   Результат: поиск рида за O(M) независимо от длины генома
#
# bwa-mem2 vs bwa:
#   bwa-mem2 — та же математика, но в 2-4× быстрее (SIMD/AVX512)
#   Требует чуть больше RAM при построении индекса (~14 GB для человека)
#   Результаты идентичны bwa mem
# =============================================================================
set -euo pipefail
conda activate mapping

DATA="data"

echo "Индексирование референса (FM-index / BWT)..."
echo "Файлы индекса: .amb .ann .bwt.2bit.64 .pac .0123"
echo ""

# -p prefix: имя префикса для файлов индекса
bwa-mem2 index -p "$DATA/ecoli_ref" "$DATA/U00096.fasta"

echo ""
echo "Индекс построен:"
ls -lh "$DATA/ecoli_ref".*

echo ""
echo "FASTA-индекс для samtools faidx и IGV:"
samtools faidx "$DATA/U00096.fasta"
echo "Создан: $DATA/U00096.fasta.fai"

# Для RNA-Seq / HISAT2 (дрожжи):
if [ -f "$DATA/GCF_000146045.2_R64_genomic.fna" ]; then
  echo ""
  echo "Индекс HISAT2 для S. cerevisiae..."
  # hisat2-build создаёт .ht2 файлы
  hisat2-build \
    "$DATA/GCF_000146045.2_R64_genomic.fna" \
    "$DATA/s_cerevisiae/genome"
fi
