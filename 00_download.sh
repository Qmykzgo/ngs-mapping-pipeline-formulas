#!/usr/bin/env bash
# =============================================================================
# 00_download.sh — Скачать референс и риды
# =============================================================================
set -euo pipefail
DATA="data"; mkdir -p "$DATA"

echo "[1/3] Референс E. coli K-12 MG1655 (U00096.3)..."
wget -q -O "$DATA/U00096.fasta" \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=U00096.3&rettype=fasta&retmode=text"

echo "[2/3] Риды SRR292770 (Illumina paired-end)..."
wget -q -O "$DATA/SRR292770_1.fastq.gz" \
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292770/SRR292770_1.fastq.gz"
wget -q -O "$DATA/SRR292770_2.fastq.gz" \
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292770/SRR292770_2.fastq.gz"

echo "[3/3] Проверяем..."
ls -lh "$DATA/"

# Расчёт ожидаемого покрытия: C = (N × L) / G
echo ""
echo "Формула покрытия: C = (N × L) / G"
N=$(zcat "$DATA/SRR292770_1.fastq.gz" | awk 'NR%4==1{c++}END{print c*2}')
echo "  N (ридов, paired) = $N"
echo "  L (длина рида)    ≈ 100 п.н."
echo "  G (геном E. coli) = 4,641,652 п.н."
echo "  C ≈ $(echo "$N * 100 / 4641652" | bc)×"
