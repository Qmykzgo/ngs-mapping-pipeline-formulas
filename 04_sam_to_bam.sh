#!/usr/bin/env bash
# =============================================================================
# 04_sam_to_bam.sh — SAM → BAM + сортировка + индекс
#
# Зачем конвертировать:
#   SAM (~200 GB для WGS 30×) → BAM (~45 GB, в 4-5× меньше)
#   BAM = BGZF-сжатый SAM (Blocked GNU Zip Format)
#   BGZF позволяет произвольный доступ к позиции — нужен .bai индекс
#
# Зачем сортировать по координатам:
#   samtools tview, IGV, GATK, bcftools требуют отсортированный BAM
#   .bai индекс позволяет перейти к любой позиции за O(log N)
#
# Альтернатива: pipe от выравнивателя (экономит место на диске):
#   bwa-mem2 mem ... | samtools sort -o out.bam
# =============================================================================
set -euo pipefail
trap 'printf "[ERROR] Ошибка в строке %s\n" "$LINENO" >&2' ERR

# ── conda ──────────────────────────────────────────────────────────────────
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate mapping
# ──────────────────────────────────────────────────────────────────────────

# FIX: DATA не был объявлен → unbound variable при set -euo pipefail
DATA="data"
MAP="results/mapping"
THREADS=4
SAMPLE="SRR292770"

printf "[STEP 1/2] SAM → BAM + сортировка по координатам...\n"
# -@ : число потоков
# -m : память на поток (для сортировки)
# -o : выходной файл
samtools sort \
    -@ "$THREADS" \
    -m 2G \
    -o "$MAP/${SAMPLE}_sorted.bam" \
    "$MAP/${SAMPLE}.sam"

printf "[STEP 2/2] Индексирование BAM (.bai)...\n"
# .bai — бинарный индекс, позволяет jump к хромосоме/позиции
samtools index -@ "$THREADS" "$MAP/${SAMPLE}_sorted.bam"

printf "\nУдаляем SAM (экономия места)...\n"
rm "$MAP/${SAMPLE}.sam"

printf "\nРезультаты:\n"
ls -lh "$MAP/${SAMPLE}_sorted.bam"*

printf "\nБазовый просмотр:\n"
samtools view "$MAP/${SAMPLE}_sorted.bam" | head -3

printf "\nВизуализация в терминале:\n"
printf "  samtools tview %s %s\n" \
    "$MAP/${SAMPLE}_sorted.bam" \
    "$DATA/U00096.fasta"
printf "  g → перейти к позиции, ? → помощь, q → выйти\n"

# ─────────────────────────────────────────────────────────────
# Yeast RNA-Seq BAM (если HISAT2 отработал в 03_align.sh)
# ─────────────────────────────────────────────────────────────
if [ -f "$MAP/yeast.sam" ]; then
    printf "\n[YEAST] SAM → BAM + сортировка (S. cerevisiae)...\n"
    samtools sort \
        -@ "$THREADS" \
        -m 2G \
        -o "$MAP/yeast_sorted.bam" \
        "$MAP/yeast.sam"
    samtools index -@ "$THREADS" "$MAP/yeast_sorted.bam"
    rm "$MAP/yeast.sam"
    printf "Yeast BAM: %s\n" "$(ls -lh "$MAP/yeast_sorted.bam" | awk '{print $5, $9}')"
fi
