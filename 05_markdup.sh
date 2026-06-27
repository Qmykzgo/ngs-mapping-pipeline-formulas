#!/usr/bin/env bash
# =============================================================================
# 05_markdup.sh — Маркировка PCR-дубликатов (Picard MarkDuplicates)
#
# Типы дубликатов:
#   PCR-дубликаты:     несколько копий одной молекулы из-за амплификации
#   Оптические дубл.:  один кластер flow cell принят за два (NovaSeq чаще)
#
# Как Picard находит дубликаты:
#   Группирует риды с одинаковыми 5'-концами (обоих ридов пары)
#   В каждой группе: лучший по суммарному качеству = оригинал
#   Остальные → FLAG 0x400 (DUP)
#
# МАРКИРОВКА (не удаление):
#   Риды помечаются, но остаются в файле
#   Большинство инструментов (GATK, samtools mpileup) их автоматически игнорируют
#   Физическое удаление: samtools view -F 1024 если нужно
# =============================================================================
set -euo pipefail
trap 'printf "[ERROR] Ошибка в строке %s\n" "$LINENO" >&2' ERR

# ── conda ──────────────────────────────────────────────────────────────────
# FIX: conda activate отсутствовал → picard/samtools не найдены
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate mapping   # picard теперь в envs/mapping.yml
# ──────────────────────────────────────────────────────────────────────────

MAP="results/mapping"
SAMPLE="SRR292770"

printf "[STEP 1/1] Маркировка PCR-дубликатов (Picard MarkDuplicates)...\n"
picard MarkDuplicates \
    I="$MAP/${SAMPLE}_sorted.bam" \
    O="$MAP/${SAMPLE}_markdup.bam" \
    M="$MAP/${SAMPLE}_dup_metrics.txt" \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    CREATE_INDEX=true
# OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 для NovaSeq (patterned flow cells)
# Для HiSeq/MiSeq используй значение по умолчанию 100

printf "\nСтатистика дубликатов:\n"
grep -A 8 "^LIBRARY" "$MAP/${SAMPLE}_dup_metrics.txt" | head -5
printf "\n  Норма:\n"
printf "    WGS/WES: PERCENT_DUPLICATION < 0.20 (< 20%%)\n"
printf "    RNA-Seq: PERCENT_DUPLICATION < 0.50 (< 50%%)\n"

# ─────────────────────────────────────────────────────────────
# RNA-Seq (yeast): физическое удаление для RSeQC
# ─────────────────────────────────────────────────────────────
if [ -f "$MAP/yeast_sorted.bam" ]; then
    printf "\nRNA-Seq (yeast): маркировка + физическое удаление для RSeQC...\n"
    picard MarkDuplicates \
        I="$MAP/yeast_sorted.bam" \
        O="$MAP/yeast_markdup.bam" \
        M="$MAP/yeast_dup_metrics.txt"
    samtools index "$MAP/yeast_markdup.bam"

    # -F 1024 = исключить FLAG DUP (0x400)
    samtools view -h -b -F 1024 "$MAP/yeast_markdup.bam" \
        > "$MAP/yeast_dedup.bam"
    samtools index "$MAP/yeast_dedup.bam"
    printf "Дедуплицированный BAM: %s\n" "$MAP/yeast_dedup.bam"
else
    printf "[INFO] yeast_sorted.bam не найден, yeast MarkDuplicates пропущен.\n"
    printf "       Сначала запусти bash scripts/03_align.sh && bash scripts/04_sam_to_bam.sh\n"
fi
