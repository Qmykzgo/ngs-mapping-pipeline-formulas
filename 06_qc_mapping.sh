#!/usr/bin/env bash
# =============================================================================
# 06_qc_mapping.sh — Метрики качества выравнивания
#
# Инструменты:
#   samtools stats    — общая статистика (% картировалось, % дубликатов)
#   mosdepth          — покрытие (быстро, рекомендуется вместо samtools depth)
#   preseq            — кривая насыщения библиотеки
#   Picard InsertSize — распределение длин вставок
#   Picard GcBias     — GC-смещение покрытия
# =============================================================================
set -euo pipefail
trap 'printf "[ERROR] Ошибка в строке %s\n" "$LINENO" >&2' ERR

# ── conda ──────────────────────────────────────────────────────────────────
# FIX: conda activate отсутствовал → все инструменты не найдены
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate mapping
# ──────────────────────────────────────────────────────────────────────────

DATA="data"
MAP="results/mapping"
QC_MAP="results/qc_mapping"
THREADS=4
mkdir -p "$QC_MAP"
SAMPLE="SRR292770"
BAM="$MAP/${SAMPLE}_markdup.bam"

printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
printf "[1/5] samtools flagstat — базовые флаги\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
# Показывает: % картировалось, % properly paired, % дубликатов
samtools flagstat "$BAM" | tee "$QC_MAP/${SAMPLE}_flagstat.txt"

printf "\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
printf "[2/5] samtools stats — детальная статистика\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
# Многострочный отчёт. Ключевые строки (SN = Summary Numbers):
#   SN reads mapped:             → число картированных ридов
#   SN non-primary alignments:   → secondary alignments (мультикартирование)
#   SN average length:           → средняя длина рида
#   SN average quality:          → среднее качество
samtools stats "$BAM" > "$QC_MAP/${SAMPLE}_stats.txt"
grep "^SN" "$QC_MAP/${SAMPLE}_stats.txt" | \
    grep -E "reads (mapped|total)|non-primary|average (length|quality)" | \
    column -t

printf "\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
printf "[3/5] mosdepth — покрытие\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
# Формула покрытия: C = (N × L) / G
# mosdepth считает реальное покрытие по BAM (учитывает CIGAR)
# -n/--no-abbrev : не записывать per-base файл (намного быстрее для больших геномов)
# -F 1024        : исключить дубликаты из подсчёта покрытия
mosdepth \
    --no-abbrev \
    --fast-mode \
    --flag 1024 \
    "$QC_MAP/${SAMPLE}" \
    "$BAM"

printf "\nРезультат покрытия:\n"
column -t "$QC_MAP/${SAMPLE}.mosdepth.summary.txt"
printf "\nИдеально для E. coli WGS: среднее покрытие > 50×\n"

printf "\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
printf "[4/5] preseq — кривая насыщения библиотеки\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
# Оценивает: если секвенировать больше, сколько получим новых уникальных молекул?
# Плато → библиотека истощена, досеквенирование бессмысленно
# Линейный рост → библиотека богатая, можно секвенировать больше
#
# FIX: preseq lc_extrap -pe требует name-sorted BAM (не coord-sorted!)
#      Делаем временную name-sorted копию, затем удаляем
BAM_NS="${QC_MAP}/${SAMPLE}_namesort_tmp.bam"
printf "  Временная name-sorted копия для preseq...\n"
samtools sort -n -@ "$THREADS" -m 2G "$BAM" -o "$BAM_NS"

preseq lc_extrap \
    -pe \
    -B "$BAM_NS" \
    -o "$QC_MAP/${SAMPLE}_complexity_curve.txt" \
    2>"$QC_MAP/${SAMPLE}_preseq.log" || \
    printf "[WARN] preseq ошибка или не установлен, см. %s\n" \
           "$QC_MAP/${SAMPLE}_preseq.log"

rm -f "$BAM_NS"

printf "\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
printf "[5/5] Picard метрики\n"
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"

# Insert size distribution
# Ожидание для WGS Illumina: медиана 150-400 п.н.
# Если < 100: библиотека коротких фрагментов (деградация?)
# Если > 1000 и ориентация RF: возможно mate-pair библиотека
picard CollectInsertSizeMetrics \
    I="$BAM" \
    O="$QC_MAP/${SAMPLE}_insert_size.txt" \
    H="$QC_MAP/${SAMPLE}_insert_size.pdf" \
    2>/dev/null

printf "\nInsert size метрики:\n"
grep -A 2 "^MEDIAN_INSERT" "$QC_MAP/${SAMPLE}_insert_size.txt" | head -3

# GC bias
# NORMALIZED_COVERAGE ≈ 1.0 для всех GC → хорошо
# NORMALIZED_COVERAGE < 0.5 при высоком GC → ПЦР-проблема
picard CollectGcBiasMetrics \
    I="$BAM" \
    O="$QC_MAP/${SAMPLE}_gc_bias.txt" \
    CHART="$QC_MAP/${SAMPLE}_gc_bias.pdf" \
    S="$QC_MAP/${SAMPLE}_gc_summary.txt" \
    R="$DATA/U00096.fasta" \
    2>/dev/null

printf "\nВсе метрики: %s/\n" "$QC_MAP"
printf "Следующий шаг: bash scripts/08_multiqc.sh\n"
