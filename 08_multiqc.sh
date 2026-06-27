#!/usr/bin/env bash
# =============================================================================
# 08_multiqc.sh — Агрегировать все QC-отчёты в один HTML
#
# MultiQC автоматически находит и парсит:
#   FastQC отчёты (*_fastqc.zip)
#   samtools stats (*_stats.txt)
#   samtools flagstat (*_flagstat.txt)
#   Picard MarkDuplicates (*_dup_metrics.txt)
#   Picard InsertSizeMetrics (*_insert_size.txt)
#   Picard GcBiasMetrics (*_gc_bias.txt)
#   mosdepth (*.mosdepth.summary.txt)
#   preseq (*_complexity_curve.txt)
#   HISAT2 summary (*hisat2*.log)
#   Qualimap (qualimap_bamqc/)
#
# Один HTML → легко поделиться с коллегами и в публикации
# =============================================================================
set -euo pipefail
trap 'printf "[ERROR] Ошибка в строке %s\n" "$LINENO" >&2' ERR

# ── conda ──────────────────────────────────────────────────────────────────
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate qc
# ──────────────────────────────────────────────────────────────────────────

REPORT_DIR="results/multiqc_report"
mkdir -p "$REPORT_DIR"

printf "Сбор всех QC-отчётов в один MultiQC HTML...\n"

# FIX: убраны inline-комментарии после \ — bash НЕ поддерживает "\ # comment"
# При "command \ # comment" continuation не работает: # начинает комментарий
# до newline, backslash не является последним символом строки
multiqc \
    results/ \
    --outdir "$REPORT_DIR" \
    --filename "multiqc_report.html" \
    --title "NGS Mapping QC — E. coli + S. cerevisiae" \
    --force \
    --verbose

printf "\nОтчёт: %s/multiqc_report.html\n" "$REPORT_DIR"
printf "\nКлючевые разделы в отчёте:\n"
printf "  General Statistics — сводная таблица всех образцов\n"
printf "  FastQC            — качество сырых ридов\n"
printf "  Picard            — дубликаты, insert size, GC-bias\n"
printf "  samtools          — %% картировалось, покрытие\n"
printf "  mosdepth          — распределение покрытия\n"
printf "  preseq            — кривая насыщения\n"

# Создать архив для отправки:
zip -r "$REPORT_DIR/multiqc_report.zip" "$REPORT_DIR/" 2>/dev/null
printf "\nАрхив: %s/multiqc_report.zip\n" "$REPORT_DIR"
