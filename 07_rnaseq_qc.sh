#!/usr/bin/env bash
# =============================================================================
# 07_rnaseq_qc.sh — RNA-Seq специфичные метрики
#
# Зачем нужны специальные метрики для RNA-Seq:
#   В отличие от WGS, у RNA-Seq есть:
#   - Сплайсинг (интроны не должны содержать ридов → если много = ДНК-контам.)
#   - Стрэндированность (риды должны идти с определённой цепи гена)
#   - Неравномерное покрытие тела гена (3'-bias = деградация РНК)
#
# Инструменты:
#   RSeQC read_distribution.py → % ридов на экзоны/интроны/UTR
#   RSeQC geneBody_coverage.py → равномерность покрытия 5'→3'
#   Qualimap rnaseq            → комплексный отчёт
# =============================================================================
set -euo pipefail
trap 'printf "[ERROR] Ошибка в строке %s\n" "$LINENO" >&2' ERR

# ── conda ──────────────────────────────────────────────────────────────────
# FIX: conda activate отсутствовал → RSeQC / Qualimap не найдены
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate qc
# ──────────────────────────────────────────────────────────────────────────

DATA="data"
MAP="results/mapping"
QC_MAP="results/qc_mapping"
mkdir -p "$QC_MAP/rnaseq"

YEAST_BAM="$MAP/yeast_dedup.bam"   # физически дедуплицированный BAM

# ─────────────────────────────────────────────────────────────
# 0. Конвертация GTF → BED (требует UCSC инструменты)
# ─────────────────────────────────────────────────────────────
# RSeQC требует BED-формат, а не GTF/GFF
# GTF → GenePred → BED
# FIX: GTF-файл переименован в соответствии с FASTA (GCF_000146045.2_R64_genomic)
GTF="$DATA/GCF_000146045.2_R64_genomic.gtf"
BED="$DATA/GCF_000146045.2_R64_genomic.bed"

if [ -f "$GTF" ] && [ ! -f "$BED" ]; then
    printf "Конвертация аннотации GTF → BED...\n"
    gtfToGenePred \
        -genePredExt \
        -ignoreGroupsWithoutExons \
        "$GTF" \
        "${GTF%.gtf}.genePred"
    genePredToBed \
        "${GTF%.gtf}.genePred" \
        "$BED"
fi

# ─────────────────────────────────────────────────────────────
# 1. Распределение ридов по частям гена
# ─────────────────────────────────────────────────────────────
printf "\n[1/3] read_distribution.py — куда падают риды?\n"
# Ожидание для poly-A RNA-Seq:
#   CDS_Exons:   ~60-70%  ← большинство ридов
#   5'UTR:       ~5-10%
#   3'UTR:       ~15-20%
#   Introns:     < 15%   ← если больше → ДНК-контаминация?
#   Intergenic:  < 5%
if [ -f "$BED" ] && [ -f "$YEAST_BAM" ]; then
    read_distribution.py \
        -r "$BED" \
        -i "$YEAST_BAM" \
        > "$QC_MAP/rnaseq/read_distribution.txt" 2>/dev/null || \
        printf "  RSeQC не установлен, пропускаем\n"
    head -20 "$QC_MAP/rnaseq/read_distribution.txt" 2>/dev/null || true
else
    printf "  BED или YEAST_BAM не найден, пропускаем\n"
fi

# ─────────────────────────────────────────────────────────────
# 2. Покрытие тела гена (5' → 3')
# ─────────────────────────────────────────────────────────────
printf "\n[2/3] geneBody_coverage.py — равномерность 5'→3'?\n"
# Идеально: плоская горизонтальная линия
# 3'-bias (правый конец высокий): деградация РНК или poly-A enrichment
# 5'-bias (левый конец высокий): cap-enriched библиотека
if [ -f "$BED" ] && [ -f "$YEAST_BAM" ]; then
    geneBody_coverage.py \
        -r "$BED" \
        -i "$YEAST_BAM" \
        -o "$QC_MAP/rnaseq/genebody" 2>/dev/null || \
        printf "  RSeQC не установлен, пропускаем\n"
fi

# ─────────────────────────────────────────────────────────────
# 3. Qualimap — комплексный QC для RNA-Seq
# ─────────────────────────────────────────────────────────────
printf "\n[3/3] Qualimap rnaseq — комплексный отчёт...\n"
# Qualimap bamqc: для WGS/WES/ChIP
# Qualimap rnaseq: специально для RNA-Seq (знает про интроны)
if [ -f "$YEAST_BAM" ] && [ -f "$GTF" ]; then
    # FIX: --java-mem-size=4G → 8G (4G вызывает OOM на реальных данных)
    qualimap rnaseq \
        -bam "$YEAST_BAM" \
        -gtf "$GTF" \
        -outdir "$QC_MAP/rnaseq/qualimap" \
        --java-mem-size=8G \
        2>/dev/null || printf "  Qualimap не установлен, пропускаем\n"
fi

# Qualimap для WGS (E. coli):
printf "\nQualimap bamqc для WGS (E. coli)...\n"
if [ -f "$MAP/SRR292770_markdup.bam" ]; then
    # FIX: --java-mem-size=4G → 8G
    qualimap bamqc \
        -bam "$MAP/SRR292770_markdup.bam" \
        -outdir "$QC_MAP/qualimap_ecoli" \
        --java-mem-size=8G \
        -c \
        2>/dev/null || printf "  Qualimap не установлен, пропускаем\n"
fi

printf "\nОтчёты: %s/rnaseq/\n" "$QC_MAP"
