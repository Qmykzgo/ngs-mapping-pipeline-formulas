# 🧬 NGS Read Mapping Pipeline

> Полный учебный репозиторий по картированию NGS-данных на референсный геном.  
> Тест-кейс: *E. coli* K-12 MG1655 (WGS) и *S. cerevisiae* (RNA-Seq).

---

## 📁 Структура репозитория

```
ngs-mapping-pipeline/
├── README.md                        ← вы здесь
├── envs/
│   ├── mapping.yml                  ← conda: выравнивание (bwa-mem2, samtools, hisat2)
│   └── qc.yml                       ← conda: контроль качества (fastqc, picard, mosdepth)
├── scripts/
│   ├── 00_download.sh               ← скачать референс и риды
│   ├── 01_qc_raw.sh                 ← FastQC сырых ридов
│   ├── 02_index.sh                  ← построить FM-индекс референса
│   ├── 03_align.sh                  ← картирование bwa-mem2 / hisat2
│   ├── 04_sam_to_bam.sh             ← SAM → BAM, сортировка, индекс
│   ├── 05_markdup.sh                ← маркировка PCR-дубликатов (Picard)
│   ├── 06_qc_mapping.sh             ← метрики выравнивания (samtools stats, mosdepth, picard)
│   ├── 07_rnaseq_qc.sh              ← RNA-Seq специфика (RSeQC, Qualimap)
│   └── 08_multiqc.sh                ← объединить все отчёты в один HTML
├── docs/
│   ├── theory_algorithms.md         ← BWT/FM-index, алгоритмы поиска
│   ├── theory_sam_format.md         ← SAM/BAM формат: флаги, CIGAR, MAPQ
│   ├── theory_qc_metrics.md         ← метрики качества: дубликаты, покрытие, insert size
│   └── formulas.md                  ← все формулы в одном месте
└── results/                         ← gitignore'd, создаётся при запуске
    ├── qc/
    ├── mapping/
    └── qc_mapping/
```

---

## 🚀 Быстрый старт

```bash
# 1. Клонировать репозиторий
git clone https://github.com/yourname/ngs-mapping-pipeline.git
cd ngs-mapping-pipeline

# 2. Создать conda-окружения
conda env create -f envs/mapping.yml
conda env create -f envs/qc.yml

# 3. Запустить пайплайн (E. coli WGS, шаг за шагом)
bash scripts/00_download.sh
bash scripts/01_qc_raw.sh
bash scripts/02_index.sh
bash scripts/03_align.sh
bash scripts/04_sam_to_bam.sh
bash scripts/05_markdup.sh
bash scripts/06_qc_mapping.sh
bash scripts/08_multiqc.sh
```

---

## 🗺 Схема пайплайна

```
Сырые риды (.fastq.gz)
        │
        ▼
[01] FastQC ──────────────────────────────► QC-отчёт (HTML)
        │
        ▼
[02] Построение FM-индекса (bwa-mem2 index / hisat2-build)
        │
        ▼
[03] Картирование (bwa-mem2 mem / hisat2)
        │
        ▼  .sam
[04] SAM → BAM + сортировка + .bai индекс
        │
        ▼  .bam
[05] Маркировка дубликатов (Picard MarkDuplicates)
        │
        ▼  _markdup.bam
[06] Метрики: samtools stats / mosdepth / Picard metrics
        │
        ▼
[08] MultiQC ────────────────────────────► Единый HTML-отчёт
```

---

## 📊 Ключевые формулы

### Формула покрытия (Lander-Waterman)
```
C = (N × L) / G

  C — среднее покрытие (×)
  N — число ридов
  L — длина рида (п.н.)
  G — размер генома (п.н.)

Пример: E. coli, N=5M ридов, L=100 п.н., G=4.64 Мб
  C = (5,000,000 × 100) / 4,640,000 ≈ 108×
```

### MAPQ — качество картирования
```
MAPQ = −10 × log₁₀( P(неверное картирование) )

MAPQ = 0   → рид картируется в ≥2 местах одинаково хорошо
MAPQ = 20  → P(ошибка) = 1%  — минимальный фильтр для вариант-коллинга
MAPQ = 60  → максимум в BWA-MEM2, уникальное вхождение
```

### Вычислительная сложность алгоритмов поиска
```
Алгоритм               │ Сложность     │ RAM (геном человека)
───────────────────────┼───────────────┼─────────────────────
Линейный поиск         │ O(N × M)      │ ~3 GB
Бинарный (suffix array)│ O(M × log N)  │ ~12 GB
Суффиксное дерево      │ O(M)          │ ~64 GB  ← слишком!
BWT + FM-index         │ O(M)          │ ~3-4 GB ← используется
```

---

## 🔬 Тест-данные

| Параметр | E. coli (WGS) | S. cerevisiae (RNA-Seq) |
|---|---|---|
| Референс | U00096.3 (NCBI) | GCF_000146045.2_R64 |
| Размер генома | 4.64 Мб | 12.2 Мб |
| SRA ID | SRR292770 | SRR36689450 |
| Тип данных | paired-end WGS | paired-end RNA-Seq |
| Выравниватель | BWA-MEM2 | HISAT2 |

---

## 🛠 Инструменты

| Инструмент | Назначение | Conda-пакет |
|---|---|---|
| FastQC | QC сырых ридов | `fastqc` |
| MultiQC | Агрегация всех отчётов | `multiqc` |
| BWA-MEM2 | Картирование ДНК (BWT) | `bwa-mem2` |
| HISAT2 | Картирование РНК (сплайсинг) | `hisat2` |
| samtools | SAM/BAM операции | `samtools` |
| Picard | Метрики, маркировка дубликатов | `picard` |
| mosdepth | Покрытие (быстро) | `mosdepth` |
| preseq | Кривая сложности библиотеки | `preseq` |
| RSeQC | RNA-Seq специфичные метрики | `rseqc` |
| Qualimap | Комплексный QC BAM | `qualimap` |

---

## 📖 Документация

| Файл | Содержание |
|---|---|
| [docs/theory_algorithms.md](docs/theory_algorithms.md) | BWT, FM-index, seed-and-extend, Нидлман-Вунш |
| [docs/theory_sam_format.md](docs/theory_sam_format.md) | SAM/BAM структура, FLAGS, CIGAR, MAPQ |
| [docs/theory_qc_metrics.md](docs/theory_qc_metrics.md) | Дубликаты, покрытие, insert size, GC-bias |
| [docs/formulas.md](docs/formulas.md) | Все формулы на одной странице |
