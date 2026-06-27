# Теория: Форматы SAM/BAM

## Содержание
1. [SAM — структура файла](#1-sam--структура-файла)
2. [Флаги SAM (FLAG)](#2-флаги-sam-flag)
3. [CIGAR — кодирование выравнивания](#3-cigar--кодирование-выравнивания)
4. [MAPQ — качество картирования](#4-mapq--качество-картирования)
5. [SAM vs BAM vs CRAM](#5-sam-vs-bam-vs-cram)
6. [Команды samtools](#6-команды-samtools)

---

## 1. SAM — структура файла

SAM (Sequence Alignment/Map) — текстовый формат хранения выравниваний.

### Заголовок (@)

```
@HD  VN:1.6  SO:coordinate        ← версия SAM и тип сортировки
@SQ  SN:U00096.3  LN:4641652      ← имя и длина хромосомы/контига
@RG  ID:SRR292770  SM:Ecoli  PL:ILLUMINA  LB:lib1   ← Read Group
@PG  ID:bwa-mem2  VN:2.2.1  CL:"bwa-mem2 mem ..."   ← программа

Типы заголовков:
  @HD — общая информация о файле
  @SQ — Sequence Dictionary (список хромосом)
  @RG — Read Group (откуда пришли риды, важно для GATK)
  @PG — Program (история обработки)
```

### Строки выравнивания — 11 обязательных полей

```
r001  99  U00096.3  1000  60  101M  =  1200  301  ATGCATGC...  IIIIII...

Поле  Название  Пример       Описание
───────────────────────────────────────────────────────────────────────
1     QNAME     r001         Имя рида
2     FLAG      99           Битовая маска (см. ниже)
3     RNAME     U00096.3     Хромосома/контиг
4     POS       1000         Позиция (1-based, включительно)
5     MAPQ      60           Качество картирования (0-60)
6     CIGAR     101M         Описание выравнивания
7     RNEXT     =            Хромосома пары (= если та же)
8     PNEXT     1200         Позиция пары
9     TLEN      301          Длина вставки (с учётом обоих ридов)
10    SEQ       ATGCATGC...  Последовательность рида
11    QUAL      IIIIII...    Phred-качество каждого нуклеотида

Опциональные поля (TAG:TYPE:VALUE):
  NM:i:2    → число несовпадений с референсом
  AS:i:145  → alignment score (больше = лучше)
  XS:i:100  → alignment score второго лучшего места
  MD:Z:10A5 → позиции мисматчей
  SA:Z:...  → supplementary alignment (химерные риды)
```

---

## 2. Флаги SAM (FLAG)

FLAG — это целое число, каждый бит которого кодирует один признак рида.

### Таблица битов

| Бит (hex) | Бит (dec) | Признак | Описание |
|---|---|---|---|
| 0x001 | 1 | PAIRED | Рид из парной библиотеки |
| 0x002 | 2 | PROPER_PAIR | Оба рида пары картируются корректно |
| 0x004 | 4 | UNMAP | Рид **не** картируется |
| 0x008 | 8 | MUNMAP | Пара **не** картируется |
| 0x010 | 16 | REVERSE | Рид на минус-цепи |
| 0x020 | 32 | MREVERSE | Пара на минус-цепи |
| 0x040 | 64 | READ1 | Первый рид пары (R1) |
| 0x080 | 128 | READ2 | Второй рид пары (R2) |
| 0x100 | 256 | SECONDARY | Secondary alignment |
| 0x200 | 512 | QCFAIL | Не прошёл фильтр QC |
| 0x400 | 1024 | DUP | PCR/оптический дубликат |
| 0x800 | 2048 | SUPPLEMENTARY | Supplementary (химерный) |

### Расшифровка примеров

```
FLAG = 99  → 64+32+2+1 = READ1 + MREVERSE + PROPER_PAIR + PAIRED
             R1, пара на минус-цепи, properly paired

FLAG = 147 → 128+16+2+1 = READ2 + REVERSE + PROPER_PAIR + PAIRED
             R2, сам на минус-цепи, properly paired

FLAG = 4   → UNMAP: рид не картируется

FLAG = 1024 → DUP: PCR-дубликат (помечен Picard)

Проверить любой флаг:
  samtools flags 99
  → Outputs binary and description of each bit
```

### Фильтрация по флагам

```bash
# Флаг -f: оставить риды У КОТОРЫХ стоят эти биты
# Флаг -F: удалить риды У КОТОРЫХ стоят эти биты

# Только картированные риды:
samtools view -F 4 file.bam

# Только properly paired:
samtools view -f 2 file.bam

# Удалить дубликаты + secondary + unmapped:
samtools view -F 1796 file.bam
# 1796 = 4 (unmapped) + 256 (secondary) + 1024 (duplicate) + 512 (qcfail)

# Только первый рид пары (R1):
samtools view -f 64 file.bam

# Риды на плюс-цепи (нет бита REVERSE):
samtools view -F 16 file.bam
```

---

## 3. CIGAR — кодирование выравнивания

CIGAR (Compact Idiosyncratic Gapped Alignment Report) — компактное описание того, как рид отличается от референса.

### Операторы CIGAR

| Оператор | Значение | Нуклеотиды в риде? | Нуклеотиды в референсе? |
|---|---|---|---|
| M | Match или Mismatch | ✓ | ✓ |
| = | Точное совпадение | ✓ | ✓ |
| X | Несовпадение (мисматч) | ✓ | ✓ |
| I | Вставка в риде | ✓ | ✗ |
| D | Делеция в риде | ✗ | ✓ |
| N | Пропуск в референсе (интрон) | ✗ | ✓ |
| S | Soft clipping (нуклеотиды сохранены) | ✓ | ✗ |
| H | Hard clipping (нуклеотиды выброшены) | ✗ | ✗ |

### Разбор примеров

```
101M
  → 101 нуклеотид выровнен (обычный WGS рид без вариантов)

50M2I49M
  → 50 совпадений, 2 нуклеотида вставки в риде, 49 совпадений
  → Рид на 2 п.н. длиннее референса в этом месте

10S91M
  → 10 п.н. soft-clipped (адаптер?), 91 п.н. выровнено
  → Нуклеотиды первых 10 п.н. сохранены в SEQ, но не выровнены

45M1000N56M
  → 45 совпадений, интрон 1000 п.н., 56 совпадений
  → Типичный RNA-Seq рид, перекрывающий сплайс-сайт

6H5M
  → 6 п.н. hard-clipped (не сохранены в SEQ), 5 выровнено
  → Используется для supplementary alignments

Как читать длину вставки из CIGAR:
  Python:
    import re
    cigar = "50M2I49M"
    insertions = sum(int(n) for n,op in re.findall(r'(\d+)([MIDNSHPX=])', cigar) if op == 'I')
```

### Позиция конца рида

```
Для рида с POS=1000 и CIGAR=50M2I49M:
  Длина рида в референсных координатах = M + D + N − I
  = 50 + 49 = 99 п.н. (вставки не занимают место в референсе)
  Конец = 1000 + 99 − 1 = 1098
```

---

## 4. MAPQ — качество картирования

```
MAPQ = −10 × log₁₀( P(неверное картирование) )

Аналогично Phred-качеству:
  MAPQ = 10 → P(ошибка) = 10%   → 1 рид из 10 картируется неверно
  MAPQ = 20 → P(ошибка) = 1%    → стандартный фильтр
  MAPQ = 30 → P(ошибка) = 0.1%
  MAPQ = 60 → P(ошибка) = 0.000001% — максимум в BWA-MEM2
```

### Как вычисляется MAPQ

```
MAPQ в BWA-MEM2 зависит от:
  - AS  (alignment score лучшего выравнивания)
  - XS  (alignment score второго лучшего выравнивания)

Если AS >> XS → рид уникально картируется → MAPQ высокий (60)
Если AS ≈ XS  → два одинаково хороших места → MAPQ = 0

Пример:
  Рид в уникальном регионе:    AS=145, XS=60   → MAPQ=60
  Рид в дублированном гене:    AS=145, XS=142  → MAPQ=3
  Рид в повторах (LINE, SINE): AS=145, XS=145  → MAPQ=0
```

### Мультикартирование и стратегии работы

```
Рид картируется в нескольких местах → несколько строк в SAM:
  - Primary alignment  (FLAG & 0x100 = 0) — основное, лучшее
  - Secondary alignment (FLAG & 0x100 = 1) — альтернативные

Стратегии:
  RNA-Seq (подсчёт экспрессии):
    → Удалять все мультикартированные: samtools view -q 20
    → Иначе один ген получит двойное/тройное покрытие

  WGS вариант-коллинг:
    → Оставлять primary, удалять secondary: samtools view -F 256
    → Иначе пропустим варианты в паралогах

  ChIP-Seq/CUT&RUN:
    → Зависит от цели: пики vs уникальные связывания
```

---

## 5. SAM vs BAM vs CRAM

| Формат | Тип | Размер (WGS 30×) | Требует референс | Команда создания |
|---|---|---|---|---|
| SAM | Текст | ~200 GB | нет | `bwa-mem2 mem ... > out.sam` |
| BAM | Бинарный | ~45 GB | нет | `samtools sort ... -o out.bam` |
| CRAM | Бинарный | ~15 GB | да | `samtools view -C -T ref.fa out.bam` |

```
BAM = BGZF (Blocked GNU Zip Format) компрессия SAM
  → Можно делать произвольный доступ с .bai индексом
  → samtools index создаёт .bai (для BAM) или .crai (для CRAM)

CRAM хранит только ОТЛИЧИЯ от референса → нужен ref при декодировании
  samtools view --reference ref.fasta out.cram
```

---

## 6. Команды samtools

### Основные операции

```bash
# Просмотр BAM в текстовом виде
samtools view file.bam | head -5
samtools view -H file.bam          # только заголовок

# Сортировка (по координатам — стандарт)
samtools sort -@ 8 -m 2G -o sorted.bam input.bam

# Сортировка по имени рида (для парных ридов)
samtools sort -n -o namesort.bam input.bam

# Индексирование
samtools index sorted.bam          # создаёт sorted.bam.bai

# Быстрый просмотр в терминале
samtools tview sorted.bam ref.fasta
# Управление: g → перейти к позиции, ? → помощь, q → выйти

# FASTA-индекс (нужен для IGV и samtools faidx)
samtools faidx ref.fasta
samtools faidx ref.fasta "chr1:1000-2000"   # извлечь регион
```

### Фильтрация и подмножества

```bash
# Конкретный регион (нужен .bai индекс)
samtools view file.bam "U00096.3:1000-5000"

# Только картированные, unique, properly paired
samtools view -b -F 1804 -q 20 file.bam > filtered.bam
# 1804 = 4 (unmapped) + 8 (mate unmapped) + 256 (secondary) + 1024 (dup) + 512 (qcfail)

# Конвертация BAM → FASTQ (обратно для перекартирования)
samtools fastq -1 R1.fq.gz -2 R2.fq.gz \
  -0 unpaired.fq.gz \
  -s singleton.fq.gz \
  sorted.bam

# SAM → BAM одной командой (pipe от выравнивателя)
bwa-mem2 mem ref R1.fq R2.fq | samtools sort -o out.bam
```

### Статистика

```bash
# Общая статистика (много метрик)
samtools flagstat file.bam

# Детальная статистика (для MultiQC)
samtools stats file.bam > file.stats

# Покрытие по контигам
samtools coverage file.bam

# Глубина покрытия по позициям
samtools depth -a file.bam > depth.txt
awk '{sum+=$3; n++} END {print "Mean:", sum/n}' depth.txt
```
