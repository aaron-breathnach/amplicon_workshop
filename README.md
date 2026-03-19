# 🧬 Amplicon Sequencing Analysis Workshop
### DADA2 + phyloseq (R-based workflow)

This repository contains materials for a 1-day workshop on amplicon sequencing analysis in R, using **DADA2** for ASV inference and **phyloseq** for downstream analysis.

## 📖 Overview

By the end of this workshop, participants should understand:

- what **ASVs (Amplicon Sequence Variants)** are
- how raw amplicon reads are processed
- how to generate an **ASV table**
- how to assign taxonomy
- how to perform basic diversity analysis in R

## 🔬 Workflow overview

```text
FASTQ reads
   ↓
Quality inspection
   ↓
Filtering and trimming
   ↓
Denoising (DADA2)
   ↓
ASV table
   ↓
Taxonomy assignment
   ↓
phyloseq analysis
   ↓
Diversity + visualisation
```

## 🧠 Key concepts

### What is an amplicon?
A short DNA region amplified using PCR, such as the 16S rRNA gene.

### What is an ASV?
An **ASV (Amplicon Sequence Variant)** is an exact DNA sequence inferred from the data.

Compared with OTUs, ASVs are:

- higher resolution
- reproducible across studies
- not dependent on arbitrary clustering thresholds

### Why DADA2?
DADA2 models sequencing errors directly and infers true biological sequences rather than clustering reads into OTUs.

## ⚙️ Workshop pipeline

### 1. Inspect read quality
```r
plotQualityProfile()
```

### 2. Filter and trim reads
```r
filterAndTrim()
```

Key parameters:

- `truncLen`: where to truncate reads
- `maxEE`: expected error threshold

### 3. Learn error rates
```r
learnErrors()
```

### 4. Infer ASVs
```r
dada()
```

### 5. Merge paired reads
```r
mergePairs()
```

### 6. Remove chimeras
```r
removeBimeraDenovo()
```

### 7. Create ASV table
```r
makeSequenceTable()
```

### 8. Assign taxonomy
```r
assignTaxonomy()
```

## 📦 What is a phyloseq object?
A **phyloseq object** combines:

- ASV count table
- taxonomy table
- sample metadata

```r
ps <- phyloseq(...)
```

## 📊 Basic analyses

### Alpha diversity
```r
plot_richness(ps)
```

### Beta diversity
```r
ordinate(ps, method = "PCoA", distance = "bray")
```

### Taxonomic composition
Relative abundance bar plots can be used to visualise dominant taxa.

## ⚠️ Common pitfalls

- poor-quality reads
- insufficient overlap between paired reads
- over-interpreting relative abundance
- contamination in low-biomass samples
- trying to do too much in one workflow

## ⚙️ Typical parameter choices

```r
truncLen = c(240, 200)
maxEE = c(2, 2)
```

These are sensible starting points for many Illumina 16S datasets, but they should be checked against the actual quality profiles and amplicon length.

## ✅ Repository contents

- `scripts/amplicon_dada2_workflow.R` — full annotated workshop script
- `docs/participant_handout.md` — participant handout
- `docs/setup.md` — pre-workshop installation/setup instructions
- `docs/troubleshooting.md` — common problems and fixes
- `data/metadata.csv` — example metadata template
- `results/` — suggested location for outputs

## ▶️ Getting started

1. Follow the setup instructions in [`docs/setup.md`](docs/setup.md).
2. Place FASTQ files in `data/fastq/`.
3. Add metadata to `data/metadata.csv`.
4. Run [`scripts/amplicon_dada2_workflow.R`](scripts/amplicon_dada2_workflow.R) in RStudio.

## 📁 Suggested repository structure

```text
amplicon_workshop/
├── README.md
├── .gitignore
├── data/
│   ├── fastq/
│   └── metadata.csv
├── docs/
│   ├── participant_handout.md
│   ├── setup.md
│   └── troubleshooting.md
├── reference/
│   └── silva_nr99_v138.2_train_set.fa.gz
├── results/
└── scripts/
    └── amplicon_dada2_workflow.R
```

## 💡 Final thought

Amplicon sequencing is a powerful and accessible method, but good interpretation matters more than complicated analysis.
