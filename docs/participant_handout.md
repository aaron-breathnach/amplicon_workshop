# Participant Handout

## What are we doing today?

In this workshop, we will learn how to analyse amplicon sequencing data in R using **DADA2** and **phyloseq**.

By the end, you should understand:

- what an **ASV** is
- how reads are processed into an ASV table
- how taxonomy is assigned
- how to perform basic diversity analysis

## Big-picture workflow

```text
FASTQ reads
   ↓
Quality filtering
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

## Key concepts

### Amplicon
A short DNA region amplified by PCR, such as 16S rRNA.

### ASV
An **Amplicon Sequence Variant** is an exact inferred biological sequence.

### Why ASVs instead of OTUs?
ASVs avoid arbitrary clustering thresholds and provide higher resolution.

## Main DADA2 steps

1. Inspect quality profiles
2. Filter and trim reads
3. Learn error rates
4. Infer ASVs
5. Merge paired reads
6. Remove chimeras
7. Assign taxonomy

## Common pitfalls

- low-quality reads
- poor overlap between forward and reverse reads
- contamination
- over-interpreting relative abundance

## Key takeaway

The most important idea is that DADA2 **models sequencing errors** and **infers exact sequence variants**.
