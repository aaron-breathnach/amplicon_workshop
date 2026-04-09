# 📚 Useful Resources

This page contains useful links and resources for learning more about amplicon sequencing and microbiome analysis.

---

## 🔬 DADA2

- Official tutorial: https://benjjneb.github.io/dada2/tutorial.html  
- DADA2 documentation: https://benjjneb.github.io/dada2/  
- Bioconductor package page: https://bioconductor.org/packages/dada2/  

---

## 📦 phyloseq

- phyloseq tutorial: https://joey711.github.io/phyloseq/  
- Bioconductor page: https://bioconductor.org/packages/phyloseq/  

---

## 🧬 Reference databases

- SILVA: https://www.arb-silva.de/  
- GTDB: https://gtdb.ecogenomic.org/  
- UNITE (fungal ITS): https://unite.ut.ee/  

---

## 🦠 Microbiome analysis tutorials

- R microbiome tutorials: https://microbiome.github.io/tutorials/  
- QIIME2 docs (useful for concepts): https://docs.qiime2.org/  

---

## 📖 Key papers

- DADA2: https://doi.org/10.1038/nmeth.3869  
- phyloseq: https://doi.org/10.1371/journal.pone.0061217  

---

## 🧠 General bioinformatics learning

- Bioconductor workflows: https://bioconductor.org/help/workflows/  
- Carpentries lessons: https://carpentries.org/  

---

## 📊 Statistical analysis

- Alpha diversity metrics (e.g. Shannon, Simpson) can be compared between groups using standard statistical tests (e.g. Wilcoxon or linear models).
- PERMANOVA (e.g. via `adonis2` in vegan) can be used to test for associations between metadata variables and overall microbiome composition (beta diversity).
  - Ordination methods (e.g. PCoA, NMDS) are commonly used to visualise differences in community structure.
- MaAsLin 3 can be used to identify differentially abundant microbiome features while accounting for covariates.
  - Alternative approaches that account for compositionality include ANCOM-BC and ALDEx2.

---

## 🤔 Rarefaction

Rarefaction is a controversial topic in microbiome research. I am not going to advise people on what strategy to use, but here are the arguments for and against rarefaction:

- For: https://journals.asm.org/doi/full/10.1128/msphere.00354-23
- Against: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531

---

## ⚠️ Important note

Many microbiome tools differ in implementation, but the core concepts are shared:

- Quality filtering  
- Error modelling / denoising  
- Taxonomic assignment  
- Diversity analysis  

Focus on understanding the concepts rather than memorising tools.

---
