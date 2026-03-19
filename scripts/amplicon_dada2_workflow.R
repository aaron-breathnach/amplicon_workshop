############################################################
# Amplicon sequencing analysis workshop
# Pipeline: DADA2 + phyloseq
# Participant-facing teaching script
############################################################

# 0. Install packages (run once only)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Uncomment if needed
# BiocManager::install(c("dada2", "phyloseq", "Biostrings"))
# install.packages(c("ggplot2", "dplyr", "readr", "tibble"))

# 1. Load libraries
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)

# 2. Define input files
path <- "data/fastq"
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz$", full.names = TRUE))

sample.names <- basename(fnFs) |>
  sub("_R1.fastq.gz$", "", x = _)

stopifnot(length(fnFs) == length(fnRs))

sample.names.r <- basename(fnRs) |>
  sub("_R2.fastq.gz$", "", x = _)

stopifnot(all(sample.names == sample.names.r))

# 3. Inspect quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# 4. Filter and trim
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fwd = fnFs,
  filt = filtFs,
  rev = fnRs,
  filt.rev = filtRs,
  truncLen = c(240, 200),
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

out_df <- as.data.frame(out) |>
  rownames_to_column("sample")

# 5. Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# 6. Dereplicate reads
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# 7. Infer ASVs
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# 8. Merge pairs
mergers <- mergePairs(
  dadaF = dadaFs,
  derepF = derepFs,
  dadaR = dadaRs,
  derepR = derepRs
)

# 9. Construct sequence table
seqtab <- makeSequenceTable(mergers)

# 10. Remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE
)

# 11. Track read retention
getN <- function(x) sum(getUniques(x))

track <- cbind(
  input = out[, "reads.in"],
  filtered = out[, "reads.out"],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim)
)

rownames(track) <- sample.names
track_df <- as.data.frame(track) |>
  rownames_to_column("sample")

write_csv(track_df, "results/read_tracking.csv")

# 12. Assign taxonomy
silva_ref <- "reference/silva_nr99_v138.2_train_set.fa.gz"

tax <- assignTaxonomy(
  seqtab.nochim,
  silva_ref,
  multithread = TRUE
)

# 13. Rename ASVs for readability
asv.seqs <- colnames(seqtab.nochim)
asv.headers <- paste0("ASV", seq_along(asv.seqs))
colnames(seqtab.nochim) <- asv.headers
rownames(tax) <- asv.headers

asv_lookup <- tibble(
  ASV = asv.headers,
  Sequence = asv.seqs
)

write_csv(asv_lookup, "results/asv_lookup.csv")
write_csv(as.data.frame(tax) |> rownames_to_column("ASV"), "results/taxonomy_table.csv")
write_csv(as.data.frame(seqtab.nochim) |> rownames_to_column("sample"), "results/asv_count_table.csv")

# 14. Import metadata
metadata <- read_csv("data/metadata.csv")
metadata_df <- metadata |>
  as.data.frame()

rownames(metadata_df) <- metadata_df[, 1]
metadata_df <- metadata_df[, -1, drop = FALSE]
metadata_df <- metadata_df[rownames(seqtab.nochim), , drop = FALSE]

# 15. Create phyloseq object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(as.matrix(tax)),
  sample_data(metadata_df)
)

saveRDS(ps, "results/phyloseq_object.rds")

# 16. Basic downstream analysis
plot_richness(ps, x = "group", measures = c("Observed", "Shannon")) +
  theme_bw()

ord <- ordinate(ps, method = "PCoA", distance = "bray")
plot_ordination(ps, ord, color = "group") +
  geom_point(size = 3) +
  theme_bw()

sessionInfo()
