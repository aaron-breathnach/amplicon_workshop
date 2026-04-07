#########################################
# Amplicon sequencing analysis workshop
#
# Pipeline: DADA2 + phyloseq
# Participant-facing teaching script
#########################################

# 0. Install packages (run once only)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Uncomment if needed
# BiocManager::install(c("dada2", "phyloseq", "Biostrings"))
# install.packages(c("ggplot2", "dplyr", "readr", "tibble", "stringr", "tidyr"))

# 1. Load libraries
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(tidyr)

# create a directory for output files
dir.create("results", showWarnings = FALSE)

# 2. Define input files
path <- "data/fastq"

## list forward reads
fnFs <- sort(
  list.files(
    path,
    pattern = "_L001_R1_001.fastq.gz$",
    full.names = TRUE
  )
)

## list reverse reads
fnRs <- sort(
  list.files(
    path,
    pattern = "_L001_R2_001.fastq.gz$",
    full.names = TRUE
  )
)

## stop early if no FASTQ files are found
if (length(fnFs) == 0 || length(fnRs) == 0) {
  stop("No FASTQ files found in data/fastq. Check that the files are in the correct folder.")
}

## check that we have the same number of forward and reverse read files
stopifnot(length(fnFs) == length(fnRs))

## read the first four lines of a fastq file
readLines(gzfile(fnFs[1]), n = 4)

## extract sample IDs from the forward read file names
sample.names <- basename(fnFs) |>
  str_replace("_L001_R1_001.fastq.gz$", "")

## extract sample IDs from the reverse read file names
sample.names.r <- basename(fnRs) |>
  str_replace("_L001_R2_001.fastq.gz$", "")

## check that forward and reverse files correspond to the same samples
stopifnot(all(sample.names == sample.names.r))

# 3. Inspect quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# 4. Filter and trim
## create a directory to write the filtered reads to
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

## names for the filtered forward reads
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

## names for the filtered reverse reads
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

## remove low-quality reads and trim low-quality tails
out <- filterAndTrim(
  fwd = fnFs,             # forward reads
  filt = filtFs,          # output filtered forward reads
  rev = fnRs,             # reverse reads
  filt.rev = filtRs,      # output filtered reverse reads
  truncLen = c(240, 160), # truncate reads at these positions (F, R)
  maxEE = c(2, 2),        # maximum expected errors
  truncQ = 2,             # truncate at the first very low-quality base
  rm.phix = TRUE,         # remove PhiX reads
  compress = TRUE,
  multithread = TRUE      # use multiple CPU cores
)

out_df <- as.data.frame(out) |>
  rownames_to_column("sample")

# 5. Learn error rates
## DADA2 learns the sequencing error profile from the data
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# 6. Dereplicate reads
## collapse identical reads together to speed up computation
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# 7. Infer ASVs
## infer true biological sequences from the error model
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# 8. Merge pairs
## merge forward and reverse reads for each sample
mergers <- mergePairs(
  dadaF = dadaFs,
  derepF = derepFs,
  dadaR = dadaRs,
  derepR = derepRs
)

# 9. Construct sequence table
## rows = samples, columns = ASVs
seqtab <- makeSequenceTable(mergers)

# 10. Remove chimeras
## remove chimeras (PCR artefacts that can inflate diversity)
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE
)

# 11. Track read retention
## useful for seeing where reads are lost during processing
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
## use a SILVA training set downloaded in advance
silva_ref <- "reference/silva_nr99_v138.2_toGenus_trainset.fa.gz"

stopifnot(file.exists(silva_ref))

tax <- assignTaxonomy(
  seqtab.nochim,
  silva_ref,
  multithread = TRUE
)

# 13. Rename ASVs for readability
## shorter ASV names are easier to read in plots and tables
asv.seqs <- colnames(seqtab.nochim)
asv.headers <- paste0("ASV", seq_along(asv.seqs))

colnames(seqtab.nochim) <- asv.headers
rownames(tax) <- asv.headers

asv_lookup <- tibble(
  ASV = asv.headers,
  Sequence = asv.seqs
)

write_csv(asv_lookup, "results/asv_lookup.csv")

write_csv(
  as.data.frame(tax) |> rownames_to_column("ASV"),
  "results/taxonomy_table.csv"
)

write_csv(
  as.data.frame(seqtab.nochim) |> rownames_to_column("sample"),
  "results/asv_count_table.csv"
)

# 14. Import metadata
metadata <- read_csv("data/metadata.csv")

metadata_df <- metadata |>
  as.data.frame()

rownames(metadata_df) <- metadata_df[, 1]
metadata_df <- metadata_df[, -1, drop = FALSE]

## check that all sequence-table samples are present in metadata
stopifnot(all(rownames(seqtab.nochim) %in% rownames(metadata_df)))

## reorder metadata so it matches the sequence table
metadata_df <- metadata_df[rownames(seqtab.nochim), , drop = FALSE]

# 15. Create phyloseq object
## combine counts, taxonomy, and metadata into one object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(as.matrix(tax)),
  sample_data(metadata_df)
)

saveRDS(ps, "results/phyloseq_object.rds")

# 16. Basic downstream analysis
meta <- as.data.frame(sample_data(ps))
meta$SampleID <- rownames(meta)

## alpha diversity analysis

## calculate alpha diversity metrics
richness <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

## add SampleID column
richness$SampleID <- rownames(richness)

## convert to long format for plotting
richness <- inner_join(richness, meta, by = "SampleID") |>
  pivot_longer(cols = -c(SampleID, time), names_to = "name", values_to = "value")

## visualise alpha diversity
ggplot(richness, aes(x = time, y = value)) +
  facet_wrap(~ name, scales = "free") +
  geom_boxplot(aes(fill = time)) +
  theme_bw()

## beta diversity analysis

## perform PCoA using Bray-Curtis distances
ord <- ordinate(ps, method = "PCoA", distance = "bray")

## extract the sample coordinates for plotting
ord_df <- as.data.frame(ord$vectors)[, 1:2]
ord_df$SampleID <- rownames(ord_df)
ord_df <- inner_join(ord_df, meta, by = "SampleID")

## percentage of variation captured by the first two PCoA axes
eig <- ord$values$Relative_eig[1:2] * 100

## each point is a sample; distances reflect community differences
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = time)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    x = paste0("PCoA 1 (", round(eig[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(eig[2], 1), "%)"),
    colour = "Timepoint"
  )