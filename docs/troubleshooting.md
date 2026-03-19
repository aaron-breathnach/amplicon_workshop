# Troubleshooting

## DADA2 will not install

Common causes:

- outdated R version
- missing Bioconductor dependencies
- internet/firewall issues

Try updating R, reinstalling `BiocManager`, and running the install command again.

## FASTQ files are not being detected

Check that:

- files are in `data/fastq/`
- file names end in `_R1.fastq.gz` and `_R2.fastq.gz`

## Forward and reverse files do not match

Make sure sample names are consistent between paired FASTQ files.

## Too many reads are lost during filtering

Possible reasons:

- `truncLen` is too aggressive
- read quality is poor
- `maxEE` is too strict

Inspect the quality profiles and adjust parameters carefully.

## Paired reads fail to merge

Possible reasons:

- insufficient overlap after truncation
- reverse reads are too poor quality
- amplicon is too long

## Taxonomy assignment fails

Check that the reference file path is correct and that the SILVA training set has been downloaded.
