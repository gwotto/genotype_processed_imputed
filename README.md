# Processing and imputation of Airwave Genotype Data

## Genotype processing

- Processing of array genotype data distributed over several
  directories in airwave project directory, using script
  `process-aray.bash`. Output in vcf format.

- Running as batch jobs with one job per chromosome (including X
  chromosome)

```
  for i in {1..22} X XY Y; do qsub -v chr=${i} -N process_array_chr${i} analysis/process-array.bash; done
``

- restricts variants to SNP, restricts to ACGT SNP, split
  pseudo-autosomal region on X, converts sample names to Airwave
  barcodes, removes sample duplicates
- annotates variants with dbSNP, replaces remaining non-rs IDs with
  position IDs
- fixes reference allele missmatches, flip mode (bcftools)
- joins multiallelic sites that are split over multiple lines
- merges data from arrays, removes duplicate samples.
