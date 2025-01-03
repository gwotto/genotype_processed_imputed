# Processing and imputation of Airwave Genotype Data

This is the workflow that was carried out to process genotype data of
the four array platforms used for Airwave participants and to imute
variants using different imputation panels.

## Genotype processing

- Processing of array genotype data distributed over several
  directories in airwave project directory, using script
  `process-aray.bash`. Output in vcf format.

- Running as batch jobs with one job per chromosome (including X
  chromosome)

```
  for i in {1..22} X XY Y; do qsub -v chr=${i} -N process_array_chr${i} analysis/process-array.bash; done
```

- restricts variants to SNP, restricts to ACGT SNP, split
  pseudo-autosomal region on X, converts sample names to Airwave
  barcodes, removes sample duplicates
- annotates variants with dbSNP, replaces remaining non-rs IDs with
  position IDs
- fixes reference allele missmatches, flip mode (bcftools)
- joins multiallelic sites that are split over multiple lines
- merges data from arrays, removes duplicate samples.

### Reformatting genotypes

- Converting vcf files to bed, bgen-1.2 and ped format

- Running as batch jobs with one job per chromosome (including X
  chromosome)

```
  for i in {1..22} X; do qsub -v chr=${i} -N process_array_chr${i} analysis/reformat-array.bash; done
```

- restricts variants to bi-allelic SNPs, includes sex variable, splits
  pseudo-autosomal region on X, treats '0/.' GT calls as missing

### Merging chromosomes

- Concatenating chromosomes in the correct order of each array,
  including the merged arrays.
  
- Runs as a batch job

```
qsub analysis/merge-array.bash
```

- Concatenates data in bed and ped format. Does not concatenate
  bgen-1.2, which needs bgenix, and vcf, because they are too large.
  
