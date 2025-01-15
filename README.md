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

### Concatenating chromosomes

- Concatenating chromosomes in the correct order of each array,
  including the merged arrays.
  
- Runs as a batch job

```
qsub analysis/merge-array.bash
```

- Concatenates data in bed and ped format. Does not concatenate
  bgen-1.2, which needs bgenix, and vcf, because they are too large.
  

## 1000 Genomes imputation

- Imputation on UMICH imputation server https://imputationserver.sph.umich.edu

  - panel 1000G Phase 3 v5 (GRCh37/hg19)
  - rsq filter off
  - beagle phasing
  - population EUR

- Imputed output in vcf format, chr*.dose.vcf.gz


### Merging arrays

- Merging imputed arrays, vcf output in ```merged``` directory

- In the merging process, duplicate sample names get a reformated
  sample ID: ```<array>:<barcode>```, e.g. 4:47781 in this case
  participant 47781 on array 4 (gsa) is a duplicate from another
  array. Such duplicates and other participants with non-standard
  5-digit barcode IDs get removed (282 samples)

```
for i in {1..22}; do qsub -v chr=${i} -N merge_umich_chr${i} analysis/process-umich.bash; done
```

### Reformatting genotypes

- Converting vcf files to bed, bgen-1.2 and ped format

- Assigning position-allele-based IDs to variants with missing variant IDs 

- Variants (indels) with position-allele based ID and allele length > 25bp are assigned the unnamed-variant ID

- Running as batch jobs with one job per chromosome 

```
  for i in {1..22}; do qsub -v chr=${i} -N process_array_chr${i} analysis/reformat-umich.bash; done
```


### Concatenating chromosomes

- Concatenating chromosomes in the correct order of each array,
  including the merged arrays.

- Concatenating data in bed format. Not concatenate bgen-1.2, which
  needs bgenix, vcf, because they are too large, and pgen, because
  concatenation gave a segmentation fault
  
- Runs as a batch job

```
qsub analysis/merge-umich.bash
```

### PCA

- PCA of concatenated chromosomes in bed format

- Pruning variants, maf threshold 0.01 and --indep-pairwise with a 100kb window and an r2 threshold 0.8

- 20 PCAs

```
qsub analysis/pca-umich.bash
```


## Topmed imputation

- Imputation on NHLBI topmed imputation server https://imputation.biodatacatalyst.nhlbi.nih.gov

  - panel version r3 
  - missmatch control vs topmed population
  - eagle phasing
  - rsq 0.3 filter
  - hg38 reference

- Imputed output in vcf format, chr*.dose.vcf.gz

### X-chromosome imputation

- Imputation server rejects X-chromosome genotypes: ChrX nonPAR region
  includes ambiguous samples (haploid and diploid positions).

- Using shapeit to phase and reformat X chromosome vcf

```
qsub analysis/phase-topmed
```



### Merging arrays

- Merging imputed arrays, vcf output in ```merged``` directory

- In the merging process, duplicate sample names get a reformated
  sample ID: ```<array>:<barcode>```, e.g. 4:47781 in this case
  participant 47781 on array 4 (exomeChip) is a duplicate from another
  array. Such duplicates and other participants with non-standard
  5-digit barcode IDs get removed (282 samples)

- exomeChip does not contain X-chromosome data

```
for i in {1..22} X; do qsub -v chr=${i} -N merge_umich_chr${i} analysis/process-topmed.bash; done
```

### Reformatting genotypes

- Converting vcf files to bed, bgen-1.2 and ped format

- Joining multiallelic sites that are split over multiple lines
- Current version of plink2 does not support multi-allelic sites in
  dosage format, restricting to bi-allelic sites
- Updating sex, splitting pseudoautosomal region on X-chromosome
  (```--split-par hg19```), settin '0/.' and similar GT values as
  missing

- Assigning position-allele-based IDs to variants with missing variant IDs 

- Variants (indels) with position-allele based ID and allele length > 25bp are assigned the unnamed-variant ID

- Running as batch jobs with one job per chromosome 

```
  for i in {1..22}; do qsub -v chr=${i} -N process_array_chr${i} analysis/reformat-topmed.bash; done
```


### Concatenating chromosomes

- Concatenating autosomes in the correct order of each array, but excluding X-chromosome.

- Concatenating data in bed format. Not concatenate bgen-1.2, which
  needs bgenix, vcf, because they are too large.
  
- Runs as a batch job

```
qsub analysis/merge-topmed.bash
```

### PCA

- PCA of concatenated autosomes in bed format, excluding X-chromosome

- Pruning variants, maf threshold 0.01 and --indep-pairwise with a 100kb window and an r2 threshold 0.8

- 20 PCAs

```
qsub analysis/pca-topmed.bash
```
