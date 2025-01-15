#!/bin/bash
#PBS -l select=1:ncpus=8:mem=100gb
#PBS -l walltime=72:00:00


## Reformation vcf files to bed and bgen-1.2 format

cd $PBS_O_WORKDIR

# chr=22

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

mkdir -p ${out_dir}/imputed-umich-1kg/affymetrix/bed
mkdir -p ${out_dir}/imputed-umich-1kg/coreExome/bed
mkdir -p ${out_dir}/imputed-umich-1kg/exomeChip/bed
mkdir -p ${out_dir}/imputed-umich-1kg/gsa/bed
mkdir -p ${out_dir}/imputed-umich-1kg/merged/bed

mkdir -p ${out_dir}/imputed-umich-1kg/affymetrix/bgen-1.2
mkdir -p ${out_dir}/imputed-umich-1kg/coreExome/bgen-1.2
mkdir -p ${out_dir}/imputed-umich-1kg/exomeChip/bgen-1.2
mkdir -p ${out_dir}/imputed-umich-1kg/gsa/bgen-1.2
mkdir -p ${out_dir}/imputed-umich-1kg/merged/bgen-1.2

mkdir -p ${out_dir}/imputed-umich-1kg/affymetrix/pgen
mkdir -p ${out_dir}/imputed-umich-1kg/coreExome/pgen
mkdir -p ${out_dir}/imputed-umich-1kg/exomeChip/pgen
mkdir -p ${out_dir}/imputed-umich-1kg/gsa/pgen
mkdir -p ${out_dir}/imputed-umich-1kg/merged/pgen


## current version of plink2 does not support multi-allelic sites, needs --max-alleles 2


array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged")


for array in "${array_list[@]}"; do

## is --max-alleles 2 necessary?
    
    ## 1 bed    
    
    # /rds/general/user/gotto/home/apps/plink2 --vcf ${out_dir}/imputed-umich-1kg/${array}/vcf/chr${chr}.dose.vcf.gz --make-bed --double-id --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out ${out_dir}/imputed-umich-1kg/${array}/bed/${array}.chr${chr}

    ## 2 bgen
    /rds/general/user/gotto/home/apps/plink2 --vcf ${out_dir}/imputed-umich-1kg/${array}/vcf/chr${chr}.dose.vcf.gz dosage=DS --export bgen-1.2 --double-id --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out ${out_dir}/imputed-umich-1kg/${array}/bgen-1.2/${array}.chr${chr}
    
    ## 3 pgen
    /rds/general/user/gotto/home/apps/plink2 --vcf ${out_dir}/imputed-umich-1kg/${array}/vcf/chr${chr}.dose.vcf.gz dosage=DS --make-pgen --double-id --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out ${out_dir}/imputed-umich-1kg/${array}/pgen/${array}.chr${chr}

done
