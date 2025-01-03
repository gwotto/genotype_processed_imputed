#!/bin/bash
#PBS -l select=1:ncpus=9:mem=124gb
#PBS -l walltime=72:00:00


## Reformating vcf files to bed and bgen-1.2 format

cd $PBS_O_WORKDIR

# chr=22

bin_dir=/rds/general/user/gotto/home/apps/

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

temp_dir=${out_dir}/temp_${PBS_JOBNAME}

mkdir -p $temp_dir

array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged")


for array in "${array_list[@]}"; do

    mkdir -p ${out_dir}/fixed/${array}/bed

    mkdir -p ${out_dir}/fixed/${array}/bgen-1.2

    mkdir -p ${out_dir}/fixed/${array}/pgen

## bed and bgen-1.2 can not handle multiallelic sites
    
# ## 1.1 bed 

    ${bin_dir}/plink2 --vcf /rds/general/project/airwave/live/genotype_processed_imputed/fixed/${array}/vcf/${array}.chr${chr}.vcf.gz --make-bed --double-id --update-sex ${out_dir}/barcode-sex.tsv --split-par hg19 --max-alleles 2 --vcf-half-call missing --out /rds/general/project/airwave/live/genotype_processed_imputed/fixed/${array}/bed/${array}.chr${chr}
    

# ## 1.2 bgen

    ${bin_dir}/plink2 --vcf /rds/general/project/airwave/live/genotype_processed_imputed/fixed/${array}/vcf/${array}.chr${chr}.vcf.gz --export bgen-1.2 --double-id --update-sex ${out_dir}/barcode-sex.tsv --split-par hg19 -max-alleles 2 --vcf-half-call missing --out /rds/general/project/airwave/live/genotype_processed_imputed/fixed/${array}/bgen-1.2/${array}.chr${chr}


# ## 1.3 pgen

    ${bin_dir}/plink2 --vcf /rds/general/project/airwave/live/genotype_processed_imputed/fixed/${array}/vcf/${array}.chr${chr}.vcf.gz --make-pgen --double-id --update-sex ${out_dir}/barcode-sex.tsv --split-par hg19 --vcf-half-call missing --out /rds/general/project/airwave/live/genotype_processed_imputed/fixed/${array}/pgen/${array}.chr${chr}

done

rm -r $temp_dir


## TODO the part below should go into a script reformat-beagle.bash

# ## 2 imputed


# ## current version of plink2 does not support multi-allelic sites in bed and bgen, needs --max-alleles 2

## regardless of that, however, using with dosage gives error when
## having multi-allelic sites in input vcf, have to be removed from
## vcf beforehand

# array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged")


# for array in "${array_list[@]}"; do


# mkdir -p ${out_dir}/imputed-beagle-1kg_v3/${array}/bed

# mkdir -p ${out_dir}/imputed-beagle-1kg_v3/${array}/bgen-1.2

# mkdir -p ${out_dir}/imputed-beagle-1kg_v3/${array}/pgen
    
# ## remove multi-allelic snps
# /rds/general/user/gotto/home/apps/bcftools/bin/bcftools view --max-alleles 2 ${out_dir}/imputed-beagle-1kg_v3/${array}/vcf/${array}.chr${chr}.1kg.vcf.gz -O b -o ${temp_dir}/${array}.chr${chr}.1kg.vcf.gz

# ## bed
# # /rds/general/user/gotto/home/apps/plink2 --vcf ${temp_dir}/${array}.chr${chr}.1kg.vcf.gz --make-bed --double-id --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out /rds/general/project/airwave/live/genotype_processed_imputed/imputed-beagle-1kg_v3/${array}/bed/${array}.chr${chr}.1kg

# ## bgen
# /rds/general/user/gotto/home/apps/plink2 --vcf ${temp_dir}/${array}.chr${chr}.1kg.vcf.gz dosage=DS --export bgen-1.2 --double-id --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out /rds/general/project/airwave/live/genotype_processed_imputed/imputed-beagle-1kg_v3/${array}/bgen-1.2/${array}.chr${chr}.1kg

# ## ped
# /rds/general/user/gotto/home/apps/plink2 --vcf ${temp_dir}/${array}.chr${chr}.1kg.vcf.gz dosage=DS --make-pgen --double-id --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out /rds/general/project/airwave/live/genotype_processed_imputed/imputed-beagle-1kg_v3/${array}/pgen/${array}.chr${chr}.1kg

# done
