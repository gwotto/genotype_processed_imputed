#!/bin/bash
#PBS -l select=1:ncpus=9:mem=128gb
#PBS -l walltime=72:00:00


## Reformating vcf files to bed bgen-1.2 and ped format

cd $PBS_O_WORKDIR


out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

bin_dir=/rds/general/user/gotto/home/apps/

temp_dir=${out_dir}/temp_topmed_chr${chr}

mkdir -p $temp_dir

## current version of plink2 does not support multi-allelic sites in
## dosage format, needs --import-max-alleles 2

array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged" "merged_then_imputed")


for array in "${array_list[@]}"; do

    ## joining multiallelic sites that are split over multiple lines
    ${bin_dir}/bcftools/bin/bcftools norm ${out_dir}/imputed-nhlbi-topmed/${array}/vcf/chr${chr}.dose.vcf.gz -m +both -O b -o ${temp_dir}/${array}.chr${chr}.vcf.gz

    
    ## 1 bed    

    mkdir -p ${out_dir}/imputed-nhlbi-topmed/${array}/bed
    
    ${bin_dir}/plink2 --vcf ${temp_dir}/${array}.chr${chr}.vcf.gz --make-bed --double-id --update-sex ${out_dir}/barcode-sex.tsv --split-par hg19 --import-max-alleles 2 --vcf-half-call missing --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out ${out_dir}/imputed-nhlbi-topmed/${array}/bed/${array}.chr${chr}

    ## 2 bgen

    mkdir -p ${out_dir}/imputed-nhlbi-topmed/${array}/bgen-1.2
    
    ${bin_dir}/plink2 --vcf ${temp_dir}/${array}.chr${chr}.vcf.gz dosage=DS --export bgen-1.2 --double-id --update-sex ${out_dir}/barcode-sex.tsv --split-par hg19 --import-max-alleles 2 --vcf-half-call missing --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out ${out_dir}/imputed-nhlbi-topmed/${array}/bgen-1.2/${array}.chr${chr}
    
    ## 3 pgen

    mkdir -p ${out_dir}/imputed-nhlbi-topmed/${array}/pgen
    
    ${bin_dir}/plink2 --vcf ${temp_dir}/${array}.chr${chr}.vcf.gz dosage=DS --make-pgen --double-id --update-sex ${out_dir}/barcode-sex.tsv --split-par hg19 --import-max-alleles 2 --vcf-half-call missing --set-missing-var-ids @:#:\$r:\$a --new-id-max-allele-len 25 missing --out ${out_dir}/imputed-nhlbi-topmed/${array}/pgen/${array}.chr${chr}

done

rm -r $temp_dir
