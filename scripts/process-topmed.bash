#!/bin/bash
#PBS -l select=1:ncpus=9:mem=128gb
#PBS -l walltime=72:00:00

module load htslib

bin_dir=/rds/general/user/gotto/home/apps/

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

mkdir -p ${out_dir}/imputed-nhlbi-topmed/merged/vcf

temp_dir=${out_dir}/temp_${PBS_JOBNAME}

mkdir -p $temp_dir/imputed-nhlbi-topmed


sex_chr=("X" "Y" "XY")

if [[ ${sex_chr[@]} =~ $chr ]]

then

    ## for X chromosome, which does not have exomeChip data
    ${bin_dir}/bcftools/bin/bcftools merge ${out_dir}/imputed-nhlbi-topmed/affymetrix/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-nhlbi-topmed/coreExome/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-nhlbi-topmed/gsa/vcf/chr${chr}.dose.vcf.gz --force-samples -O b -o ${temp_dir}/imputed-nhlbi-topmed/chr${chr}.dose.temp.vcf.gz

else
    
    ## for autosomes with exomeChip data, prioritising duplicates in
    ## order, i.e. if a sample is duplicated in gsa and exomeChip, the
    ## exomeChip sample name will get a 4:<sample-id>, which will be
    ## filtered out in the next step
    ${bin_dir}/bcftools/bin/bcftools merge ${out_dir}/imputed-nhlbi-topmed/affymetrix/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-nhlbi-topmed/coreExome/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-nhlbi-topmed/gsa/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-nhlbi-topmed/exomeChip/vcf/chr${chr}.dose.vcf.gz --force-samples -O b -o ${temp_dir}/imputed-nhlbi-topmed/chr${chr}.dose.temp.vcf.gz

fi

${bin_dir}/bcftools/bin/bcftools query -l ${temp_dir}/imputed-nhlbi-topmed/chr${chr}.dose.temp.vcf.gz | grep -Pv '^\d{5}$' > ${out_dir}/imputed-nhlbi-topmed/merged/vcf/chr${chr}-samples-removed.txt

${bin_dir}/bcftools/bin/bcftools view ${temp_dir}/imputed-nhlbi-topmed/chr${chr}.dose.temp.vcf.gz --samples-file ^${out_dir}/imputed-nhlbi-topmed/merged/vcf/chr${chr}-samples-removed.txt -O b -o ${out_dir}/imputed-nhlbi-topmed/merged/vcf/chr${chr}.dose.vcf.gz

tabix -p vcf ${out_dir}/imputed-nhlbi-topmed/merged/vcf/chr${chr}.dose.vcf.gz

rm -r $temp_dir
