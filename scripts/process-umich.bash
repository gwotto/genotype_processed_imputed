#!/bin/bash
#PBS -l select=1:ncpus=8:mem=100gb
#PBS -l walltime=72:00:00

module load htslib


bin_dir=/rds/general/user/gotto/home/apps/

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

mkdir -p ${out_dir}/imputed-umich-1kg/merged/vcf

temp_dir=${out_dir}/temp_${PBS_JOBNAME}

mkdir -p $temp_dir/imputed-umich-1kg



${bin_dir}/bcftools/bin/bcftools merge ${out_dir}/imputed-umich-1kg/affymetrix/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-umich-1kg/coreExome/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-umich-1kg/exomeChip/vcf/chr${chr}.dose.vcf.gz ${out_dir}/imputed-umich-1kg/gsa/vcf/chr${chr}.dose.vcf.gz --force-samples -O b -o ${temp_dir}/imputed-umich-1kg/chr${chr}.dose.temp.vcf.gz

${bin_dir}/bcftools/bin/bcftools query -l ${temp_dir}/imputed-umich-1kg/chr${chr}.dose.temp.vcf.gz | grep -Pv '^\d{5}$' > ${out_dir}/imputed-umich-1kg/merged/vcf/chr${chr}-samples-removed.txt

${bin_dir}/bcftools/bin/bcftools view ${temp_dir}/imputed-umich-1kg/chr${chr}.dose.temp.vcf.gz --samples-file ^${out_dir}/imputed-umich-1kg/merged/vcf/chr${chr}-samples-removed.txt -O b -o ${out_dir}/imputed-umich-1kg/merged/vcf/chr${chr}.dose.vcf.gz

tabix -p vcf ${out_dir}/imputed-umich-1kg/merged/vcf/chr${chr}.dose.vcf.gz


rm -r $temp_dir
