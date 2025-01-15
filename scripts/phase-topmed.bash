#!/bin/bash
#PBS -l select=1:ncpus=9:mem=128gb
#PBS -l walltime=72:00:00


## Reformating vcf files to bed bgen-1.2 and ped format

cd $PBS_O_WORKDIR


out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

bin_dir=/rds/general/user/gotto/home/apps/

temp_dir=${out_dir}/temp_shapeit_chr${chr}

mkdir -p $temp_dir


array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged")


for array in "${array_list[@]}"; do

    ## filter and convert vcf files to vcf 4.2 format
    ${bin_dir}/plink2 --export vcf-4.2 bgz --mind 0.05 --geno 0.05 --vcf-half-call missing --split-par hg19 --out ${temp_dir}/${array}.chr${chr}.ed --update-sex ${out_dir}/barcode-sex.tsv  --vcf ${out_dir}/fixed/${array}/vcf/${array}.chr${chr}.vcf.gz 

    ## phase genotypes using shapeit 
    if [ "$chr" == "X" ]; then
	${bin_dir}/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -V ${temp_dir}/${array}.chr${chr}.ed.vcf.gz -O ${temp_dir}/${array}.chr${chr}.phased --chrX --thread 8
    else
	
	echo "chr is not 'X'. Program will not run."
	## put code here for autosomes (not done)
	break
    fi

    ## convert output to vcf
    ${bin_dir}/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert --input-haps ${temp_dir}/${array}.chr${chr}.phased --output-vcf ${temp_dir}/${array}.chr{chr}.phased.vcf

    ${bin_dir}/bcftools/bin/bcftools view ${temp_dir}/${array}.chr{chr}.phased.vcf -O b -o ${out_dir}/fixed//${array}/vcf/${array}.chr${chr}.phased.vcf.gz

    tabix -p vcf ${out_dir}/fixed/${array}/vcf/${array}.chr${chr}.phased.vcf.gz

done

rm -r $temp_dir
