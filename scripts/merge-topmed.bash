#!/bin/bash
#PBS -l select=1:ncpus=9:mem=128gb
#PBS -l walltime=72:00:00
#PBS -N merge_topmed

## concatenating chromosomes

cd $PBS_O_WORKDIR

bin_dir=/rds/general/user/gotto/home/apps/

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged" "merged_then_imputed")


for array in "${array_list[@]}"; do

    ## bed
    for bed_file in $(ls ${out_dir}/imputed-nhlbi-topmed/${array}/bed/${array}.chr[1-9]*.bed | sort -V); do  echo ${bed_file%.bed}; done > ${out_dir}/imputed-nhlbi-topmed/${array}/bed/merge-list-${array}.txt

    ${bin_dir}/plink2 --pmerge-list ${out_dir}/imputed-nhlbi-topmed/${array}/bed/merge-list-${array}.txt bfile --make-bed --out ${out_dir}/imputed-nhlbi-topmed/${array}/bed/${array}.autosomes


## bgen
## plink does not merge bgen files directly. Use bgenix
    # for bgen_file in $(ls ${out_dir}/imputed-nhlbi-topmed/${array}/bgen-1.2/${array}.chr[1-9]*.bgen | sort -V); do  echo ${bgen_file%.bgen}; done > ${out_dir}/imputed-nhlbi-topmed/${array}/bgen-1.2/merge-list-${array}.txt
		   
    # ${bin_dir}/plink2 --pmerge-list ${out_dir}/imputed-nhlbi-topmed/${array}/bgen-1.2/merge-list-${array}.txt --export bgen-1.2 --out ${out_dir}/imputed-nhlbi-topmed/${array}/bgen-1.2/${array}.autosomes

## pgen
    for pgen_file in $(ls ${out_dir}/imputed-nhlbi-topmed/${array}/pgen/${array}.chr[1-9]*.pgen | sort -V); do  echo ${pgen_file%.pgen}; done > ${out_dir}/imputed-nhlbi-topmed/${array}/pgen/merge-list-${array}.txt
		   
    ${bin_dir}/plink2 --pmerge-list ${out_dir}/imputed-nhlbi-topmed/${array}/pgen/merge-list-${array}.txt --out ${out_dir}/imputed-nhlbi-topmed/${array}/pgen/${array}.autosomes

done
