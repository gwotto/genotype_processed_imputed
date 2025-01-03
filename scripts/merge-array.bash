#!/bin/bash
#PBS -l select=1:ncpus=9:mem=128gb
#PBS -l walltime=72:00:00
#PBS -N merge_array

## concatenating chromosomes

cd $PBS_O_WORKDIR

bin_dir=/rds/general/user/gotto/home/apps/

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged")

for array in "${array_list[@]}"; do

    ## bed
    for bed_file in $(ls ${out_dir}/fixed/${array}/bed/${array}.chr*.bed | sort -V); do  echo ${bed_file%.bed}; done > ${out_dir}/fixed/${array}/bed/merge-list-${array}.txt

    ${bin_dir}/plink2 --pmerge-list ${out_dir}/fixed/${array}/bed/merge-list-${array}.txt bfile --make-bed --out ${out_dir}/fixed/${array}/bed/${array}.combined


## bgen
## plink does not merge bgen files directly. Use bgenix
    for bgen_file in $(ls ${out_dir}/fixed/${array}/bgen-1.2/${array}.chr*.bgen | sort -V); do  echo ${bgen_file%.bgen}; done > ${out_dir}/fixed/${array}/bgen-1.2/merge-list-${array}.txt
		   
    ${bin_dir}/plink2 --pmerge-list ${out_dir}/fixed/${array}/bgen-1.2/merge-list-${array}.txt --export bgen-1.2 --out ${out_dir}/fixed/${array}/bgen-1.2/${array}.combined

## pgen
    for pgen_file in $(ls ${out_dir}/fixed/${array}/pgen/${array}.chr*.pgen | sort -V); do  echo ${pgen_file%.pgen}; done > ${out_dir}/fixed/${array}/pgen/merge-list-${array}.txt
		   
    ${bin_dir}/plink2 --pmerge-list ${out_dir}/fixed/${array}/pgen/merge-list-${array}.txt --out ${out_dir}/fixed/${array}/pgen/${array}.combined

done
