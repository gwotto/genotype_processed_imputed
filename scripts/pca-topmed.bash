#!/bin/bash
#PBS -l select=1:ncpus=9:mem=128gb
#PBS -l walltime=72:00:00
#PBS -N pca_topmed

## concatenating chromosomes

cd $PBS_O_WORKDIR

bin_dir=/rds/general/user/gotto/home/apps/

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

array_list=("affymetrix" "coreExome" "exomeChip" "gsa" "merged" "merged_then_imputed")


for array in "${array_list[@]}"; do

    ${bin_dir}/plink2 --bfile ${out_dir}/imputed-nhlbi-topmed/${array}/bed/${array}.autosomes --maf 0.01 --indep-pairwise 100kb 0.8 --make-bed --rm-dup force-first list --out ${out_dir}/imputed-nhlbi-topmed/${array}/bed/${array}.autosomes.pruned 


    ${bin_dir}/plink2 --threads 8 --bfile ${out_dir}/imputed-nhlbi-topmed/${array}/bed/${array}.autosomes.pruned --pca 20 --out ${out_dir}/imputed-nhlbi-topmed/${array}/bed/${array}.autosomes

done
