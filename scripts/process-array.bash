
#!/bin/bash
#PBS -l select=1:ncpus=8:mem=100gb
#PBS -l walltime=72:00:00

cd $PBS_O_WORKDIR

# chr=21

bin_dir=/rds/general/user/gotto/home/apps/

export BCFTOOLS_PLUGINS=${bin_dir}/bcftools/libexec/bcftools/

out_dir=/rds/general/project/airwave/live/genotype_processed_imputed/

mkdir -p ${out_dir}/fixed/merged/vcf

temp_dir=${out_dir}/temp_chr${chr}

mkdir -p $temp_dir

module load htslib


# == convert bed files to vcf, replace affy IDs with barcode ==


# == affymetrix == 

## vcf file of chr ${chr}, affymetrix, raw genotype, --snps-only just-acgt, deletes any illegal bases
${bin_dir}/plink2 --bfile /rds/general/project/airwave/live/affymetrix_genotype/all.clean-base --chr ${chr} --allow-no-sex --recode vcf bgz --snps-only just-acgt --split-par hg19 --out ${temp_dir}/affymetrix.chr${chr}.id

## convert affymetrix sample IDs to barcode IDs, remove sample that can not be converted
${bin_dir}/bcftools/bin/bcftools reheader --samples data/affy_sample_mapping.txt ${temp_dir}/affymetrix.chr${chr}.id.vcf.gz | ${bin_dir}/bcftools/bin/bcftools view -s ^SP0001IC1028560Z_a550484-4279740-040717-255_A09 -O b -o ${temp_dir}/affymetrix.chr${chr}.barcode.vcf.gz

tabix -p vcf ${temp_dir}/affymetrix.chr${chr}.barcode.vcf.gz


# == coreExome ===

## vcf file of chr ${chr}, coreExome, raw genotype, --recode vcf id-paste=iid bgz --snps-only just-acgt
${bin_dir}/plink2 --bfile /rds/general/project/airwave/live/coreExome_genotype/all.clean-base --chr ${chr} --allow-no-sex --recode vcf id-paste=iid bgz --snps-only just-acgt --split-par hg19 --out ${temp_dir}/coreExome.chr${chr}.id

## remove duplicated samples
${bin_dir}/bcftools/bin/bcftools view -s ^11421_rep,46750_rep  ${temp_dir}/coreExome.chr${chr}.id.vcf.gz -O b -o ${temp_dir}/coreExome.chr${chr}.barcode.vcf.gz

tabix -p vcf ${temp_dir}/coreExome.chr${chr}.barcode.vcf.gz

# == exomeChip ==

## vcf file of chr ${chr}, exomeChip, raw genotype, --recode vcf id-paste=iid bgz --snps-only just-acgt

## samples "47781" and  "61960" are duplicated between exomeChip and gsa, see below in qc

${bin_dir}/plink2 --bfile /rds/general/project/airwave/live/exomeChip_genotype/zcalled_QCed_Final --chr ${chr} --allow-no-sex --recode vcf id-paste=iid bgz --snps-only just-acgt --split-par hg19 --out ${temp_dir}/exomeChip.chr${chr}.barcode

tabix -p vcf ${temp_dir}/exomeChip.chr${chr}.barcode.vcf.gz


# == gsa ==

## vcf file of chr ${chr}, gsa, raw genotype --recode vcf id-paste=iid bgz --snps-only just-acgt
${bin_dir}/plink2 --bfile /rds/general/project/airwave/live/QC/gsa_zcall_b37_maf_001_2021/all.clean-base --chr ${chr} --allow-no-sex --recode vcf id-paste=iid bgz --snps-only just-acgt --split-par hg19 --out ${temp_dir}/gsa.chr${chr}.id

## convert IDs to barcode IDs (reheader output is uncompressed)
${bin_dir}/bcftools/bin/bcftools reheader --samples data/gsa-reheader.txt ${temp_dir}/gsa.chr${chr}.id.vcf.gz | ${bin_dir}/bcftools/bin/bcftools view -O b -o ${temp_dir}/gsa.chr${chr}.barcode.vcf.gz


tabix -p vcf ${temp_dir}/gsa.chr${chr}.barcode.vcf.gz



array_list=("affymetrix" "coreExome" "exomeChip" "gsa")


for array in "${array_list[@]}"; do

    mkdir -p ${out_dir}/fixed/${array}/vcf

    # == annotate with dbSNP, then replace the remaining non-rs IDs with position IDs ==

    ${bin_dir}/bcftools/bin/bcftools annotate -a /rds/general/project/uk-biobank-2020/live/resources/reference/GRCh37/dbsnp/00-All.vcf.gz -c ID ${temp_dir}/${array}.chr${chr}.barcode.vcf.gz | ${bin_dir}/bcftools/bin/bcftools annotate --include 'ID!~"^rs"' --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' --keep-sites -Ob --output ${temp_dir}/${array}.chr${chr}.annot.vcf.gz

    tabix -p vcf ${temp_dir}/${array}.chr${chr}.annot.vcf.gz


    ## fixing reference allele missmatches, flip mode

    ${bin_dir}/bcftools/bin/bcftools +fixref ${temp_dir}/${array}.chr${chr}.annot.vcf.gz -Ob -o ${temp_dir}/${array}.chr${chr}.flip.vcf.gz -- -d -f /rds/general/project/uk-biobank-2020/live/resources/reference/GRCh37/ensembl/Homo_sapiens.GRCh37.dna.toplevel.fa -i /rds/general/project/uk-biobank-2020/live/resources/reference/GRCh37/dbsnp/00-All.vcf.gz -m flip

    ## joining multiallelic sites that are split over multiple lines
    ${bin_dir}/bcftools/bin/bcftools norm ${temp_dir}/${array}.chr${chr}.flip.vcf.gz -m +both -O b -o ${out_dir}/fixed/${array}/vcf/${array}.chr${chr}.vcf.gz
    
    tabix -p vcf ${out_dir}/fixed/${array}/vcf/${array}.chr${chr}.vcf.gz

## 
    
done


# == merge fixed files

${bin_dir}/bcftools/bin/bcftools merge ${out_dir}/fixed/affymetrix/vcf/affymetrix.chr${chr}.vcf.gz ${out_dir}/fixed/coreExome/vcf/coreExome.chr${chr}.vcf.gz ${out_dir}/fixed/exomeChip/vcf/exomeChip.chr${chr}.vcf.gz ${out_dir}/fixed/gsa/vcf/gsa.chr${chr}.vcf.gz --force-samples -O z -o ${temp_dir}/merged.chr${chr}.tmp.vcf.gz


## detect sample IDs not in barcode format and remove (e.g. duplicates)
${bin_dir}/bcftools/bin/bcftools query -l ${temp_dir}/merged.chr${chr}.tmp.vcf.gz | grep -Pv '^\d{5}$' > ${out_dir}/fixed/merged/vcf/chr${chr}-samples-removed.txt

${bin_dir}/bcftools/bin/bcftools view ${temp_dir}/merged.chr${chr}.tmp.vcf.gz --samples-file ^${out_dir}/fixed/merged/vcf/chr${chr}-samples-removed.txt -O b -o ${out_dir}/fixed/merged/vcf/merged.chr${chr}.vcf.gz

tabix -p vcf ${out_dir}/fixed/merged/vcf/merged.chr${chr}.vcf.gz

rm -r $temp_dir
