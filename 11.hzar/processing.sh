module load plink
module load plink2
module load bcftools
module load htslib

### Plink tped and tfam files to vcf conversion
plink --noweb --tfile 01.output/hzar_plink_allinds --make-bed --out 01.output/vcf_hzar_allinds
plink2 --bfile 01.output/vcf_hzar_allinds --recode vcf --out 01.output/hzar_allinds
conda activate bcftools
bgzip 01.output/hzar_allinds.vcf
bcftools index 01.output/hzar_allinds.vcf.gz

### rename inds using a list of all new sample names 
bcftools reheader -s allinds.txt 01.output/hzar_allinds.vcf.gz -o 01.output/hzar_allinds_renamed.vcf.gz

#####make structure file: need plink1.9 for this!!
## conda activate plink1.9
gzip -d 01.output/hzar_allinds.vcf.gz
plink --vcf 01.output/hzar_allinds_renamed.vcf --recode structure --out 01.output/hzar_allinds_renamed

### vcftools allele frequencies
for pop in GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC; do
    vcftools --vcf 01.output/hzar_allinds_renamed.vcf --keep pops/pop_${pop}.txt --positions snps.txt --freq --max-missing 0 --out 01.output/vcf_af_${pop}
done

### make full vcf for maria
vcftools --vcf 01.output/hzar_allinds_renamed.vcf --keep pops/por_transect.txt --positions snps.txt --recode --max-missing 0 --out 01.output/vcf_por_transect
