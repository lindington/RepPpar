### turning snpselection outfile into sites file 
for l in linked random; do 
    for pop in parent all; do
        sed 's/ /\n/g' ${l}snps_${pop}.txt > trylink
        awk 'BEGIN{FS=OFS="\t"}{print "chr1", $0,$1}' trylink > ${l}snps_${pop}.sites
        rm trylink
        rm ${l}snps_${pop}.txt
    done
done

module load plink
module load plink2
module load bcftools
module load htslib

### Plink tped and tfam files to vcf conversion
cd 01.output
for l in linked random; do 
    for pop in parent all; do
        plink --noweb --tfile vcf_dadi_${l}_${pop} --make-bed --out vcf_dadi_${l}_${pop}
        plink2 --bfile vcf_dadi_${l}_${pop} --recode vcf --out dadi_${l}_${pop}
        bgzip dadi_${l}_${pop}.vcf
        bcftools index dadi_${l}_${pop}.vcf.gz
    done
done