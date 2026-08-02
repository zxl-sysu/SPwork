cat  ./Norsid.txt | cut -f 1,2,3,4 | awk -v OFS='\t' 'NR>0' - |  sort -k1n -k2n | awk 'NR>0' > chrposAle.txt


.

#CHR POS REF ALT
cut -f 1,2 chrposAle.txt |  bcftools query -T -  -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/RS\n" ./anno_reference/share/VCF/SNP/sites.All.norm.vcf.gz | awk -v OFS='\t' '{if (length($3) == 1  || NR==1) print $1,$2,$3,$4,"rs"$5}' > QC/anno_37_rsid.txt
#CHR POS REFANNO ALTANNO SNP REFRAW ALTRAW
awk -v OFS='\t' 'NR==FNR {a[$1","$2]=$3"\t"$4; next} {print $0, a[$1","$2]}' chrposAle.txt   QC/anno_37_rsid.txt > QC/anno_37id_Ale.txt



awk '( ($3==toupper($7) && $4==toupper($6) ) || ($4==toupper($7) && $3==toupper($6)  )) ' QC/anno_37id_Ale.txt  > QC/True_snp_gzAnno.txt

cut -f 1,2,5,6,7  QC/True_snp_gzAnno.txt | awk -v OFS='\t' '!seen[$3]++' - > QC/AnnoSNP

