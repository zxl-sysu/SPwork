
#CHR POS37 CHR POS38
cat TSS_chrPOS38.txt | awk -v OFS='\t' '{ print "chr"$1, $2, $2+1, "chr"$1":"$2}'| ./liftOver  stdin ./hg38ToHg19.over.chain /dev/stdout unmapped.txt | awk '{gsub(/^chr/, "", $1);gsub(/^chr/, "", $4);gsub(/:/, "\t", $4); print $1 "\t" $2 "\t" $4 }'  > QC/lifted_TSSto37.txt


