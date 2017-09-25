#!/usr/bin/bash

## note that in contrast with Javi's scripts, this works for both UC and IBD (just specify the trait parameter)

trait=$1
chr=$2;
main_snp_pos=$3

utrait=${trait^^} #upper case
ltrait=${trait,,} # lower case

tmp=$(cat "/nfs/users/nfs_l/lm12/IBD_conditional/step2.ALL.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.auto.with_corr_pvals.txt" | grep "$main_snp_pos" |awk '{print $4,$20,$21}')
echo "$tmp"
loci=($(cat ~lm12/IBD_conditional/step2.ALL.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.auto.with_corr_pvals.txt | awk -v MAIN="${main_snp_pos}" '{if($4==MAIN){print $20,$21}}'))
echo "${loci[0]} and ${loci[1]}"
snp1=$(cat "/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.${utrait}.1e_05.filtered.txt" | awk -v POS1="${loci[0]}" '{if($3==POS1){print $2}}')

snp2=$(cat "/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.${utrait}.1e_05.filtered.txt" | awk -v POS2="${loci[1]}" '{if($3==POS2){print $2}}')


conditional_snps=" ${snp1}  ${snp2} "
if [ "${loci[1]}" = "NA" ];then
conditional_snps=" ${snp1} "
fi
if [ "${loci[0]}" = "NA" ];then
conditional_snps=" ${snp2}  "
fi

echo "Condition on : ${conditional_snps}"

