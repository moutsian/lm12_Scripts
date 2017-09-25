#!/usr/bin/bash

## note that in contrast with Javi's scripts, this works for both UC and IBD (just specify the trait parameter)

trait=$1
chr=$2;
main_snp_pos=$3

utrait=${trait^^} #upper case
ltrait=${trait,,} # lower case

#tmp=$(cat "/nfs/users/nfs_l/lm12/IBD_conditional/step2.ALL.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.auto.with_corr_pvals.txt" | awk '{print $4,$20,$21}')
#eval "$tmp"
loci=($(cat ~lm12/IBD_conditional/step2.ALL.NOVEL.round2.with_info_on_previous_hits_from_5PAPERS.auto.with_corr_pvals.txt | awk -v MAIN="${main_snp_pos}" '{if($4==MAIN){print $20,$21}}'))
echo "${loci[0]} and ${loci[1]}"
snp1=$(cat "/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.ALL.1e_05.filtered.txt" | awk -v POS1="${loci[0]}" '{if($3==POS1){print $2}}'|sed -n '1,1p')

snp2=$(cat "/lustre/scratch115/teams/anderson/ibd_conditional/final_meta_analysis_C_results/final_meta_analysis_C.ALL.1e_05.filtered.txt" | awk -v POS2="${loci[1]}" '{if($3==POS2){print $2}}'|sed -n '1,1p')


conditional_snps=" ${snp1}  ${snp2} "
if [ "${loci[1]}" = "NA" ];then
conditional_snps=" ${snp1} "
fi
if [ "${loci[0]}" = "NA" ];then
conditional_snps=" ${snp2}  "
fi

echo "Condition on : ${conditional_snps}"

startLoci=$(($main_snp_pos - 5))
endLoci=$(($main_snp_pos + 5))

opltrait="uc"
if [ "$ltrait" = "uc" ]; then
opltrait="cd"
fi


#IBDseq
study="IBDseq"
Javi_dir="/lustre/scratch113/teams/anderson/users/jga/001_Projects/Project2_FineMappingSeqData"
Yang_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL"
Data_dir="/lustre/scratch113/projects/crohns/RELEASE/v2/$study"
Sample_dir="/lustre/scratch113/projects/crohns/RELEASE/v1/IBDseq"
results_dir="/lustre/scratch115/teams/anderson/ibd_conditional"

#echo "/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "${conditional_snps}" -range $startLoci-$endLoci -exclude_samples $Javi_dir/${opltrait}SamplesAndFailSamples.23Jul2015.txt -exclude_snps $Javi_dir/SNPsPassFilters/${chr}.${utrait}.snps.NoPassFilters.OnlyIDs.txt"


if [ "$ltrait" = "ibd" ]; then
command=$(echo "/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "${conditional_snps}" -range $startLoci-$endLoci -exclude_samples /lustre/scratch113/projects/crohns/2013Aug07/assoc/beagle-v2/fail-qc-inds.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt")
eval "$command"
if [ "${run_unconditional}" = "YES" ]; then
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_samples /lustre/scratch113/projects/crohns/2013Aug07/assoc/beagle-v2/fail-qc-inds.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
fi
else
command=$(echo "/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "${conditional_snps}"  -range $startLoci-$endLoci -exclude_samples $Javi_dir/${opltrait}SamplesAndFailSamples.23Jul2015.txt -exclude_snps $Javi_dir/SNPsPassFilters/${chr}.${utrait}.snps.NoPassFilters.OnlyIDs.txt")
eval "$command"
if [ "${run_unconditional}" = "YES" ]; then
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_samples $Javi_dir/${opltrait}SamplesAndFailSamples.23Jul2015.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
fi
fi

#GWAS3
study="GWAS3"
Javi_dir="/lustre/scratch113/teams/anderson/users/jga/001_Projects/Project2_FineMappingSeqData"
Yang_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL"
Sample_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc"
Data_dir="/lustre/scratch113/projects/crohns/RELEASE/v2/${study}"
results_dir="/lustre/scratch115/teams/anderson/ibd_conditional"
if [ "$ltrait" = "ibd" ]; then
command=$(echo "/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "$conditional_snps" -range $startLoci-$endLoci -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt")
eval "$command"
if [ "${run_unconditional}" = "YES" ]; then
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
fi
else
command=$(echo "/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "$conditional_snps" -range $startLoci-$endLoci -exclude_samples $Javi_dir/${opltrait}SamplesAndFailSamples.23Jul2015.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt")
eval "$command"
if [ "${run_unconditional}" = "YES" ]; then
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_samples $Javi_dir/${opltrait}SamplesAndFailSamples.23Jul2015.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
fi
fi


if [ "$ltrait" = "uc" ] || [ "$ltrait" = "ibd" ] ; then
#GWAS2
study="GWAS2"
Sample_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc"
results_dir="/lustre/scratch115/teams/anderson/ibd_conditional"
Data_dir="/lustre/scratch113/projects/crohns/RELEASE/v2/${study}"
Javi_dir="/lustre/scratch113/teams/anderson/users/jga/001_Projects/Project2_FineMappingSeqData"
command=$(echo "/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous  -condition_on "${conditional_snps}" -range $startLoci-$endLoci -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt")
eval "$command"
if [ "${run_unconditional}" = "YES" ]; then
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
fi
fi

if [ "$ltrait" = "cd" ] || [ "$ltrait" = "ibd" ] ; then
#GWAS1
study="GWAS1"
Javi_dir="/lustre/scratch113/teams/anderson/users/jga/001_Projects/Project2_FineMappingSeqData"
Yang_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL"
Data_dir="/lustre/scratch113/projects/crohns/RELEASE/v2/$study"
Sample_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc"
results_dir="/lustre/scratch115/teams/anderson/ibd_conditional"
command=$(echo "/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "$conditional_snps" -range $startLoci-$endLoci -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt")
eval "$command"
if [ "${run_unconditional}" = "YES" ]; then
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}_with_pcs.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
fi
fi


