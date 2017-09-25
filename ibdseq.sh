#!/usr/bin/bash

## note that in contrast with Javi's scripts, this works for both UC and IBD (just specify the trait parameter)

trait=$1
chr=$2;
main_snp_pos=$3
#6:69009744_T_C
conditional_snp1=6:32806553_G_A
conditional_snp2=6:90931858_C_T
#conditional_snp1=$4
#conditional_snp2=$5

startLoci=$(())
startLoci=$(($main_snp_pos - 50))
endLoci=$(($main_snp_pos + 50))

utrait=${trait^^} #upper case
ltrait=${trait,,} # lower case

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

if [ "$ltrait" = "ibd" ]; then
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "$conditional_snp1" "$conditional_snp2" -range $startLoci-$endLoci -exclude_samples /lustre/scratch113/projects/crohns/2013Aug07/assoc/beagle-v2/fail-qc-inds.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$ltrait.snps.NoPassFilters.OnlyIDs.txt
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_samples /lustre/scratch113/projects/crohns/2013Aug07/assoc/beagle-v2/fail-qc-inds.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$ltrait.snps.NoPassFilters.OnlyIDs.txt
else
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on "$conditional_snp1" "$conditional_snp2" -range $startLoci-$endLoci -exclude_samples $Javi_dir/${opltrait}SamplesAndFailSamples.23Jul2015.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/${study}.sample -o $results_dir/$study.$utrait.$chr.$startLoci.$endLoci.unconditional.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -range $startLoci-$endLoci -exclude_samples $Javi_dir/${opltrait}SamplesAndFailSamples.23Jul2015.txt -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$utrait.snps.NoPassFilters.OnlyIDs.txt
fi


