
#!/bin/bash
assoc=$1
#chr=$2
#start=$3
#end=$3


#unconditional
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5 \
-data /lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.gen /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample \
-o "/lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.uncond" \
-frequentist 1 \
-method score \
-pheno bin1 \
-exclude_samples "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/gwas3-$assoc-assoc-sample-exclusion.txt" \
-cov_all_continuous \
#-exclude_snps "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/info/GWAS3_new/$chr-fail.list" \
#-range "$start"-"$end"

#conditional
cond_snp="2:187502846_T_C"
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5 \
-data /lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.gen /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample \
-o "/lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.cond_on.${cond_snp}" \
-frequentist 1 \
-method score \
-pheno bin1 \
-exclude_samples "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/gwas3-${assoc}-assoc-sample-exclusion.txt" \
-cov_all_continuous \
-condition_on "$cond_snp" \
#-exclude_snps "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/info/GWAS3_new/$chr-fail.list" \
#-range "$start"-"$end"

#conditional 2
cond_snp="7:20577298_G_A"
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5 \
-data /lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.gen /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample \
-o "/lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.cond_on.${cond_snp}" \
-frequentist 1 \
-method score \
-pheno bin1 \
-exclude_samples "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/gwas3-${assoc}-assoc-sample-exclusion.txt" \
-cov_all_continuous \
-condition_on "$cond_snp" \
#-exclude_snps "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/info/GWAS3_new/$chr-fail.list" \
#-range "$start"-"$end"


#unconditional - general model
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5 \
-data /lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.gen /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample \
-o "/lustre/scratch115/teams/anderson/ibd_conditional/ITGAV_ITGB8.uncond.general" \
-frequentist 4 \
-method score \
-pheno bin1 \
-exclude_samples "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/gwas3-$assoc-assoc-sample-exclusion.txt" \
-cov_all_continuous \
#-exclude_snps "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/info/GWAS3_new/$chr-fail.list" \
#-range "$start"-"$end"

