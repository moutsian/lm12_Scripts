
#!/bin/bash
assoc=$1
chr=$2
start=$3
end=$4
cond_snp_nosp="16_50763778_50756540_50745926_50816078_50762771_50694011_50746199"
#/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5 \
#-data /lustre/scratch113/projects/crohns/RELEASE/v1/GWAS3/$chr.gen.gz /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample \
#-o "/lustre/scratch115/teams/anderson/ibd_conditional/gwas3.$chr.$start.uncond" \
#-frequentist 1 \
#-method score \
#-pheno bin1 \
#-exclude_samples "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/gwas3-$assoc-assoc-sample-exclusion.txt" \
#-cov_all_continuous \
#-exclude_snps "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/info/GWAS3_new/$chr-fail.list" \
#-range "$start"-"$end"

#conditional
/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5 \
-data /lustre/scratch113/projects/crohns/RELEASE/v1/GWAS3/$chr.gen.gz /lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/GWAS3.ibd.sample \
-o "/lustre/scratch115/teams/anderson/ibd_conditional/gwas3.$chr.$assoc.$start.$end.cond_on.${cond_snp}.fast" \
-frequentist 1 \
-method score \
-pheno bin1 \
-exclude_samples "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc_NEW/GWAS3/refs/gwas3-${assoc}-assoc-sample-exclusion.txt" \
-cov_all_continuous \
-condition_on "16:50763778_G_GC" "16:50756540_G_C" "16:50745926_C_T" "16:50816078_T_A" "16:50762771_T_C" "16:50694011_A_G" "16:50746199_G_A"  \
-exclude_snps "/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/info/GWAS3_new/$chr-fail.list" \
-range "$start"-"$end"
