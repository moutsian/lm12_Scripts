#!/usr/bin/bash

## note that in contrast with Javi's scripts, this works for both UC and IBD (just specify the trait parameter)

trait=$1
locus=$2
utrait=${trait^^} #upper case
ltrait=${trait,,} # lower case

lociFile="/lustre/scratch110/sanger/lm12/ibd/LM_list_to_condition_on_${utrait}_with_updated_IDs.txt"
conditional_res_file="/lustre/scratch110/sanger/lm12/ibd/meta_output/${utrait}.1.meta.tophits.txt"
round=1

study="GWAS2"
Javi_dir="/lustre/scratch113/teams/anderson/users/jga/001_Projects/Project2_FineMappingSeqData"
Yang_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc/ALL"
Data_dir="/lustre/scratch113/projects/crohns/RELEASE/v2/${study}"
Sample_dir="/lustre/scratch113/projects/crohns/2015jan20/GWAS_imputation/assoc"

results_dir="/lustre/scratch110/sanger/lm12/ibd/snptest_output/reverse/${study}"
log_dir="/lustre/scratch110/sanger/lm12/ibd/logs"
echo "$results_dir"
#mkdir "$results_dir"

while IFS=$'\t' read loci;
do
        teas=($loci)  #### all i need is to make a file containing these information. I need to take the output from the top association files  (we generated in step1 ) and the additional SNPs to condition on 
        loci=${teas[0]} 
        chr=${teas[1]}
        pos=${teas[2]}
        startLoci=${teas[3]}
        endLoci=${teas[4]}
		topSNP=${teas[6]} #the array without its first 5 elements.
		#echo $teas
#echo "$loci $topSNP"

if [ "$locus" == "$loci" ]; then
#SNPTEST
echo "Running for $loci"
bsub -q normal -J $loci.$utrait.$study -R"select[mem>3000] rusage[mem=3000]" -M3000 -o $results_dir/$study.$utrait.$loci.$chr.$startLoci.$endLoci.reverse.out -e $log_dir/$study.$utrait.$loci.$chr.$startLoci.$endLoci.reverse.err 	/nfs/team152/javi/ToolBox/snptest_v2.5_linux_x86_64_static/snptest_v2.5  -data $Data_dir/$chr.gen.gz  $Sample_dir/GWAS2_with_pcs.sample -o $results_dir/$study.$utrait.$loci.$chr.$startLoci.$endLoci.$round.reverse.assoc -frequentist 1 -method score -pheno bin1 -cov_all_continuous -condition_on $topSNP -range $startLoci-$endLoci -exclude_snps $Javi_dir/SNPsPassFilters/$chr.$ltrait.snps.NoPassFilters.OnlyIDs.txt 
fi
done < $conditional_res_file

