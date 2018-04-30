# Pipeline for exome sequencing data filtering and analysis for the TILI cohort, using bcftools and vcftools
# Dec 2017 - iteration 3, with the addition of the other TILI exome seq data as controls. 

exomedir="/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/"
bcftoolsdir="/software/team152/bcftools/bcftools1.6/bcftools/"
basename=${exomedir}"tili.i3.pb"
##########################################################################################################################################################
#### 	STEP1: Remove variants which don't have PASS for GATK VQSR. 																				  ####
#### 	Remove monomorphic variants. Keep only biallelic variants. Also keep only the samples in the list provided by Exeter.						  ####
##########################################################################################################################################################

# Remove special CASES (note that only 5 of them are actually sequenced): 
#PRED3825
#PRED5022 - NS
#PRED4004
#PRED4374
#PRED4292
#PRED3715
#PRED5008 - NS

# Iteration 3 addition: Include additional samples to keep as controls from Gareth's list, send in an email on the  7th of Dec. Note that 
# this list contains people a large number of people with other thiopurine adverse effects (e.g. TIM, THR).  
# We may want to remove these, but they constitute about half of the total list.
cat Lukas20171207_additional_exomeseq_ctrls.txt |grep -v VCF_ID|awk '{print $1}' > list_of_samples_to_keep.iteration3.txt
cat list_of_samples_to_keep.both.nospecial.txt >> list_of_samples_to_keep.iteration3.txt
cat list_of_samples_to_keep.iteration3.txt |sort|uniq> list_of_samples_to_keep.iteration3.unique.txt

bsub -q normal -J vcf_i3 -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}vcf_i3.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}list_of_samples_to_keep.iteration3.unique.txt --force-samples -O v \
 -o ${basename}.vcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}v31_ahmad.vcf.gz

# Make sure you are using the latest bcftools version (1.6 or more)

##########################################################################################################################################################
#### 	STEP2: Set genotypes with GQ<30 to missing and also filter by allelic imbalance																  ####
##########################################################################################################################################################

#**** ! NOTE THAT THE GQ genotype filter  DOES NOT WORK PROPERLY when outputting it as bcf file (at least not with vcftools 0.1.14)
# Just output it as vcf.
#vcftools --vcf ${basename}.vcf  --minGQ 30  --recode --recode-INFO-all --out ${basename}.gq30 #used in i2 version
${bcftoolsdir}./bcftools +setGT --threads 0 -Oz -o ${basename}.gq30.vcf.gz ${basename}.vcf -- -n. -t q -i 'FMT/GQ<30'
# Filter by allelic imbalance - then missingness, case-control differential missingness, then deviations from HWE
# Try the filtering by allele-balance that Hilary uses (note that Hilary is using a cutoff of p=1e-03)
${bcftoolsdir}./bcftools +setGT --threads 0 -Oz -o ${basename}.gq30.ab1e03.vcf.gz  ${basename}.gq30.vcf.gz --  -n . -t 'b:AD<1e-03'

##########################################################################################################################################################
#### 	STEP3: Filter variants for missingness																										  ####
##########################################################################################################################################################

#filter variants for missingness
vcftools --gzvcf ${basename}.gq30.ab1e03.vcf.gz  --max-missing 0.9 --recode --recode-INFO-all --out ${basename}.gq30.ab1e03.vcf.miss10pc 
#After filtering, we kept 404348 out of a possible 557065 Sites  -> The vast majority of them are due to the GQ>=30 threshold. Without this, less than 35-40K variants would be filtered out.
vcftools --gzvcf ${basename}.gq30.ab1e03.vcf.gz  --missing-site --out ${basename}.gq30.ab1e03
vcftools --vcf ${basename}.gq30.ab1e03.vcf.miss10pc.recode.vcf  --missing-site --out ${basename}.gq30.ab1e03.miss10pc

#${bcftoolsdir}./bcftools +fill-AN-AC ${basename}.gq30.ab1e03.miss10pc.vcf

mv ${basename}.gq30.ab1e03.vcf.miss10pc.recode.vcf ${basename}.gq30.ab1e03.miss10pc.vcf
#replace missing IDs '.' with a name containing the chromosome and position.
${bcftoolsdir}./bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ${basename}.gq30.ab1e03.miss10pc.vcf -Ob -o ${basename}.gq30.ab1e03.miss10.rs1.bcf

#now also update IDs to rsIDs where possible (from 1KG) 
#the following works, but make sure that the -a file is tab-delimited, bgzipped and tabixed
tabix -s1 -b2 -e2 chr${i}.1KG_b37.alleles #see tabixme.sh
#we need to use bcftools annotate across all chromosomes, so i have done so in a small script:
sh ${exomedir}update_rsIDs.sh

#remove duplicate (in terms of position) variants here because they seem to give trouble downstream with multiallelic SNPs
${bcftoolsdir}./bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER\n' ${basename}.gq30.ab1e03.miss10.rsID.bcf > ${basename}.gq30.ab1e03.miss10.rsID.bcf.variants
python handle_duplicates_in_tili_exomeseq.py --all ${basename}.gq30.ab1e03.miss10.rsID.bcf.variants --out duplicate_entries_to_remove_from_python_script.i3.txt #new for i3 (from R. much faster now)
#note that below we not just remove duplicate entries but also the sex chroms
${bcftoolsdir}./bcftools index ${basename}.gq30.ab1e03.miss10.rsID.bcf
${bcftoolsdir}./bcftools view -T ^${exomeseq}duplicate_entries_to_remove_from_python_script.i3.txt -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
-Ob -o ${basename}.gq30.ab1e03.miss10.rsID.nodup.bcf ${basename}.gq30.ab1e03.miss10.rsID.bcf

#now extract the cases and controls, to then get HWE stats separately
bsub -q normal -J vcfilt -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}vcfextract_ctrls.log ${bcftoolsdir}./bcftools view \
-S ^${exomedir}list_of_samples_to_keep.cases.txt --force-samples -O v -o ${basename}.gq30.ab1e03.miss10.rsID.ctrls.vcf --threads 4 \
${basename}.gq30.ab1e03.miss10.rsID.nodup.bcf

bsub -q normal -J vcfilt -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}vcfextract_cases.log ${bcftoolsdir}./bcftools view \
-S ${exomedir}list_of_samples_to_keep.cases.txt --force-samples -O v -o ${basename}.gq30.ab1e03.miss10.rsID.cases.vcf --threads 4 \
${basename}.gq30.ab1e03.miss10.rsID.nodup.bcf

#get HWE stats
vcftools --vcf ${basename}.gq30.ab1e03.miss10.rsID.cases.vcf  --hardy --out ${basename}.gq30.ab1e03.miss10.rsID.cases
vcftools --vcf ${basename}.gq30.ab1e03.miss10.rsID.ctrls.vcf  --hardy --out ${basename}.gq30.ab1e03.miss10.rsID.ctrls

cat ${basename}.gq30.ab1e03.miss10.rsID.ctrls.hwe|awk '{if($6<1e-08){print $1"\t"$2}}'> ${exomedir}snps_to_remove_after_hwe_check.ctrls.i3.txt 
cat ${basename}.gq30.ab1e03.miss10.rsID.cases.hwe|awk '{if($6<1e-08){print $1"\t"$2}}'> ${exomedir}snps_to_remove_after_hwe_check.cases.i3.txt 
cat ${exomedir}snps_to_remove_after_hwe_check.ctrls.i3.txt  ${exomedir}snps_to_remove_after_hwe_check.cases.i3.txt > ${exomedir}snps_to_remove_after_hwe_check.i3.txt

#before removing due to hwe deviations, also check for differential missingness, then remove together.
#here convert to plink and save the case/control information, then merge and give in to plink with -test-missing command 
plink --vcf ${basename}.gq30.ab1e03.miss10.rsID.cases.vcf --make-bed --out ${basename}.gq30.ab1e03.miss10.rsID.cases
plink --vcf ${basename}.gq30.ab1e03.miss10.rsID.ctrls.vcf --make-bed --out ${basename}.gq30.ab1e03.miss10.rsID.ctrls
cat ${basename}.gq30.ab1e03.miss10.rsID.ctrls.fam|awk '{print $1,$2,$3,$4,$5,1}' > ${basename}.gq30.ab1e03.miss10.rsID.ctrls.corr.fam
cat ${basename}.gq30.ab1e03.miss10.rsID.cases.fam|awk '{print $1,$2,$3,$4,$5,2}' > ${basename}.gq30.ab1e03.miss10.rsID.cases.corr.fam
mv ${basename}.gq30.ab1e03.miss10.rsID.ctrls.corr.fam  ${basename}.gq30.ab1e03.miss10.rsID.ctrls.fam
mv ${basename}.gq30.ab1e03.miss10.rsID.cases.corr.fam  ${basename}.gq30.ab1e03.miss10.rsID.cases.fam
plink --allow-no-sex --bfile ${basename}.gq30.ab1e03.miss10.rsID.ctrls --bmerge \
${basename}.gq30.ab1e03.miss10.rsID.cases --make-bed --out ${basename}.gq30.ab1e03.miss10.rsID
#now we can run the differential missingness test
plink --allow-no-sex --bfile  ${basename}.gq30.ab1e03.miss10.rsID --test-missing midp --out ${basename}.gq30.ab1e03.miss10.rsID
cat ${basename}.gq30.ab1e03.miss10.rsID.missing|awk '{if($5<1e-05){print $2}}'> ${exomeseq}snps_to_remove_due_to_diffmiss.1e-05.i3.txt

#now get chrom and position from the rsIDs, so then bcftools is used to extract the relevant variants 
python from_rsID_to_varpos.py -V ${basename}.gq30.ab1e03.miss10.rsID.bcf.variants --rs ${exomeseq}snps_to_remove_due_to_diffmiss.1e-05.i3.txt --out snps_to_remove_due_to_diffmiss.varpos.i3.txt
#now remove
cat ${exomedir}snps_to_remove_due_to_diffmiss.varpos.i3.txt  ${exomedir}snps_to_remove_after_hwe_check.i3.txt > ${exomedir}snps_to_remove_after_hwe_and_diffmiss_check.i3.txt

#the following should work (Petr's suggestion):
${bcftoolsdir}./bcftools view -T ^${exomeseq}snps_to_remove_after_hwe_and_diffmiss_check.i3.txt -Ob -o ${exomedir}tili.i3.site_QC.bcf \
${exomedir}tili.i3.pb.gq30.ab1e03.miss10.rsID.nodup.bcf

#list of the remaining variants:
${bcftoolsdir}./bcftools view ${exomedir}tili.i3.site_QC.bcf |grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${exomedir}tili.i3.site_QC.variants

##################################################################################################################################################################################
############################################# END OF FIRST ROUND OF VARIANT QC 						##############################################################################
##################################################################################################################################################################################

##################################################################################################################################################################################
############################################# 			SAMPLE QC1 				##################################################################################################
##################################################################################################################################################################################

basename2=${exomedir}"tili.i3.site_QC"
${bcftoolsdir}./bcftools view ${basename2}.bcf -Ov -o ${basename2}.vcf
## At this stage we have a vcf file with standard site-filtering (VQSR, GQ, missingness, HWE deviations) applied to it. Now we can do some sample filtering too
# a) Start by looking into individuals with high missingness  
vcftools --vcf ${basename2}.vcf --missing-indv --out ${basename2}

# b) Heterozygosity 
vcftools --vcf ${basename2}.vcf --het --out ${basename2}
# b2) heterozygosity at low frequency variants
vcftools --vcf ${basename2}.vcf --het --out ${basename2}.lowfreq --max-maf 0.01
#c) number of singletons
vcftools --vcf ${basename2}.vcf  --indv-freq-burden --out ${basename2}
Rscript het_and_miss_filter.R tili.i3.site_QC
#note that using 3 sd for number of singletons seems to only remove samples that would have already been removed due to heterozygosity and/or high missingness.

# CAREFUL WITH THE FORMAT OF THE samples_to_remove file, or bcftools will have issues.
${bcftoolsdir}./bcftools view -S ^${basename2}.samples_to_remove.txt --force-samples -O v -o ${basename2}.sample_qc1.vcf \
--threads 4 ${basename2}.vcf
#check this worked OK
${bcftoolsdir}./bcftools query -l  ${basename2}.vcf | wc -l
${bcftoolsdir}./bcftools query -l  ${basename2}.sample_qc1.vcf | wc -l #number of samples reduced compared to previous file.

# d) Relatedness

#relatedness, measure 2, only using variants with MAF>0.5% 
vcftools --vcf ${basename2}.sample_qc1.vcf --relatedness2 --maf 0.005 --out ${basename2}.sample_qc1

#get pairs of related individuals
cat ${basename2}.sample_qc1.relatedness2|awk '{if($7>0.09375 && $1!=$2){print $0}}' > ${basename2}.sample_qc1.highrelatedness.txt
#now pepare file of IDs to remove based on which individual in each pair has more missing data.
python ${exomeseq}prep_high_ibs_samples_to_remove_exomeseq.py --miss  ${basename2}.imiss --ibs ${basename2}.sample_qc1.highrelatedness.txt --out ${exomeseq}high_ibs.it3.samples.toremove.txt

${bcftoolsdir}./bcftools view -S ^${exomedir}high_ibs.it3.samples.toremove.txt --force-samples -O v -o ${basename2}.sample_qc2.vcf --threads 4 ${basename2}.sample_qc1.vcf
${bcftoolsdir}./bcftools query -l  ${basename2}.sample_qc2.vcf | wc -l #number of samples reduced compared to previous file.

##################################################################################################################################################################################
############################################# 			PCA 				######################################################################################################
##################################################################################################################################################################################
basename3=${exomedir}"tili.i3.site_QC.sample_qc2"
DATADIR="/lustre/scratch115/projects/crohns/exome/TIH/"
plink --vcf ${basename3}.vcf --exclude ${DATADIR}aux_and_intermediate_files/high_ld.txt \
--indep-pairwise 1000 50 0.2 --maf 0.01 --chr 1-22 --out ${basename3}
plink --vcf ${basename3}.vcf --make-bed --out  ${basename3}
plink --allow-no-sex --bfile ${basename3} --extract ${basename3}.prune.in  --make-bed --out ${basename3}.pruned

# add case /ctrl status
python annotate_case_ctrl.py --ctrl tili.i3.pb.gq30.ab1e03.miss10.rsID.ctrls.fam --case tili.i3.pb.gq30.ab1e03.miss10.rsID.cases.fam -f tili.i3.site_QC.sample_qc2.pruned.fam
grep "\-9\b" ${basename3}.pruned.fam.updated  #to double check it worked
#now move it and continue
mv ${basename3}.pruned.fam.updated ${basename3}.pruned.fam
grep 'rs' ${basename3}.pruned.bim|awk '{print $2}' > ${exomedir}exomeseq.snps_to_extract_from_1KG.txt
##now calculate PCs using the pruned merged dataset we just created.
sh get_1KG_pc_exomeseq.sh "${exomedir}" "exomeseq.snps_to_extract_from_1KG.txt"
#merge into one file
plink --bfile ${exomedir}1KG.chr1.forPCA.exomeseq --merge-list ${exomedir}1KG.list --make-bed --out ${exomedir}1KG_for_PCA_allchr.i3

#Remove SNPs which have different alleles between 1KG and exome seq data.
python check_allele_agreement.py ${exomeseq}1KG_for_PCA_allchr.i3.bim ${basename3}.pruned.bim
plink --bfile ${exomeseq}1KG_for_PCA_allchr.i3 --exclude ${exomedir}snps_with_different_alleles.txt  --make-bed --out ${exomeseq}1KG_for_PCA_allchr.i3_ready
cat ${exomeseq}1KG_for_PCA_allchr.i3_ready.bim|awk '{print $2}'> ${exomeseq}snps_present_in_both.txt
plink --allow-no-sex --bfile ${basename3}.pruned --bmerge ${exomeseq}1KG_for_PCA_allchr.i3_ready --make-bed \
--extract ${exomeseq}snps_present_in_both.txt  --out ${exomeseq}TILI_1KG_forPCA.i3
#now run PCs - **REMEMBER** that smartpca misbehaves (ignores samples) when they have no case/ctrl status assigned, so do this for 1KG prior to running
python annotate_case_ctrl.py --ctrl tili.i3.pb.gq30.ab1e03.miss10.rsID.ctrls.fam --case tili.i3.pb.gq30.ab1e03.miss10.rsID.cases.fam -f TILI_1KG_forPCA.i3.fam
mv TILI_1KG_forPCA.i3.fam.updated TILI_1KG_forPCA.i3.fam
plink --allow-no-sex --bfile TILI_1KG_forPCA.i3 --recode --out TILI_1KG_forPCA.i3
cat   TILI_1KG_forPCA.i3.ped |awk '{print $1,$2,$3,$4,$5,$6}' > TILI_1KG_forPCA.i3.pedind
smartpca -p TILI_1KG_exomeseq.par > TILI_1KG_exomeseq.Sout  #note that this isn't automated and I should prepare an automated generation of this file at some point.

##now visually observe and make removal lists in R using pca_TILI_with1KG_exomeseq.R
cat ${exomeseq}ctrls_to_remove_due_to_PCA_exomeseq.i3.txt ${exomeseq}cases_to_remove_due_to_PCA_exomeseq.i3.txt |awk '{print $1}' > ${exomeseq}all_to_remove_due_to_PCA_exomeseq.i3.txt
${bcftoolsdir}./bcftools view -S ^${exomedir}all_to_remove_due_to_PCA_exomeseq.i3.txt --force-samples -O b -o ${basename2}.sample_QC.j18.bcf --threads 4 ${basename3}.bcf

#############################################################################################################################################################
#############################################################################################################################################################
################################## 			i3 iteration addition 					#########################################################################
## 	Here, because of the large number of samples being removed, I rerun the snp filtering process starting from the filtered set of samples.			   ##
##################################							   						#########################################################################
#############################################################################################################################################################
#############################################################################################################################################################
cat ${exomeseq}all_to_remove_due_to_PCA_exomeseq.i3.txt ${exomeseq}high_ibs.it3.samples.toremove.txt tili.i3.site_QC.samples_to_remove.txt > tili.i3.allsamples_to_remove.txt
${bcftoolsdir}./bcftools view -S ^${exomedir}tili.i3.allsamples_to_remove.txt --force-samples -O b -o ${basename}.sample_QC.j18.bcf ${basename}.vcf
# I did this part through a script, var_re_qc_after_sample_qc.sh. Check there if there are any issues 
sh ${exomedir}var_re_qc_after_sample_qc.sh

##################################################################################################################################################################################
############################################# END OF SECOND ROUND OF VARIANT QC 						##############################################################################
##################################################################################################################################################################################

#After the above is done, we can continue with annotation and the final round of variant QC, that of binomial test of frequencies vs 1kg and gnomAD.
#Then we can move to the analysis 
#Output file from the above process is tili.i3.site_QC.sample_QC.bcf
vcftools --vcf ${exomedir}tili.i3.site_QC.sample_QC.vcf  --missing-site --out ${exomedir}tili.i3.site_QC.sample_QC #just a check

#perl INSTALL.pl --CACHEDIR /software/team152/ensembl-vep/cache/ -s homo_sapiens_merged --NO_TEST #this gives 23/39 failed tests (26/858 failed subtests)
##I have downloaded the cache for GRCh37 in /software/team/ensembl-vep/cache/ To use it, since it's a merged Ensembl & Refseq, I need to add --merged.
****edw 11/01

${bcftoolsdir}./bcftools view tili.i3.site_QC.sample_QC.bcf -Ov -o tili.i3.site_QC.sample_QC.vcf
###here we should do annotation and then  filter using the binomial test.
###Only then we should move into sample qc.
 module add hgi/systems/R/3.1.2
${bcftoolsdir}./bcftools view ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.cases.vcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.cases.variants
${bcftoolsdir}./bcftools view ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.ctrls.vcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.ctrls.variants
#running vep on the farm (worked fine)
bsub -q long -J veprunXmas -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep_siteQC_sampleQC.log \
 /software/team152/ensembl-vep/./vep -i ${exomedir}tili.i3.site_QC.sample_QC.vcf  --offline  \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  ${exomedir}tili.i3.site_QC.sample_QC.vep.stats --use_given_ref \
--output_file ${exomedir}tili.i3.site_QC.sample_QC.vcf.annot --force_overwrite \
--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

# the above produced a separate file, now I will work on incorporating some of these info into the vcf
##note that the fields that are added to the vcf are specified explicitly. 
##note that the *.annot file can be checked for more information
bsub -q long -J veprun2Xmas -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep2_siteQC_sampleQC.log \
/software/team152/ensembl-vep/./vep -i ${exomedir}tili.i3.site_QC.sample_QC.vcf --cache --offline \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --use_given_ref --vcf \
--output_file ${exomedir}tili.i3.site_QC.sample_QC.annot.vcf --force_overwrite --per_gene \
--sift b -poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz --symbol \
--fields Consequence,SYMBOL,CANONICAL,SIFT,PolyPhen,CADD_PHRED,EUR_AF,gnomAD_AMR_AF  

***the below are to be run after the annotation is done (think about the exact filter in the awk line, this is currently quite strict and seems to be removing ca. 2.5% of the post-QC variants)
python2 freq_binom_test_v3.py --ann tili.i3.site_QC.sample_QC.vcf.annot --var tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.ctrls.variants --out tili_binom_test_results.i3.site_QC.sample_QC.txt
cat tili_binom_test_results.i3.site_QC.sample_QC.txt|awk '{if( ($8<1e-05 && $6>0 && $6<1) || ($8=="variant_not_in_gnomAD" && $9=="variant_not_in_1KG" && ($4/$5)>=0.01 ) ){print $1"\t"$2}}'> vars_to_remove_binom_test_1e-05.gnomAD.i3.txt
cat tili_binom_test_results.i3.site_QC.sample_QC.txt|awk '{if( ($9<1e-09 && $7>0 && $7<1)  ){print $1"\t"$2}}'> vars_to_remove_binom_test_1e-05.1KG.i3.txt
cat vars_to_remove_binom_test_1e-05.gnomAD.i3.txt vars_to_remove_binom_test_1e-05.1KG.i3.txt | sort -n -k 1,1 -k 2,2 |uniq > vars_to_remove_binom_test_1e-05.i3.correct.txt
${bcftoolsdir}./bcftools view -T ^${exomeseq}vars_to_remove_binom_test_1e-05.i3.correct.txt -Ov -o ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf \
${exomedir}tili.i3.site_QC.sample_QC.annot.vcf
##################################################################################################################################################################################
############################################# END OF VARIANT ANNOTATION 						######################################################################
##################################################################################################################################################################################
basename5=${exomedir}"tili.i3.site_QCplus_binom.sample_QC.annot"

#GETTING PCs for filtered dataset 
${bcftoolsdir}./bcftools query -l ${basename5}.vcf|wc -l
#now get the PCs from this dataset to use as covariates
plink --vcf ${basename5}.vcf --exclude /lustre/scratch115/projects/crohns/exome/TIH/aux_and_intermediate_files/high_ld.txt \
 --maf 0.01 --chr 1-22 --indep-pairwise 1000 50 0.2 --out ${basename5}
plink --allow-no-sex --vcf ${basename5}.vcf --extract ${basename5}.prune.in --make-bed --out ${basename5}.pruned
python annotate_case_ctrl.py --ctrl tili.i3.pb.sample_QC.gq30.ab1e03.miss10.rsID.ctrls.fam --case tili.i3.pb.sample_QC.gq30.ab1e03.miss10.rsID.cases.fam -f tili.i3.site_QCplus_binom.sample_QC.annot.pruned.fam
mv ${basename5}.pruned.fam.updated ${basename5}.pruned.fam
plink --bfile ${basename5}.pruned --recode --out ${basename5}.pruned
cat   ${basename5}.pruned.ped |awk '{print $1,$2,$3,$4,$5,$6}' > ${basename5}.pruned.pedind
smartpca -p tili.i3.postQC.par > tili.i3.postQC.pca.Sout

**** here we should have all we need for downstream analysis.
########################################
## Single variant test using SNPTEST2 ##
########################################
#tili.i3.site_QCplus_binom.sample_QC.annot.vcf
#SNPTEST2,2 PCs as covariates (Only 2 PCs with p<0.05 in Tracy Widom test)
plink --vcf ${basename5}.vcf --recode oxford --out ${basename5}
#now add covariates and case/ctrl info to the .sample file
python2 prepare_sample_for_snptest.py -e tili.i3.site_QCplus_binom.sample_QC.annot.pruned.30.evec  -i tili.i3.site_QCplus_binom.sample_QC.annot.pruned.fam -s tili.i3.site_QCplus_binom.sample_QC.annot.sample
#now run snptest using the updated sample file
snptest -data  ${basename5}.gen ${basename5}.sample.updated \
-o ${exomedir}output/tili.i3.output.score.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method score 
#run snptest2 with em as well
snptest -data  ${basename5}.gen ${basename5}.sample.updated \
-o ${exomedir}output/tili.i3.output.em.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method em 


###############################		trying GEMMA for association testing with univariate LMMs		###########################################################################

Skipping it for i3 (for now at least)
# plink --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf --make-bed --out ${exomedir}tili.i2.qc.sample_qc3.rsID.annot
# python ${exomedir}annotate_case_ctrl.py --ctrl tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.fam --case tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.fam -f tili.i2.qc.sample_qc3.rsID.annot.fam
# mv tili.i2.qc.sample_qc3.rsID.annot.fam.updated tili.i2.qc.sample_qc3.rsID.annot.fam
# /software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc.sample_qc3.rsID.annot -gk 2 -o ${exomedir}tili.i2.qc.sample_qc3.rsID.annot #generate standardised relatedness matrix

# /software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc.sample_qc3.rsID.annot -k ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.sXX.txt \
# -eigen -o tili.i2.qc.sample_qc3.rsID.annot

# /software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc.sample_qc3.rsID.annot \
# -d ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.eigenD.txt \
# -u ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.eigenU.txt \
# -hwe 0.00000001 -maf 0.005 -lmm 4 -o tili.i2.qc.sample_qc3.rsID.annot.gemma_output


###################################################################################################################################
################################################## Gene-based tests using EPACTS ##################################################
###################################################################################################################################
 basename6=${exomedir}"tili.i3.site_QCplus_binom.sample_QC"

#For file preparation look at 
prepare_group_file_for_EPACTS.sh

#after done with the above step, prepare a ped file 
python prepare_sample_for_epacts.py --sample tili.i3.site_QCplus_binom.sample_QC.annot.sample.updated
EPACTS_DIR="/software/team152/EPACTS/" 
CHR=21
#vcf file must be tabixed
bgzip tili.i3.site_QCplus_binom.sample_QC.annot.vcf
tabix tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient.log \
${EPACTS_DIR}bin/epacts group --vcf ${basename6}.annot.vcf.gz \
  --groupf ${exomedir}tili.i3.site_QC.sample_QC.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/tili.i3.site_QC.sample_QC.${CHR}.PD.lenient.0.01.skato.nocov \
  --ped ${basename6}.annot.sample.updated.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${basename6}.annot.vcf.gz \
  --groupf ${exomedir}tili.i3.site_QC.sample_QC.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/tili.i3.site_QC.sample_QC.${CHR}.PD.strict.0.01.skato.nocov \
  --ped ${basename6}.annot.sample.updated.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.fc.log \
${EPACTS_DIR}bin/epacts group --vcf ${basename6}.annot.vcf.gz \
  --groupf ${exomedir}tili.i3.site_QC.sample_QC.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/tili.i3.site_QC.sample_QC.${CHR}.fc.skato.nocov \
  --ped ${basename6}.annot.sample.updated.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done

###		also run with covariates		###
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient.log \
${EPACTS_DIR}bin/epacts group --vcf ${basename6}.annot.vcf.gz \
  --groupf ${exomedir}tili.i3.site_QC.sample_QC.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/tili.i3.site_QC.sample_QC.${CHR}.PD.lenient.0.01.skato.withcov \
  --ped ${basename6}.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 \
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${basename6}.annot.vcf.gz \
  --groupf ${exomedir}tili.i3.site_QC.sample_QC.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/tili.i3.site_QC.sample_QC.${CHR}.PD.strict.0.01.skato.withcov \
  --ped ${basename6}.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 \
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.i3.fc.log \
${EPACTS_DIR}bin/epacts group --vcf ${basename6}.annot.vcf.gz \
  --groupf ${exomedir}tili.i3.site_QC.sample_QC.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/tili.i3.site_QC.sample_QC.${CHR}.fc.skato.withcov \
  --ped ${basename6}.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 \
   --pheno DISEASE --test skat --skat-o --run 2
done


##############		Try emmaxVT		##############
# basename6=${exomedir}"tili.i3.site_QCplus_binom.sample_QC"

${EPACTS_DIR}bin/epacts make-kin  --vcf ${basename6}.annot.vcf.gz \
--ped  ${basename6}.annot.sample.updated.for_epacts.ped  --min-maf 0.01 -min-callrate 0.95  \
--out ${exomedir}output/tili.i3.site_QCplus_binom.sample_QC.for_epacts.kin --run 2

annot="fc"
annot="PD.lenient"
annot="PD.strict"

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.${annot}.i3.emmaxVT.log \
${EPACTS_DIR}bin/epacts group --vcf ${basename6}.annot.vcf.gz \
  --groupf ${exomedir}tili.i3.site_QC.sample_QC.${CHR}.annot.${annot}.genes.epacts \
  --kin ${exomedir}output/tili.i3.site_QCplus_binom.sample_QC.for_epacts.kin \
  --out ${exomedir}output/tili.i3.${CHR}.${annot}.0.01.emmaxVT \
  --ped ${basename6}.annot.sample.updated.for_epacts.ped --max-maf 0.01\
   --pheno DISEASE --test emmaxVT --run 2
done


###################***************************************************************************************************************************************************################
######################################################################################################################################################################################
##########################################			 now have a look at subgroups too 			###################################################################################### 
######################################################################################################################################################################################
###################***************************************************************************************************************************************************################

#cholestatic and mixed vs controls 
#a) get controls

exomedir="/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/"
bcftoolsdir="/software/team152/bcftools/bcftools1.6/bcftools/"
cat hepmix.txt cholmix.txt > hepcholmix.txt
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}bcftools.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ^${exomedir}hepcholmix.txt --force-samples -O b \
 -o ${exomedir}tili.i3.site_QCplus_binom.sample_QC.ctrls.bcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz

#b) get cholestatic and mixed cases: (rememeber to remove special cases)
#special cases:
#PRED3825
#PRED5022 - NS
#PRED4004
#PRED4374
#PRED4292
#PRED3715
#PRED5008 - NS

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}bcftoolscholmix.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}cholmix.txt --force-samples -O b \
 -o ${exomedir}tili.i3.site_QCplus_binom.sample_QC.cholmix.bcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz
 
#c) get hepatocellular and mixed cases (again, remember to remove special cases):

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}bcftoolscholmix.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}hepmix.txt --force-samples -O b \
 -o ${exomedir}tili.i3.site_QCplus_binom.sample_QC.hepmix.bcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz

basename6=${exomedir}"tili.i3.site_QCplus_binom.sample_QC"
 
${bcftoolsdir}./bcftools query -l  ${basename6}.hepmix.bcf | wc -l
${bcftoolsdir}./bcftools query -l  ${basename6}.cholmix.bcf | wc -l
${bcftoolsdir}./bcftools query -l  ${basename6}.ctrls.bcf | wc -l

 #now merge each of the subgroups with ctrls, to prepare the files for the analyses. 
 # Ensure that only sites present in the subgroup and in the controls are kept.
${bcftoolsdir}./bcftools view ${basename6}.ctrls.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${basename6}.ctrls.variants
${bcftoolsdir}./bcftools view ${basename6}.cholmix.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${basename6}.cholmix.variants
${bcftoolsdir}./bcftools view ${basename6}.hepmix.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${basename6}.hepmix.variants

python2 list_shared_vars.py -1 ${basename6}.cholmix.variants -2 ${basename6}.ctrls.variants --out shared_vars_between_cholmix_and_ctrls.i3.txt
python2 list_shared_vars.py -1 ${basename6}.hepmix.variants -2 ${basename6}.ctrls.variants --out shared_vars_between_hepmix_and_ctrls.i3.txt

#bgzip  -f ${exomedir}tili.i2.qc.sample_qc3.rsID.ctrls.vcf 
#bgzip  -f tili.i2.qc.sample_qc3.rsID.cholmix.vcf
#bgzip  -f tili.i2.qc.sample_qc3.rsID.hepmix.vcf
#tabix tili.i2.qc.sample_qc3.rsID.hepmix.vcf.gz
#tabix tili.i2.qc.sample_qc3.rsID.cholmix.vcf.gz
#tabix tili.i2.qc.sample_qc3.rsID.ctrls.vcf.gz
${bcftoolsdir}./bcftools index  ${basename6}.cholmix.bcf
${bcftoolsdir}./bcftools index  ${basename6}.ctrls.bcf
${bcftoolsdir}./bcftools index  ${basename6}.hepmix.bcf

bcftools merge  ${basename6}.ctrls.bcf ${basename6}.cholmix.bcf -O b -o ${basename6}.cholmix_vs_ctrls.bcf --threads 3 --force-samples
bcftools merge  ${basename6}.ctrls.bcf ${basename6}.hepmix.bcf  -O b -o ${basename6}.hepmix_vs_ctrls.bcf  --threads 3 --force-samples
#now do the site filtering
${bcftoolsdir}./bcftools view -T ${exomeseq}shared_vars_between_cholmix_and_ctrls.i3.txt -Ov -o ${basename6}.cholmix_vs_ctrls.qc.vcf  \
${basename6}.cholmix_vs_ctrls.bcf
${bcftoolsdir}./bcftools view -T ${exomeseq}shared_vars_between_hepmix_and_ctrls.i3.txt -Ov -o ${basename6}.hepmix_vs_ctrls.qc.vcf  \
${basename6}.hepmix_vs_ctrls.bcf


***edw Paraskevi 19/01 ****
plink --vcf ${basename6}.hepmix_vs_ctrls.qc.vcf --recode oxford --out ${basename6}.hepmix_vs_ctrls.qc
plink --vcf ${basename6}.cholmix_vs_ctrls.qc.vcf --recode oxford --out ${basename6}.cholmix_vs_ctrls.qc
#now add covariates and case/ctrl info to the .sample file
python2 prepare_sample_for_snptest.py -e tili.i3.site_QCplus_binom.sample_QC.annot.pruned.30.evec  -i tili.i3.site_QCplus_binom.sample_QC.annot.pruned.fam -s tili.i3.site_QCplus_binom.sample_QC.hepmix_vs_ctrls.qc.sample
python2 prepare_sample_for_snptest.py -e tili.i3.site_QCplus_binom.sample_QC.annot.pruned.30.evec  -i tili.i3.site_QCplus_binom.sample_QC.annot.pruned.fam -s tili.i3.site_QCplus_binom.sample_QC.cholmix_vs_ctrls.qc.sample
#now run snptest using the updated sample file
snptest -data  ${basename6}.hepmix_vs_ctrls.qc.gen ${basename6}.hepmix_vs_ctrls.qc.sample.updated \
-o ${exomedir}output/tili.i3.hepmix_vs_ctrls.score.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method score 
snptest -data  ${basename6}.cholmix_vs_ctrls.qc.gen ${basename6}.cholmix_vs_ctrls.qc.sample.updated \
-o ${exomedir}output/tili.i3.cholmix_vs_ctrls.score.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method score 


#after merging, continue with the individual analyses.
#SNPTEST2,3 PCs as covariates (Only 2 PCs with p<0.05 in Tracy Widom test)
plink --vcf ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.vcf  --recode oxford --out ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls
plink --vcf ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.vcf --recode oxford --out ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls
#now add covariates and case/ctrl info to the .sample file
python prepare_sample_for_snptest.py -e tili.i2.qc.sample_qc3.pruned.30.evec  -i tili.i2.qc.sample_qc3.pruned.fam -s tili.i2.qc3binom.hepmix_vs_ctrls.sample
python prepare_sample_for_snptest.py -e tili.i2.qc.sample_qc3.pruned.30.evec  -i tili.i2.qc.sample_qc3.pruned.fam -s tili.i2.qc3binom.cholmix_vs_ctrls.sample

#now run snptest using the updated sample file
snptest -data  ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.gen ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.sample.updated \
-o ${exomedir}output/tili.i2.qc3binom.hepmix_vs_ctrls.output.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method score 
snptest -data  ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.gen ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.sample.updated \
-o ${exomedir}output/tili.i2.qc3binom.cholmix_vs_ctrls.output.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method score 


###############################		 GEMMA for association testing with univariate LMMs		###########################################################################

plink --vcf ${exomedir}tili.i2.qc3.cholmix_vs_ctrls.vcf --make-bed --out ${exomedir}tili.i2.qc3.cholmix_vs_ctrls
python ${exomedir}annotate_case_ctrl.py --ctrl tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.fam --case tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.fam -f tili.i2.qc3.cholmix_vs_ctrls.fam
mv tili.i2.qc3.cholmix_vs_ctrls.fam.updated tili.i2.qc3.cholmix_vs_ctrls.fam
/software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc3.cholmix_vs_ctrls -gk 2 -o tili.i2.qc3.cholmix_vs_ctrls #generate standardised relatedness matrix

/software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc3.cholmix_vs_ctrls -k ${exomedir}output/tili.i2.qc3.cholmix_vs_ctrls.sXX.txt \
-eigen -o tili.i2.qc3.cholmix_vs_ctrls

/software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc3.cholmix_vs_ctrls \
-d ${exomedir}output/tili.i2.qc3.cholmix_vs_ctrls.eigenD.txt \
-u ${exomedir}output/tili.i2.qc3.cholmix_vs_ctrls.eigenU.txt \
-lmm 4 -o tili.i2.qc3.cholmix_vs_ctrls.0.1.gemma_output -hwe 0.00000001 -maf 0.0005 -miss 0.1

###############################		 Annotation		###########################################################################

module add hgi/systems/R/3.1.2
#running vep on the farm (worked fine)
bsub -q normal -J veprun -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep.log \
 /software/team152/ensembl-vep/./vep -i ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.vcf  --offline  \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.vep.stats --use_given_ref \
--output_file ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.vcf.annot --force_overwrite \
--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

# the above produced a separate file, now I will work on incorporating some of these info into the vcf
##note that the fields that are added to the vcf are specified explicitly. 
##note that the *.annot file can be checked for more information
bsub -q normal -J veprun2 -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep2.log \
/software/team152/ensembl-vep/./vep -i ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.vcf --cache --offline \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --use_given_ref --vcf \
--output_file ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.annot.vcf --force_overwrite --per_gene \
--sift b -poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz --symbol \
--fields Consequence,SYMBOL,CANONICAL,SIFT,PolyPhen,CADD_PHRED,EUR_AF,gnomAD_AMR_AF  



###################################################################################################################################
################################################## Gene-based tests using EPACTS ##################################################
###################################################################################################################################

#For file preparation look at 
prepare_group_file_for_EPACTS.sh



##############
#########################
#####################################
###############################################
#################################################################### QUICK CHECK IF BINOM IMPROVES GENE BASED RESULTS FOR THE FULL DATASET
###############################################################################################################################################################
###############################################
#####################################
#########################
##############

${bcftoolsdir}./bcftools view -T ^${exomeseq}vars_to_remove_binom_test_strict.txt -Ov -o ${exomedir}tili.i2.qc3binom_strict.sample_qc3.rsID.annot.vcf \
${exomedir}tili.i2.qc.sample_qc3.rsID.vcf
###################################################################################################################################
################################################## Gene-based tests using EPACTS ##################################################
###################################################################################################################################
#tili.i2.qc3binom_strict.sample_qc3.rsID.annot.vcf
#module add hgi/systems/R/3.1.2
#running vep on the farm (worked fine)
bsub -q long -J veprun -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep.log \
 /software/team152/ensembl-vep/./vep -i ${exomedir}tili.i2.qc3binom_strict.sample_qc3.rsID.annot.vcf  --offline  \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  ${exomedir}tili.i2.qc3binom_strict.sample_qc3.rsID.vep.stats --use_given_ref \
--output_file ${exomedir}tili.i2.qc3binom_strict.sample_qc3.rsID.vcf.annot --force_overwrite \
--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

1****edw**** me to strict *****1
#For file preparation look at 
prepare_group_file_for_EPACTS.sh

#after done with the above step, prepare a ped file 
python prepare_sample_for_epacts.py --sample tili.i2.qc.sample_qc3.rsID.annot.sample.updated
EPACTS_DIR="/software/team152/EPACTS/" 
CHR=21
#vcf file must be tabixed
bgzip tili.i2.qc3binom_strict.sample_qc3.rsID.annot.vcf
tabix tili.i2.qc3binom_strict.sample_qc3.rsID.annot.vcf.gz

tili.i2.qc3binom.sample_qc3.rsID.22.annot.fc.genes.epacts
###		also run with covariates		###
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict1.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.sample_qc3.rsID.${CHR}.PD.strict.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient1.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.sample_qc3.rsID.${CHR}.PD.lenient.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.fc1.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.sample_qc3.rsID.${CHR}.fc.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done


# check for AF between cases and controls
exomedir="/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/"
bcftoolsdir="/software/team152/bcftools/bcftools1.6/bcftools/"
${bcftoolsdir}./bcftools stats   --threads 4 ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.cases.vcf >  ${exomedir}tili.i3.site_QCplus_binom.sample_QC.cases.stats
${bcftoolsdir}./bcftools stats  --threads 4 ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.ctrls.vcf >  ${exomedir}tili.i3.site_QCplus_binom.sample_QC.ctrls.stats 

 ${bcftoolsdir}/misc/plot-vcfstats ${exomedir}tili.i3.site_QCplus_binom.sample_QC.ctrls.stats -p ${exomedir}/ctrl_plots/
 ${bcftoolsdir}/misc/plot-vcfstats ${exomedir}tili.i3.site_QCplus_binom.sample_QC.cases.stats -p ${exomedir}/case_plots/

vcftools --vcf ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.cases.vcf --indv-freq-burden --out tili.i3.site_QCplus_binom.sample_QC.cases
vcftools --vcf ${exomedir}tili.i3.pb.sample_QC.j18.gq30.ab1e03.miss10.rsID.ctrls.vcf --indv-freq-burden --out tili.i3.site_QCplus_binom.sample_QC.ctrls
vcftools --gzvcf ${exomedir}tili.i3.site_QCplus_binom.sample_QC.annot.vcf.gz --indv-freq-burden --out tili.i3.site_QCplus_binom.sample_QC