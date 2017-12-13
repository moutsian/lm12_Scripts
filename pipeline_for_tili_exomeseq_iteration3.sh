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

mv ${basename}.gq30.ab1e03.vcf.miss10pc.recode.vcf ${basename}.gq30.ab1e03.miss10pc.vcf
#replace missing IDs '.' with a name containing the chromosome and position.
${bcftoolsdir}./bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ${basename}.gq30.ab1e03.miss10pc.vcf -Ob -o ${basename}.gq30.ab1e03.miss10.rs1.bcf

#remove duplicate (in terms of position) variants here because they seem to give trouble downstream with multiallelic SNPs
${bcftoolsdir}./bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER\n'  ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.bcf > ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.bcf.variants
#(I used R,file: handle_duplicates_in_tili_exomeseq.R)
${bcftoolsdir}./bcftools view -T ^${exomeseq}duplicate_snps_to_remove.txt -Ov -o ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.vcf \
${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.bcf

#now extract the cases and controls, to then get HWE stats separately
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}vcfextract_ctrls.log ${bcftoolsdir}./bcftools view \
-S ${exomedir}list_of_samples_to_keep.ctrls.txt --force-samples -O v -o ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.vcf --threads 4 \
${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.vcf

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}vcfextract_cases.log ${bcftoolsdir}./bcftools view \
-S ${exomedir}list_of_samples_to_keep.cases.txt --force-samples -O v -o ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.vcf --threads 4 \
${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.vcf

#get HWE stats
vcftools --vcf ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.vcf  --hardy --out ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases
vcftools --vcf ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.vcf  --hardy --out ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls

#to remove a subset we need their IDs, and need to make sure they have rsID or otherwise something in that column.
#for now I have used the chr_pos_ref_alt format but I'd rather get proper dbSNP rsIDs.
cat ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.hwe|awk '{if($6<1e-08){print $1"\t"$2}}'> ${exomedir}snps_to_remove_after_hwe_check.ctrls.txt 
cat ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.hwe|awk '{if($6<1e-08){print $1"\t"$2}}'> ${exomedir}snps_to_remove_after_hwe_check.cases.txt 
cat ${exomedir}snps_to_remove_after_hwe_check.ctrls.txt  ${exomedir}snps_to_remove_after_hwe_check.cases.txt > ${exomedir}snps_to_remove_after_hwe_check.txt

#before removing due to hwe deviations, also check for differential missingness, then remove together.

#Also remove sites with differential missingness between cases and controls. to do this we also have to remove duplicates for plink to work.
#here convert to plink and save the case/control information, then merge and give in to plink with -test-missing command 
plink --vcf ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.vcf --make-bed --out ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases
plink --vcf ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.vcf --make-bed --out ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls
cat ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.fam|awk '{print $1,$2,$3,$4,$5,1}' > ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.corr.fam
cat ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.fam|awk '{print $1,$2,$3,$4,$5,2}' > ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.corr.fam
mv ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.corr.fam  ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.fam
mv ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.corr.fam  ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.fam
plink --allow-no-sex --bfile ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls --chr 1-22 --bmerge \
${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.cases --make-bed --out tili.i2.pb.gq30.ab1e03.miss10.rs1
#now we can run the differential missingness test
plink --allow-no-sex --bfile  ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1 --test-missing midp -chr 1-22 --out ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1
cat ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.missing|awk '{if($5<1e-08){print $2}}'> ${exomeseq}snps_to_remove_due_to_diffmiss.txt

#now get chrom and position from the rsIDs, so then bcftools is used to extract the relevant variants 
python from_rsID_to_varpos.py

#now remove
cat ${exomedir}snps_to_remove_due_to_diffmiss.varpos.i2.txt  ${exomedir}snps_to_remove_after_hwe_check.txt > ${exomedir}snps_to_remove_after_hwe_and_diffmiss_check.i2.txt

#the following should work (Petr's suggestion):
${bcftoolsdir}./bcftools view -T ^${exomeseq}snps_to_remove_after_hwe_and_diffmiss_check.i2.txt -Ov -o ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.nodiffmiss.HWE.vcf \
${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.vcf

/software/team152/vcflib/bin/vcfstats ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.nodiffmiss.HWE.vcf > ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.nodiffmiss.HWE.vcfstats &
#list of the remaining variants:
${bcftoolsdir}./bcftools view ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.nodiffmiss.HWE.vcf |grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${exomedir}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.nodiffmiss.HWE.variants

###here we should do annotation and then  filter using the binomial test.
###Only then we should move into sample qc.
python freq_binom_test_v2.py --ann tili.i2.qc.sample_qc3.rsID.vcf.annot --var tili.i2.qc.sample_qc3.rsID.ctrls.variants --out tili_binom_test_results_v3.txt
 
##################################################################################################################################################################################
############################################# END OF FIRST ROUND OF VARIANT QC 						##############################################################################
##################################################################################################################################################################################

##################################################################################################################################################################################
############################################# 			SAMPLE QC1 				##################################################################################################
##################################################################################################################################################################################


## At this stage we have a vcf file with standard site-filtering (VQSR, GQ, missingness, HWE deviations) applied to it. Now we can do some sample filtering too
#rename for brevity:
cp ${exomeseq}tili.i2.pb.gq30.ab1e03.miss10.rs1.nodup.nodiffmiss.HWE.vcf ${exomeseq}tili.i2.qc.vcf
# a) Start by looking into individuals with high missingness  
vcftools --vcf ${exomeseq}tili.i2.qc.vcf --missing-indv --out \
${exomeseq}tili.i2.qc

# b) Heterozygosity 
vcftools --vcf ${exomeseq}tili.i2.qc.vcf --het \
--out ${exomeseq}tili.i2.qc

# b2) heterozygosity at low frequency variants
vcftools --vcf ${exomeseq}tili.i2.qc.vcf --het \
--out ${exomeseq}tili.i2.qc.lowfreq --max-maf 0.01

#c) number of singletons
vcftools --vcf ${exomeseq}tili.i2.qc.vcf  --indv-freq-burden \
--out ${exomeseq}tili.i2.qc

####Following in R will make a list of samples to remove due to being outliers for heterozygosity or having high missingness or high number of singletons    ###################################
inputpathname = "tili.i2.qc"

hetlinesd = 3
hetlinesd2= 4

missingnessline = 0.05
missingnessline2 = 0.1

print(paste("plotting with sd =", hetlinesd, "and miss =s", missingnessline))

imiss = read.table(paste(inputpathname, ".imiss", sep=""), h=T)
imiss$logF_MISS = log10(imiss[,5])


het = read.table(paste(inputpathname, ".het", sep=""), h=T)
het_lowfreq = read.table(paste(inputpathname, ".lowfreq.het", sep=""), h=T)

het$meanHet = (het$N_SITES - het$O.HOM.) / het$N_SITES
het_upper=mean(het$meanHet) + (hetlinesd*sd(het$meanHet))
het_lower=mean(het$meanHet) - (hetlinesd*sd(het$meanHet))

het_lowfreq$meanHet = (het_lowfreq$N_SITES - het_lowfreq$O.HOM.) / het_lowfreq$N_SITES
het_lowfreq_upper=mean(het_lowfreq$meanHet) + (hetlinesd*sd(het_lowfreq$meanHet))
het_lowfreq_lower=mean(het_lowfreq$meanHet) - (hetlinesd*sd(het_lowfreq$meanHet))

n_singl = read.table(paste(inputpathname, ".ifreqburden", sep=""), h=T)
n_singl_upper=mean(n_singl[,3]) + (hetlinesd*sd(n_singl[,3]))
n_singl_lower=mean(n_singl[,3]) - (hetlinesd*sd(n_singl[,3]))

samples_to_remove=unique(c(which(n_singl[,3]>n_singl_upper),which(het_lowfreq$meanHet>het_lowfreq_upper),which(het_lowfreq$meanHet<het_lowfreq_lower),which(het$meanHet>het_upper),which(het$meanHet<het_lower),which(imiss$F_MISS>missingnessline2)))
write.table(het[samples_to_remove,1],"/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/samples_to_remove.i2.qc1.txt",quote=F,col.names=F,row.names=F)
#END    ################################################################################################################################################################################

#note that using 3 sd for number of singletons seems to only remove samples that would have already been removed due to heterozygosity and/or high missingness.

# CAREFUL WITH THE FORMAT OF THE samples_to_remove file, or bcftools will have issues.
${bcftoolsdir}./bcftools view -S ^${exomedir}samples_to_remove.i2.qc1.txt --force-samples -O v -o ${exomedir}tili.i2.qc.sample_qc1.vcf \
--threads 4 ${exomedir}tili.i2.qc.vcf
#check this worked OK
${bcftoolsdir}./bcftools query -l  ${exomedir}tili.i2.qc.vcf | wc -l
${bcftoolsdir}./bcftools query -l  ${exomedir}tili.i2.qc.sample_qc1.vcf | wc -l #number of samples reduced compared to previous file.


# d) Relatedness

#relatedness, measure 1, only using variants with MAF>0.5% (first iteration I had it 1%)
vcftools --vcf ${exomeseq}tili.i2.qc.sample_qc1.vcf --relatedness --maf 0.005 \
--out ${exomeseq}tili.i2.qc.sample_qc1 
#relatedness, measure 2, only using variants with MAF>0.5% (first iteration I had it 1%)
vcftools --vcf ${exomeseq}tili.i2.qc.sample_qc1.vcf --relatedness2 --maf 0.005 \
--out ${exomeseq}tili.i2.qc.sample_qc1

cat ${exomeseq}tili.i2.qc.sample_qc1.relatedness2|awk '{if($7>0.09375 && $1!=$2){print $0}}' #this will give the related individuals

#get pairs of related individuals
cat ${exomeseq}tili.i2.qc.sample_qc1.relatedness2|awk '{if($7>0.09375 && $1!=$2){print $0}}' \
 > ${exomeseq}tili.i2.qc.sample_qc1.highrelatedness.txt
#now pepare file of IDs to remove based on which individual in each pair has more missing data.
python ${exomeseq}prep_high_ibs_samples_to_remove_exomeseq.py --miss  tili.i2.qc.imiss --ibs tili.i2.qc.sample_qc1.highrelatedness.txt --out high_ibs.it2.samples.toremove.txt

${bcftoolsdir}./bcftools view -S ^${exomedir}high_ibs.it2.samples.toremove.txt --force-samples -O v \
-o ${exomedir}tili.i2.qc.sample_qc2.vcf --threads 4 ${exomedir}tili.i2.qc.sample_qc1.vcf

${bcftoolsdir}./bcftools query -l  ${exomedir}tili.i2.qc.sample_qc2.vcf | wc -l #number of samples reduced compared to previous file.

##################################################################################################################################################################################
############################################# 			PCA 				######################################################################################################
##################################################################################################################################################################################
DATADIR="/lustre/scratch115/projects/crohns/exome/TIH/"
plink --vcf ${exomeseq}tili.i2.qc.sample_qc2.vcf --exclude ${DATADIR}aux_and_intermediate_files/high_ld.txt \
--indep-pairwise 1000 50 0.2 --maf 0.01 --chr 1-22 --out ${exomeseq}tili.i2.qc.sample_qc2
plink --vcf ${exomeseq}tili.i2.qc.sample_qc2.vcf --make-bed --out  ${exomeseq}tili.i2.qc.sample_qc2
plink --allow-no-sex --bfile ${exomeseq}tili.i2.qc.sample_qc2 --extract ${exomeseq}tili.i2.qc.sample_qc2.prune.in \
  --make-bed --out ${exomeseq}tili.i2.qc.sample_qc2.pruned

# add case /ctrl status
python annotate_case_ctrl.py --ctrl tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.fam --case tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.fam -f tili.i2.qc.sample_qc2.pruned.fam
grep "\-9\b" tili.i2.qc.sample_qc2.pruned.fam.updated  #to double check it worked
#now move it and continue
mv tili.i2.qc.sample_qc2.pruned.fam.updated tili.i2.qc.sample_qc2.pruned.fam
grep 'rs' tili.i2.qc.sample_qc2.pruned.bim|awk '{print $2}' > ${exomeseq}exomeseq.snps_to_extract_from_1KG.txt
##now calculate PCs using the pruned merged dataset we just created.
sh get_1KG_pc_exomeseq.sh

#merge into one file
plink --bfile ${exomeseq}1KG.chr1.forPCA.exomeseq --merge-list ${exomeseq}1KG.list --make-bed --out ${exomeseq}1KG_for_PCA_allchr.i2
#Remove SNPs which have different alleles between 1KG and exome seq data.
plink --bfile ${exomeseq}tili.i2.qc.sample_qc2.pruned --exclude ${exomeseq}snps_to_remove_before_merging_1KG_exomeseq.txt \
--make-bed --out ${exomeseq}tili.ex.QC2.pruned
plink --bfile ${exomeseq}1KG_for_PCA_allchr.i2 --exclude ${exomeseq}snps_to_remove_before_merging_1KG_exomeseq.txt  --make-bed --out ${exomeseq}1KG_for_PCA_allchr_ready
cat ${exomeseq}1KG_for_PCA_allchr_ready.bim|awk '{print $2}'> ${exomeseq}snps_present_in_both.txt
python check_allele_agreement.py 1KG_for_PCA_allchr_ready.bim tili.ex.QC2.pruned.bim
#the above function also made a list of the SNPs which are in common but differ in the two alleles they have, and which will also have to be removed before merging.
plink --allow-no-sex --bfile ${exomeseq}tili.ex.QC2.pruned --exclude  ${exomeseq}snps_with_different_alleles.txt --make-bed --out ${exomeseq}tili.ex.QC2.pruned.i2
plink --allow-no-sex --bfile ${exomeseq}1KG_for_PCA_allchr_ready --exclude  ${exomeseq}snps_with_different_alleles.txt --make-bed --out ${exomeseq}1KG_for_PCA_allchr_ready.i2

plink --allow-no-sex --bfile ${exomeseq}tili.ex.QC2.pruned.i2 --bmerge ${exomeseq}1KG_for_PCA_allchr_ready.i2 --make-bed \
--extract ${exomeseq}snps_present_in_both.txt  --out ${exomeseq}TILI_1KG_forPCA.i2
#now rerun PCs - **REMEMBER** that smartpca misbehaves (ignores samples) when they have no case/ctrl status assigned, so do this for 1KG prior to running
python annotate_case_ctrl.py --ctrl tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.fam --case tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.fam -f TILI_1KG_forPCA.i2.fam
mv TILI_1KG_forPCA.i2.fam.updated TILI_1KG_forPCA.i2.fam
sed -i 's/-9/1/g' TILI_1KG_forPCA.i2.fam
plink --allow-no-sex --bfile TILI_1KG_forPCA.i2 --recode --out TILI_1KG_forPCA.i2
cat   TILI_1KG_forPCA.i2.ped |awk '{print $1,$2,$3,$4,$5,$6}' > TILI_1KG_forPCA.i2.pedind
smartpca -p TILI_1KG_exomeseq.par > TILI_1KG_exomeseq.Sout

##now visually observe and make removal lists in R using pca_TILI_with1KG_exomeseq.R
cat ${exomeseq}ctrls_to_remove_due_to_PCA_exomeseq.i2.txt ${exomeseq}cases_to_remove_due_to_PCA_exomeseq.i2.txt |awk '{print $1}' > ${exomeseq}all_to_remove_due_to_PCA_exomeseq.i2.txt
${bcftoolsdir}./bcftools view -S ^${exomedir}all_to_remove_due_to_PCA_exomeseq.i2.txt --force-samples -O v -o ${exomedir}tili.i2.qc.sample_qc3.vcf \
--threads 4 ${exomedir}tili.i2.qc.sample_qc2.vcf
${bcftoolsdir}./bcftools query -l tili.i2.qc.sample_qc3.vcf|wc -l
#now get the PCs from this dataset to use as covariates
plink --vcf ${exomedir}tili.i2.qc.sample_qc3.vcf --exclude /lustre/scratch115/projects/crohns/exome/TIH/aux_and_intermediate_files/high_ld.txt \
 --maf 0.01 --chr 1-22 --indep-pairwise 1000 50 0.2 --out ${exomedir}tili.i2.qc.sample_qc3.vcf
 #problem is with the length of SNPs IDs. I will just keep variants with rsIDs for now, (~125K out of ~136K)
grep rs ${exomedir}tili.qc3.prune.in > ${exomedir}tili.qc3.prune.rs.in
plink --allow-no-sex --vcf ${exomedir}tili.i2.qc.sample_qc3.vcf --extract ${exomedir}tili.qc3.prune.rs.in \
--make-bed --out ${exomedir}tili.i2.qc.sample_qc3.pruned
python annotate_case_ctrl.py --ctrl tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.fam --case tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.fam -f tili.i2.qc.sample_qc3.pruned.fam
mv tili.i2.qc.sample_qc3.pruned.fam.updated tili.i2.qc.sample_qc3.pruned.fam
plink --bfile ${exomedir}tili.i2.qc.sample_qc3.pruned --recode --out ${exomedir}tili.i2.qc.sample_qc3.pruned
cat   ${exomedir}tili.i2.qc.sample_qc3.pruned.ped |awk '{print $1,$2,$3,$4,$5,$6}' > ${exomedir}tili.i2.qc.sample_qc3.pruned.pedind
smartpca -p tili.i2.qc3.par > tili.i2.qc3.pca.Sout

##################################################################################################################################################################################
############################################# 			ANNOTATION 				##################################################################################################
##################################################################################################################################################################################

##Now taking ${exomeseq}tili.i2.qc.sample_qc3.vcf downstream, adding rsIDs to it where available then adding annotation.

#Once I have a file in bed format, I can use bcftools annotate to do that

#the following works, but make sure that the -a file is tab-delimited, bgzipped and tabixed
tabix -s1 -b2 -e2 chr${i}.1KG_b37.alleles #see tabixme.sh
#we need to use bcftools annotate across all chromosomes, so i have done so in a small script:
sh ${exomedir}update_rsIDs.sh


${bcftoolsdir}./bcftools view ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.i2.qc.sample_qc3.rsID.variants
${bcftoolsdir}./bcftools view ${exomedir}tili.i2.qc.sample_qc3.rsID.cases.vcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.i2.qc.sample_qc3.rsID.cases.variants
# Once we have the filtered file from above, we can add annotation. 

perl INSTALL.pl --CACHEDIR /software/team152/ensembl-vep/cache/ -s homo_sapiens_merged --NO_TEST #this gives 23/39 failed tests (26/858 failed subtests)
##I have downloaded the cache for GRCh37 in /software/team/ensembl-vep/cache/ To use it, since it's a merged Ensembl & Refseq, I need to add --merged.

############################################################################################################################################################
## for EPACTS installation and trying to make ensembl-vep work, see earlier pipeline_for_tili_exomeseq* files											 ###
############################################################################################################################################################
module add hgi/systems/R/3.1.2

#running vep on the farm (worked fine)
bsub -q normal -J veprun -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep.log \
 /software/team152/ensembl-vep/./vep -i ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf  --offline  \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  ${exomedir}tili.i2.qc.sample_qc3.rsID.vep.stats --use_given_ref \
--output_file ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf.annot --force_overwrite \
--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

# the above produced a separate file, now I will work on incorporating some of these info into the vcf
##note that the fields that are added to the vcf are specified explicitly. 
##note that the *.annot file can be checked for more information
bsub -q normal -J veprun2 -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep2.log \
/software/team152/ensembl-vep/./vep -i ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf --cache --offline \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --use_given_ref --vcf \
--output_file ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf --force_overwrite --per_gene \
--sift b -poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz --symbol \
--fields Consequence,SYMBOL,CANONICAL,SIFT,PolyPhen,CADD_PHRED,EUR_AF,gnomAD_AMR_AF  


##################################################################################################################################################################################
############################################# END OF FIRST ROUND OF VARIANT ANNOTATION 						######################################################################
##################################################################################################################################################################################


########################################
## Single variant test using SNPTEST2 ##
########################################

#SNPTEST2,3 PCs as covariates (Only 2 PCs with p<0.05 in Tracy Widom test)
plink --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf --recode oxford --out ${exomedir}tili.i2.qc.sample_qc3.rsID.annot
#now add covariates and case/ctrl info to the .sample file
python prepare_sample_for_snptest.py -e tili.i2.qc.sample_qc3.pruned.30.evec  -i tili.i2.qc.sample_qc3.pruned.fam -s tili.i2.qc.sample_qc3.rsID.annot.sample
#now run snptest using the updated sample file
snptest -data  ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.gen ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated \
-o ${exomedir}output/tili.i2.qc.sample_qc3.output.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method score 
#run snptest2 with em as well
snptest -data  ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.gen ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated \
-o ${exomedir}output/tili.i2.qc.sample_qc3.output.em.txt -frequentist 1  -cov_names pc1 pc2 -pheno phenotype -method em 


###############################		trying GEMMA for association testing with univariate LMMs		###########################################################################
plink --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf --make-bed --out ${exomedir}tili.i2.qc.sample_qc3.rsID.annot
python ${exomedir}annotate_case_ctrl.py --ctrl tili.i2.pb.gq30.ab1e03.miss10.rs1.ctrls.fam --case tili.i2.pb.gq30.ab1e03.miss10.rs1.cases.fam -f tili.i2.qc.sample_qc3.rsID.annot.fam
mv tili.i2.qc.sample_qc3.rsID.annot.fam.updated tili.i2.qc.sample_qc3.rsID.annot.fam
/software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc.sample_qc3.rsID.annot -gk 2 -o ${exomedir}tili.i2.qc.sample_qc3.rsID.annot #generate standardised relatedness matrix

/software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc.sample_qc3.rsID.annot -k ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.sXX.txt \
-eigen -o tili.i2.qc.sample_qc3.rsID.annot

/software/team152/./gemma.linux -bfile ${exomedir}tili.i2.qc.sample_qc3.rsID.annot \
-d ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.eigenD.txt \
-u ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.eigenU.txt \
-hwe 0.00000001 -maf 0.005 -lmm 4 -o tili.i2.qc.sample_qc3.rsID.annot.gemma_output


###################################################################################################################################
################################################## Gene-based tests using EPACTS ##################################################
###################################################################################################################################

#For file preparation look at 
prepare_group_file_for_EPACTS.sh

#after done with the above step, prepare a ped file 
python prepare_sample_for_epacts.py --sample tili.i2.qc.sample_qc3.rsID.annot.sample.updated
EPACTS_DIR="/software/team152/EPACTS/" 
CHR=21
#vcf file must be tabixed
bgzip tili.i2.qc.sample_qc3.rsID.annot.vcf
tabix tili.i2.qc.sample_qc3.rsID.annot.vcf.gz
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/tili.i2.${CHR}.PD.lenient.0.01.skato.nocov \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/tili.i2.${CHR}.PD.strict.0.01.skato.nocov \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.fc.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/tili.i2.${CHR}.fc.skato.nocov \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done

###		also run with covariates		###
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/tili.i2.${CHR}.PD.lenient.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/tili.i2.${CHR}.PD.strict.0.01.withcov.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 \
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.fc.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/tili.i2.${CHR}.fc.0.01.withcov.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 \
   --pheno DISEASE --test skat --skat-o --run 2
done

##############		Try emmaxVT		##############
${EPACTS_DIR}bin/epacts make-kin  --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
--ped  ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped  --min-maf 0.01 -min-callrate 0.95  \
--out ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.for_epacts.kin --run 2

annot="fc"
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.${annot}.emmaxVT.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.${annot}.genes.epacts \
  --kin ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.for_epacts.kin \
  --out ${exomedir}output/tili.i2.${CHR}.${annot}.0.01.emmaxVT \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01\
   --pheno DISEASE --test emmaxVT --run 2
done

annot="PD.lenient"
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.${annot}.emmaxVT.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.${annot}.genes.epacts \
  --kin ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.for_epacts.kin \
  --out ${exomedir}output/tili.i2.${CHR}.${annot}.0.01.emmaxVT \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01\
   --pheno DISEASE --test emmaxVT --run 2
done 

annot="PD.strict"
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.${annot}.emmaxVT.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.vcf.${CHR}.annot.${annot}.genes.epacts \
  --kin ${exomedir}output/tili.i2.qc.sample_qc3.rsID.annot.for_epacts.kin \
  --out ${exomedir}output/tili.i2.${CHR}.${annot}.0.01.emmaxVT \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01\
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

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}bcftools.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}list_of_samples_to_keep.ctrls.txt --force-samples -O v \
 -o ${exomedir}tili.i2.qc.sample_qc3.rsID.ctrls.vcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf

#b) get cholestatic and mixed cases: (remmeber to remove special cases)
#special cases:
#PRED3825
#PRED5022 - NS
#PRED4004
#PRED4374
#PRED4292
#PRED3715
#PRED5008 - NS

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}bcftoolscholmix.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}cholmix.txt --force-samples -O v \
 -o ${exomedir}tili.i2.qc.sample_qc3.rsID.cholmix.vcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf
 
#c) get hepatocellular and mixed cases (again, remember to remove special cases):

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}bcftoolscholmix.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}hepmix.txt --force-samples -O v \
 -o ${exomedir}tili.i2.qc.sample_qc3.rsID.hepmix.vcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf
 
${bcftoolsdir}./bcftools query -l  ${exomedir}tili.i2.qc.sample_qc3.rsID.hepmix.vcf | wc -l
${bcftoolsdir}./bcftools query -l  ${exomedir}tili.i2.qc.sample_qc3.rsID.cholmix.vcf | wc -l

 #now merge each of the subgroups with ctrls, to prepare the files for the analyses. 
 # Ensure that only sites present in the subgroup and in the controls are kept.
${bcftoolsdir}./bcftools view ${exomedir}tili.i2.qc.sample_qc3.rsID.ctrls.vcf.gz|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.i2.qc.sample_qc3.rsID.ctrls.variants
${bcftoolsdir}./bcftools view ${exomedir}tili.i2.qc.sample_qc3.rsID.cholmix.vcf.gz|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.i2.qc.sample_qc3.rsID.cholmix.variants
${bcftoolsdir}./bcftools view ${exomedir}tili.i2.qc.sample_qc3.rsID.hepmix.vcf.gz|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.i2.qc.sample_qc3.rsID.hepmix.variants
python list_shared_vars.py -1 tili.i2.qc.sample_qc3.rsID.cholmix.variants -2 tili.i2.qc.sample_qc3.rsID.ctrls.variants --out shared_vars_between_cholmix_and_ctrls.txt
python list_shared_vars.py -1 tili.i2.qc.sample_qc3.rsID.hepmix.variants -2 tili.i2.qc.sample_qc3.rsID.ctrls.variants --out shared_vars_between_hepmix_and_ctrls.txt
bgzip  -f ${exomedir}tili.i2.qc.sample_qc3.rsID.ctrls.vcf 
bgzip  -f tili.i2.qc.sample_qc3.rsID.cholmix.vcf
bgzip  -f tili.i2.qc.sample_qc3.rsID.hepmix.vcf
tabix tili.i2.qc.sample_qc3.rsID.hepmix.vcf.gz
tabix tili.i2.qc.sample_qc3.rsID.cholmix.vcf.gz
tabix tili.i2.qc.sample_qc3.rsID.ctrls.vcf.gz
bcftools merge  ${exomedir}tili.i2.qc.sample_qc3.rsID.ctrls.vcf.gz tili.i2.qc.sample_qc3.rsID.cholmix.vcf.gz -O v -o tili.i2.qc3.cholmix_vs_ctrls.vcf --threads 3 --force-samples
bcftools merge  ${exomedir}tili.i2.qc.sample_qc3.rsID.ctrls.vcf.gz tili.i2.qc.sample_qc3.rsID.hepmix.vcf.gz  -O v -o tili.i2.qc3.hepmix_vs_ctrls.vcf  --threads 3 --force-samples

#now do the site filtering
${bcftoolsdir}./bcftools view -T ${exomeseq}shared_vars_between_cholmix_and_ctrls.txt -Ov -o ${exomedir}tili.i2.qc4.cholmix_vs_ctrls.vcf  \
${exomedir}tili.i2.qc3.cholmix_vs_ctrls.vcf
${bcftoolsdir}./bcftools view -T ${exomeseq}shared_vars_between_hepmix_and_ctrls.txt -Ov -o ${exomedir}tili.i2.qc4.hepmix_vs_ctrls.vcf  \
${exomedir}tili.i2.qc3.hepmix_vs_ctrls.vcf

##	now work on the site filtering due to frequency in controls being significantly lower than that in gnomAD NFE ###############################################

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}bcftools.log \
 /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}list_of_samples_to_keep.cases.txt --force-samples -O v \
 -o ${exomedir}tili.i2.qc.sample_qc3.rsID.cases.vcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}tili.i2.qc.sample_qc3.rsID.vcf

python freq_binom_test.py --ann tili.i2.qc.sample_qc3.rsID.vcf.annot --var tili.i2.qc.sample_qc3.rsID.ctrls.variants --out tili_binom_test_results.txt
python freq_binom_test_v2.py --ann tili.i2.qc.sample_qc3.rsID.vcf.annot --var tili.i2.qc.sample_qc3.rsID.ctrls.variants --out tili_binom_test_results_v3.txt

cat tili_binom_test_results_v3.txt|awk '{if( ($7<1e-05 && $6>0 && $6<1) || ($7=="variant_not_in_gnomAD" && ($4/$5)>=0.01 ) ){print $1"\t"$2}}'> vars_to_remove_binom_test_strict.txt

##	##############################################################################################################################################################

${bcftoolsdir}./bcftools view -T ^${exomeseq}vars_to_remove_binom_test_strict.txt -Ov -o ${exomedir}tili.i2.qc3binom_strict.sample_qc3.rsID.annot.vcf \
${exomedir}tili.i2.qc.sample_qc3.rsID.vcf

${bcftoolsdir}./bcftools view -T ^${exomeseq}vars_to_remove_binom_test.txt -Ov -o ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.vcf \
${exomedir}tili.i2.qc3.hepmix_vs_ctrls.vcf
${bcftoolsdir}./bcftools view -T ^${exomeseq}vars_to_remove_binom_test.txt -Ov -o ${exomedir}tili.i2.qc3binom.cholmix_vs_ctrls.vcf \
${exomedir}tili.i2.qc3.cholmix_vs_ctrls.vcf

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

#after done with the above step, prepare a ped file 
python prepare_sample_for_epacts.py --sample tili.i2.qc3binom.hepmix_vs_ctrls.sample.updated
EPACTS_DIR="/software/team152/EPACTS/" 
CHR=21
#vcf file must be tabixed
bgzip tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf
tabix tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.gz

tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.3.annot.fc.genes.epacts
###		also run with covariates		###
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.hepmix_vs_ctrls.${CHR}.PD.lenient.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.hepmix_vs_ctrls.${CHR}.PD.strict.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.fc.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.annot.vcf.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.hepmix_vs_ctrls.${CHR}.fc.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc3binom.hepmix_vs_ctrls.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done


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
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.sample_qc3.rsID.${CHR}.PD.strict.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient1.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.sample_qc3.rsID.${CHR}.PD.lenient.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.fc1.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.annot.vcf.gz \
  --groupf ${exomedir}tili.i2.qc3binom.sample_qc3.rsID.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/tili.i2.qc3binom.sample_qc3.rsID.${CHR}.fc.withcov.0.01.skato \
  --ped ${exomedir}tili.i2.qc.sample_qc3.rsID.annot.sample.updated.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2\
   --pheno DISEASE --test skat --skat-o --run 2
done
