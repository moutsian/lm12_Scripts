# Pipeline for exome sequencing data filtering and analysis for the TILI cohort, using bcftools and vcftools
# Oct-Nov 2017

exomedir="/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/"
bcftoolsdir="/software/team152/bcftools/bcftools1.6/bcftools/"

#Remove variants which don't have PASS for GATK VQSR. Remove monomorphic variants. Keep only biallelic variants. Also keep only the samples in the list provided by Exeter.
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}vcffilterpoly_newbcftools.log /software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -S ${exomedir}list_of_samples_to_keep.both.txt --force-samples -O b -o ${exomedir}tili.poly.biallelic.v5.bcf --threads 4 -c1:minor -m2 -M2 -f PASS ${exomedir}v31_ahmad.vcf.gz

#something was wrong and bftools aborted (glibc exception) -  in version 1.4.1. 
# Make sure you are using the latest bcftools version (1.6 or more)

#to have a look at the variants in the file:
bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER\n'  /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/v31_ahmad.vcf.gz > /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/v31_ahmad.vcf.variants
${bcftoolsdir}./bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER\n'  ${exomedir}tili.poly.biallelic.v5.bcf > ${exomedir}tili.poly.biallelic.v5.bcf.variants

##Set variants with GQ<30 to missing

#************* !! NOTE THAT THE GQ genotype filter  DOES NOT WORK PROPERLY when outputting it as bcf file (at least not with vcftools 0.1.14)
# Just output it as vcf.
${bcftoolsdir}bcftools view tili.poly.biallelic.v5.bcf -Ov -o tili.poly.biallelic.v5.vcf
vcftools --vcf tili.poly.biallelic.v5.vcf  --minGQ 30  --recode --recode-INFO-all --out tili.poly.biallelic.v5.gq30
# trying also with vcffilter:
vcffilter -f "GQ > 30" tili.poly.biallelic.v5.vcf > tili.poly.biallelic.v5.gq30.vcffilter.vcf
##vcffilter works. You want to keep the spaces - formatting is a bit funny. The problem is that it doesn't work with the AD format.. 

#filter variants for missingness
vcftools --vcf tili.poly.biallelic.v5.gq30.vcffilter.vcf --max-missing 0.9 --recode --recode-INFO-all --out tili.poly.biallelic.v5.gq30.miss10pc.vcffilter
vcftools --vcf tili.poly.biallelic.v5.gq30.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out tili.poly.biallelic.v5.gq30.miss10pc

#After filtering, we kept 407613 out of a possible 559130 Sites  -> This is due to the GQ>=30 threshold. Without this, less than 35-40K variants would be filtered out.

/software/team152/vcflib/bin/vcfstats ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.recode.vcf

## Now consider looking for allelic imbalance or filter on the basis of alternate allele read depth. Then move on to remove sites with missingness below a certain threshold.  
# Note that I am not filtering for allelic imbalance at present because I don't have the relevant field from gatk and it's not feasible atm with bcftools, vcftools or vcffilter tools.
# Also not filtering for alternate allele read depth for similar reasons

#get stats of site missingness to confirm
vcftools --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.recode.vcf  --missing-site --out ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.recode

#replace missing IDs '.' with a name containing the chromosome and position.
${bcftoolsdir}./bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.recode.vcf -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.rs1.bcf
#now extract the controls, to then get HWE stats for controls

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}vcfextract_ctrls.log ${bcftoolsdir}./bcftools view \
-S ${exomedir}list_of_samples_to_keep.ctrls.txt --force-samples -O v -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls.vcf --threads 4  -m2 -M2 -f PASS \
${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.rs1.bcf

#Here's how to test for allelic imbalance from P. Danacek (given to me by Hilary, 26/10, haven't run it yet ---update Nov 6th: ran it further down, have a look) 
#nfs/users/nfs_h/hcm/bcftools/bcftools +setGT --threads 0 -O b $inputbcf --  -n . -t 'b:AD<'0.001

#get HWE stats
vcftools --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls.vcf  --hardy --out ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls

#to remove a subset we need their IDs, and need to make sure they have rsID or otherwise something in that column.
#for now I have used the chr_pos_ref_alt format but I'd rather get proper dbSNP rsIDs.
${bcftoolsdir}./bcftools view ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls.vcf |awk '{print $1,$2,$3,$4,$5,$6,$7,$8}'|grep -v '#'> tili.poly.biallelic.v5.gq30.miss10pc.ctrls.variants
python create_snplist_to_extract.py #this outputs a file with rsIDs, which should be fine for vcftools, though it does not work in practice
cat ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls.hwe|awk '{if($6>1e-08){print $1"\t"$2}}'> ${exomedir}snps_to_keep_after_hwe_check.txt #this is for bcftools -R option, though also doesn't work atm
# 406134 snps_to_keep_after_hwe_check.txt

#bcftools complains that the bcf file wasn't zipped with bgzip. To overcome this, I use the -Ob and index commands of bcftools as follows:
${bcftoolsdir}./bcftools view  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls.vcf -I -Ob -o  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls.bcf  #i) bgzip
${bcftoolsdir}./bcftools index ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.ctrls.bcf #ii) indexing

#I am also removing duplicate entries because they seem to give trouble downstream with multiallelic SNPs 
(I used R,file: handle_duplicates_in_tili_exomeseq.R)

#/software/team152/bcftools/bin/./bcftools concat --rm-dup all -a -o antegamisou_rmvdup.bcf  tili.poly.biallelic.v5.gq30.miss10pc.ctrls.bcf
#/software/team152/bcftools/bin/./bcftools view -R ${exomedir}snps_to_keep_after_hwe_check_nodup.txt --o ${exomedir}tili.poly.biallelic.qc1.bcf ${exomedir}antegamisou.bcf 
#the following should work (Petr's suggestion):
${bcftoolsdir}./bcftools view -T ${exomedir}snps_to_keep_after_hwe_check_nodup_updated.txt -c1:minor -m2 -M2  -Ov -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.rs1.bcf

/software/team152/vcflib/bin/vcfstats ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcf > ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcfstats &

#list of the remaining variants:
${bcftoolsdir}./bcftools view ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcf |grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.variants

##################################################################################################################################################################################
############################################# END OF FIRST ROUND OF VARIANT QC 						##############################################################################
##################################################################################################################################################################################

##################################################################################################################################################################################
############################################# 			SAMPLE QC1 				##################################################################################################
##################################################################################################################################################################################


## At this stage we have a vcf file with standard site-filtering (VQSR, GQ, missingness, HWE deviations) applied to it. Now we can do some sample filtering too

# a) Start by looking into individuals with high missingness  
vcftools --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcf --missing-indv --out \
${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass

# b) Heterozygosity 
vcftools --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcf --het \
--out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass

# b2) heterozygosity at low frequency variants
vcftools --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcf --het \
--out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.lowfreq --max-maf 0.01

####Following in R will make a list of samples to remove due to being outliers for heterozygosity or having high missingness or high number of singletons    ###################################
inputpathname = "tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.rsID"

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
write.table(het[samples_to_remove,c(1:2)],"/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/samples_to_remove.qc1.txt",quote=F,col.names=T,row.names=F)
#END    ################################################################################################################################################################################

#note that using 3 sd for number of singletons seems to only remove samples that would have already been removed due to heterozygosity and/or high missingness.

${bcftoolsdir}./bcftools view -S ^${exomedir}samples_to_remove.qc1.txt --force-samples -O v -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc1.vcf --threads 4 ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.vcf
#check this worked OK
${bcftoolsdir}./bcftools query -l  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc1.vcf | wc -l #number of samples reduced compared to previous file.


# d) Relatedness

#relatedness, measure 1, only using variants with MAF>1%
vcftools --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc1.vcf --relatedness --maf 0.01 \
--out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc1 

#relatedness, measure 2, only using variants with MAF>1%
vcftools --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc1.vcf --relatedness2 --maf 0.01 \
--out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc1 
cat tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.relatedness2|awk '{if($7>0.09375 && $1!=$2){print $0}}' #this will give the related individuals

#get pairs of related individuals
cat ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.relatedness2|awk '{if($7>0.09375 && $1!=$2){print $0}}' \
 > ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.highrelatedness.txt
#now prepare file of IDs to remove based on which individual in each pair has more missing data.
python ${exomeseq}prep_high_ibs_samples_to_remove_exomeseq.py
${bcftoolsdir}./bcftools view -S ^${exomedir}high_ibs.qc1.samples.toremove.txt --force-samples -O v -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.vcf --threads 4 ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc1.vcf
${bcftoolsdir}./bcftools query -l  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.vcf | wc -l #number of samples reduced compared to previous file.

# FOR THE VERSION WHERE WE WERE DOING THE CHECKS WITH GENO DATA HAVE A LOOK AT pipeline_for_tili_exomeseq_clean26Oct.sh

#Remove sites with differential missingness between cases and controls. ***For next run move this to the site QC part.
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomeseq}ctrls.log ${bcftoolsdir}./bcftools view \
 -S ${exomeseq}list_of_samples_to_keep.ctrls.txt --force-samples -O v -o ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.vcf \
 --threads 4   ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.vcf
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomeseq}cases.log ${bcftoolsdir}./bcftools view \
 -S ${exomeseq}list_of_samples_to_keep.cases.txt --force-samples -O v -o ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.vcf \
 --threads 4   ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.vcf

#here convert to plink and save the case/control information, then merge and give in to plink with -test-missing command 
plink --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.vcf --make-bed --out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases
plink --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.vcf --make-bed --out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls
cat ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.fam|awk '{print $1,$2,$3,$4,$5,1}' > ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.corr.fam
cat ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.fam|awk '{print $1,$2,$3,$4,$5,2}' > ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.corr.fam
mv ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.corr.fam  ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.fam
mv ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.corr.fam  ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.fam
plink --allow-no-sex --bfile ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls --bmerge \
${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases --make-bed --out tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2
#now we can run the differential missingness test
plink --allow-no-sex --bfile  ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2 --test-missing midp -chr 1-22 --out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2
cat ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.missing|awk '{if($5<1e-08){print $2}}'> ${exomeseq}snps_to_remove_due_to_diffmiss.txt

#now get chrom and position from the rsIDs, so then bcftools is used to extract the relevant variants 
python from_rsID_to_varpos.py
${bcftoolsdir}./bcftools view -T ^${exomeseq}snps_to_remove_due_to_diffmiss.varpos.txt -Ov -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.vcf \
${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.vcf


##################################################################################################################################################################################
############################################# 			PCA 				######################################################################################################
##################################################################################################################################################################################

DATADIR="/lustre/scratch115/projects/crohns/exome/TIH/"
plink --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.vcf --exclude ${DATADIR}aux_and_intermediate_files/high_ld.txt \
--indep-pairwise 1000 50 0.2 --maf 0.01 --chr 1-22 --out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2
plink --vcf ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.vcf --make-bed --out  ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2
plink --allow-no-sex --bfile ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2 --extract ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.prune.in \
  --make-bed --out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.pruned

# add case /ctrl status
python annotate_case_ctrl.py
grep "\-9\b" tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.pruned.fam.updated  #to double check it worked
#now move it and continue
mv tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.pruned.fam.updated tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.pruned.fam
grep 'rs' tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.pruned.bim|awk '{print $2}' > ${exomeseq}exomeseq.snps_to_extract_from_1KG.txt
##now calculate PCs using the pruned merged dataset we just created.
sh get_1KG_pc_exomeseq.sh


#merge into one file
plink --bfile ${exomeseq}1KG.chr1.forPCA.exomeseq --merge-list ${exomeseq}1KG.list --make-bed --out ${exomeseq}1KG_for_PCA_allchr
#Remove SNPs which have different alleles between 1KG and exome seq data.
plink --bfile ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.pruned --exclude ${exomeseq}snps_to_remove_before_merging_1KG_exomeseq.txt \
--make-bed --out ${exomeseq}tili.ex.QC2.pruned
plink --bfile ${exomeseq}1KG_for_PCA_allchr --exclude ${exomeseq}snps_to_remove_before_merging_1KG_exomeseq.txt  --make-bed --out ${exomeseq}1KG_for_PCA_allchr_ready
cat ${exomeseq}1KG_for_PCA_allchr_ready.bim|awk '{print $2}'> ${exomeseq}snps_present_in_both.txt
python check_allele_agreement.py
python check_allele_agreement.py 1KG_for_PCA_allchr_ready.bim tili.ex.QC2.pruned.bim
plink --allow-no-sex --bfile ${exomeseq}tili.ex.QC2.pruned --bmerge ${exomeseq}1KG_for_PCA_allchr_ready --make-bed \
--extract ${exomeseq}snps_present_in_both.txt  --out ${exomeseq}TILI_1KG_forPCA
#plink --allow-no-sex --bfile TILI_1KG_forPCAtmp --extract snps_to_keep_for_joint_PCA.txt --out TILI_1KG_forPCA --make-bed 
#now rerun PCs - **REMEMBER** that smartpca misbehaves (ignores samples) when they have no case/ctrl status assigned, so do this for 1KG prior to running
python annotate_case_ctrl.py --ctrl tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.fam --case tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.fam -f TILI_1KG_forPCA.fam
mv TILI_1KG_forPCA.fam.updated TILI_1KG_forPCA.famsed -i 's/-9/1/g' TILI_1KG_forPCA.fam
plink --allow-no-sex --bfile TILI_1KG_forPCA --recode --out TILI_1KG_forPCA
cat   TILI_1KG_forPCA.ped |awk '{print $1,$2,$3,$4,$5,$6}' > TILI_1KG_forPCA.pedind
#sed -i 's/urn:wtsi//g' TILI_1KG_forPCA.ped*
smartpca -p TILI_1KG_exomeseq.par > TILI_1KG_exomeseq.Sout

#file to use: ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.vcf
cat ${exomeseq}ctrls_to_remove_due_to_PCA_exomeseq.txt ${exomeseq}cases_to_remove_due_to_PCA_exomeseq.txt |awk '{print $1}' > ${exomeseq}all_to_remove_due_to_PCA_exomeseq.txt
${bcftoolsdir}./bcftools view -S ^${exomedir}all_to_remove_due_to_PCA_exomeseq.txt --force-samples -O v -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vcf \
--threads 4 ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.vcf
#shorten filenames so smartpca doesn't complain
cp ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vcf ${exomedir}tili.qc3.vcf
#now get the PCs from this dataset to use as covariates
plink --vcf ${exomedir}tili.qc3.vcf --exclude /lustre/scratch115/projects/crohns/exome/TIH/aux_and_intermediate_files/high_ld.txt \
 --maf 0.01 --chr 1-22 --indep-pairwise 1000 50 0.2 --out ${exomedir}tili.qc3
 #problem is with the length of SNPs IDs. I will just keep variants with rsIDs for now, (~125K out of ~136K)
grep rs ${exomedir}tili.qc3.prune.in > ${exomedir}tili.qc3.prune.rs.in
plink --allow-no-sex --vcf ${exomedir}tili.qc3.vcf --extract ${exomedir}tili.qc3.prune.rs.in \
--make-bed --out ${exomedir}tili.qc3.pruned
python annotate_case_ctrl.py --ctrl tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.fam --case tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.fam -f tili.qc3.pruned.fam
mv tili.qc3.pruned.fam.updated tili.qc3.pruned.fam
plink --bfile ${exomedir}tili.qc3.pruned --recode --out ${exomedir}tili.qc3.pruned
cat   ${exomedir}tili.qc3.pruned.ped |awk '{print $1,$2,$3,$4,$5,$6}' > ${exomedir}tili.qc3.pruned.pedind
smartpca -p tili.qc3.par > tili.qc3.pruned.pca.Sout

##################################################################################################################################################################################
############################################# 			ANNOTATION 				##################################################################################################
##################################################################################################################################################################################

##Now taking ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.nodiffmiss.qc2.vcf downstream, adding rsIDs to it where available then adding annotation.

#Once I have a file in bed format, I can use bcftools annotate to do that

#the following works, but make sure that the -a file is tab-delimited, bgzipped and tabixed
tabix -s1 -b2 -e2 chr${i}.1KG_b37.alleles #see tabixme.sh
#bcftools annotate -c CHROM,POS,ID,REF,ALT -a /lustre/scratch115/projects/ibdgwas/aux_files/chr${i}.1KG_b37.alleles.gz   -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.bcf   /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.bcf
#we need to use bcftools annotate across all chromosomes, so i have done so in a small script:
sh ${exomedir}update_rsIDs.sh

${bcftoolsdir}./bcftools view ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.rsID.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.rsID.variants
# Once we have the filtered file from above, we can add annotation. 

perl INSTALL.pl --CACHEDIR /software/team152/ensembl-vep/cache/ -s homo_sapiens_merged --NO_TEST #this gives 23/39 failed tests (26/858 failed subtests)
##I have downloaded the cache for GRCh37 in /software/team/ensembl-vep/cache/ To use it, since it's a merged Ensembl & Refseq, I need to add --merged.


## EPACTS installation ###
module add hgi/systems/R/3.1.2
git clone https://github.com/statgen/EPACTS.git
cd EPACTS
./configure CC=/software/gcc-4.9.2/bin/gcc CXX=/software/gcc-4.9.2/bin/g++ prefix=/software/team152/EPACTS/
aclocal
automake
make
make install
## 09/10/17: The above worked fine. I will add R.3.1.2 in the system so the module loads automatically.

#trying to make ensembl-vep work######################################################################################################
## -example below works
#./vep -i examples/homo_sapiens_GRCh37.vcf --cache --force_overwrite --merged --dir /software/team152/ensembl-vep/cache/ --port 3337

##this still creates a separate file, this time with a ton of information (Sift, Polyphen, CADD, AFs). Again it worked.
#/software/team152/ensembl-vep/./vep -i examples/homo_sapiens_GRCh37.vcf --cache --force_overwrite --merged --dir /software/team152/ensembl-vep/cache/ --offline --canonical  \
#  --stats_file  /software/team152/ensembl-vep/cache/vep.stats --output_file /software/team152/ensembl-vep/cache/homo_sapiens_GRCh37.vcf.annot  --sift p --poly p --af_1kg --af_gnomad\
# --regulatory --use_given_ref --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz
########################################################################################################################################

# now for my actual file - make a full version as above, and also see if we can make an updated vcf automatically with a small number of fields.
#tab-delimited
##first convert bcf file to vcf, for input to vep
#${bcftoolsdir}./bcftools view -O v -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.rsID.vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.rsID.bcf

##now run vep 
#/software/team152/ensembl-vep/./vep -i /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf  --offline  \
#--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  /software/team152/ensembl-vep/cache/vep.stats --use_given_ref \
#--output_file /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf.annot --force_overwrite \
#--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

#running vep on the farm (worked fine)
bsub -q normal -J veprun -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep.log \
 /software/team152/ensembl-vep/./vep -i ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vcf  --offline  \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  /${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vep.stats --use_given_ref \
--output_file ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vcf.annot --force_overwrite \
--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

# the above produced a separate file, now I will work on incorporating some of these info into the vcf
##note that the fields that are added to the vcf are specified explicitly. 
##note that the *.annot file can be checked for more information
bsub -q normal -J veprun2 -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}vep2.log \
/software/team152/ensembl-vep/./vep -i ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vcf --cache --offline \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --use_given_ref --vcf \
--output_file ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf --force_overwrite --per_gene \
--sift b -poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz --symbol \
--fields Consequence,SYMBOL,CANONICAL,SIFT,PolyPhen,CADD_PHRED,EUR_AF,gnomAD_AMR_AF  


##################################################################################################################################################################################
############################################# END OF FIRST ROUND OF VARIANT ANNOTATION 						######################################################################
##################################################################################################################################################################################


########################################
## Single variant test using SNPTEST2 ##
########################################

#SNPTEST2,3 PCs as covariates (3 PCs with p<0.05 in Tracy Widom test)
plink --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf --recode oxford --out ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot
#now add covariates and case/ctrl info to the .sample file
python prepare_sample_for_snptest.py -e TILI.qc3.30.evec  -i tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.fam -s tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample
#now run snptest using the updated sample file
snptest -data  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.gen ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.updated \
-o ${exomedir}output/tili.qc3_output.txt -frequentist 1  -cov_names pc1 pc2 pc3 -pheno phenotype -method score 

##too many significant hits. Check differential missingness again
plink --allow-no-sex -data ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot --test-missing midp -chr 1-22 --out ${exomeseq}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot
plink --allow-no-sex  --bfile tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases --missing --out tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases
plink --allow-no-sex  --bfile tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls --missing --out tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls
grep -v '#' tili.qc3_output.txt|awk '{if($42<5e-08 && $30>0.005 && $31>0.005){print $1,$2,$4,$30,$31"\t"$42,$44}}' |more #this shows some gw significant hits for which the risk alleles are present in both cases and controls with MAF>0.5% in both
# For downstream analysis I will also remove variants with MAF < 0.5% in either cases or controls. I need to remove variants which are far out of HWE in cases too.

#run snptest2 with em as well ,to double check 
snptest -data  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.gen ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.updated \
-o ${exomedir}output/tili.qc3_output.em.txt -frequentist 1  -cov_names pc1 pc2 pc3 -pheno phenotype -method em

###############################		trying GEMMA for association testing with univariate LMMs		###########################################################################
plink --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vcf --make-bed --out ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3
python ${exomedir}annotate_case_ctrl.py --ctrl tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.ctrls.fam --case tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.cases.fam -f tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.fam
mv tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.fam.updated tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.fam
/software/team152/./gemma.linux -bfile ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3 -gk 2 -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3 #generate standardised relatedness matrix

/software/team152/./gemma.linux -bfile ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3 -k ${exomedir}output/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.sXX.txt \
-eigen -o tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3

/software/team152/./gemma.linux -bfile ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3 \
-d ${exomedir}output/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.eigenD.txt \
-u ${exomedir}output/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.eigenU.txt \
-hwe 0.00000001 -maf 0.005 -lmm 4 -o tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.gemma_output


###################################################################################################################################
################################################## Gene-based tests using EPACTS ##################################################
###################################################################################################################################

#For file preparation look at prepare_group_file_for_EPACTS.sh

#after done with the above step, prepare a ped file 
python prepare_sample_for_epacts.py --sample tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample
EPACTS_DIR="/software/team152/EPACTS/" 
CHR=21
#vcf file must be tabixed
bgzip tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf
tabix tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/test.gene.${CHR}.PD.lenient.0.01.skato \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/test.gene.${CHR}.PD.strict.0.01.skato \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01 \
   --pheno DISEASE --test skat --skat-o --run 2
done
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.fc.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/test.gene.${CHR}.fc.skato \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.05 \
   --pheno DISEASE --test skat --skat-o --run 2
done

###		also run with covariates		###
for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.lenient.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.PD.lenient.genes.epacts --out ${exomedir}output/test.gene.${CHR}.PD.lenient.withcov.0.01.skato \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 --cov pc3 \
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.PD.strict.genes.epacts --out ${exomedir}output/test.gene.${CHR}.PD.strict.0.01.withcov.skato \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 --cov pc3\
   --pheno DISEASE --test skat --skat-o --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.strict.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.fc.genes.epacts --out ${exomedir}output/test.gene.${CHR}.fc.0.01.withcov.skato \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01 --cov pc1 --cov pc2 --cov pc3\
   --pheno DISEASE --test skat --skat-o --run 2
done
##############		Try emmaxVT		##############
${EPACTS_DIR}bin/epacts make-kin  --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
--ped  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped  --min-maf 0.01 -min-callrate 0.95  \
--out ${exomedir}output/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.kin --run 2

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.emmaxVT.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.PD.strict.genes.epacts \
  --kin ${exomedir}/output/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.kin \
  --out ${exomedir}output/test.gene.${CHR}.PD.strict.0.01.emmaxVT \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01\
   --pheno DISEASE --test emmaxVT --run 2
done

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.emmaxVT.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.PD.lenient.genes.epacts \
  --kin ${exomedir}/output/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.kin \
  --out ${exomedir}output/test.gene.${CHR}.PD.lenient.0.01.emmaxVT \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01\
   --pheno DISEASE --test emmaxVT --run 2
done 

for((CHR=1;CHR<=22;CHR++));do
bsub -q normal -J epactsrun -G team152 -n2 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o ${exomedir}epacts${CHR}.PD.emmaxVT.log \
${EPACTS_DIR}bin/epacts group --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.gz \
  --groupf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf.${CHR}.annot.fc.genes.epacts \
  --kin ${exomedir}/output/tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.kin \
  --out ${exomedir}output/test.gene.${CHR}.fc.0.01.emmaxVT \
  --ped ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.for_epacts.ped --max-maf 0.01\
   --pheno DISEASE --test emmaxVT --run 2
done 

 
 ##too many hits...
cat tili.qc3_output.txt|awk '{if($1==9 && $4>86452019 && $4<86530228 && $42!="NA"){print $1,$2,$4,$30,$31,$42,$44}}' 

###################################################################################################################################
################################ Variant QC, part 2 ###############################################################################
###################################################################################################################################
exomedir="/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/"
bcftoolsdir="/software/team152/bcftools/bcftools1.6/bcftools/"

# Filter by allelic imbalance - then missingness again, then case control differential missingness check, then deviations from HWE in *cases*
# Try the filtering by allele-balance that Hilary uses (note that Hilary is using a cutoff of p=1e-03)
${bcftoolsdir}./bcftools +setGT --threads 0 -Ov -o ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc4.ab.1e03.vcf  \
${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.vcf --  -n . -t 'b:AD<1e-03'

# missingness filter
vcftools --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc4.ab1e03.vcf --max-missing 0.9 --recode --recode-INFO-all --out ${exomedir}tili.qc4.ab1e03.miss10pc.vcf
mv ${exomedir}tili.qc4.ab1e03.miss10pc.vcf.recode.vcf ${exomedir}tili.qc4.ab1e03.miss10pc.vcf
#kept 393709 out of a possible 394808 sites for p=1e-05, 393242 for p=1e-03.
# hwe in cases
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomedir}cases.log ${bcftoolsdir}./bcftools view \
 -S ${exomedir}list_of_samples_to_keep.cases.txt --force-samples -O v -o ${exomedir}tili.qc4.ab1e03.miss10pc.cases.vcf \
 --threads 4   ${exomedir}tili.qc4.ab1e03.miss10pc.vcf

vcftools --vcf ${exomedir}tili.qc4.ab1e03.miss10pc.cases.vcf  --hardy --out ${exomedir}tili.qc4.ab1e03.miss10pc.cases
${bcftoolsdir}./bcftools view ${exomedir}tili.qc4.ab1e03.miss10pc.cases.vcf |awk '{print $1,$2,$3,$4,$5,$6,$7,$8}'|grep -v '#'> tili.qc4.ab1e03.miss10pc.cases.variants
cat ${exomedir}tili.qc4.ab1e03.miss10pc.cases.hwe|awk '{if($6>1e-08){print $1"\t"$2}}' > ${exomedir}snps_to_keep_after_hwe_check_cases.2.txt

#Interestingly, only 4 SNPs seem to now be out of HWE in cases - this could be an effect of allelic imbalance filter capturing the "errors"
#remember to remove header before running bcftools
${bcftoolsdir}./bcftools view -T ${exomedir}snps_to_keep_after_hwe_check_cases.2.txt -Ov -o ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.vcf ${exomedir}tili.qc4.ab1e03.miss10pc.vcf

#Remove sites with differential missingness between cases and controls. ***For next run move this to the site QC part.
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomeseq}ctrls.log ${bcftoolsdir}./bcftools view \
 -S ${exomeseq}list_of_samples_to_keep.ctrls.txt --force-samples -O v -o ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls.vcf \
 --threads 4   ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.vcf
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o ${exomeseq}cases.log ${bcftoolsdir}./bcftools view \
 -S ${exomeseq}list_of_samples_to_keep.cases.txt --force-samples -O v -o ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases.vcf \
 --threads 4   ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.vcf


#here convert to plink and save the case/control information, then merge and give in to plink with -test-missing command 
plink --vcf ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases.vcf --make-bed --out ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases
plink --vcf ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls.vcf --make-bed --out ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls
cat ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases.fam|awk '{print $1,$2,$3,$4,$5,2}' > ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases.corr.fam
cat ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls.fam|awk '{print $1,$2,$3,$4,$5,1}' > ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls.corr.fam
mv ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls.corr.fam  ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls.fam
mv ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases.corr.fam  ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases.fam
plink --allow-no-sex --bfile ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.ctrls --bmerge \
${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.cases --make-bed --out ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE

#now we can run the differential missingness test
plink --allow-no-sex --bfile  ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE --test-missing midp -chr 1-22 --out ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE
cat ${exomeseq}tili.qc4.ab1e03.miss10pc.HWE.missing|awk '{if($5<1e-08){print $2}}'> ${exomeseq}snps_to_remove_due_to_diffmiss.3.txt 
#now get chrom and position from the rsIDs, so then bcftools is used to extract the relevant variants 
${bcftoolsdir}./bcftools view ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.vcf |awk '{print $1,$2,$3,$4,$5,$6,$7,$8}'|grep -v '#'> tili.qc4.ab1e03.miss10pc.HWE.variants
python from_rsID_to_varpos.py
${bcftoolsdir}./bcftools view -T ^${exomeseq}snps_to_remove_due_to_diffmiss.varpos.3.txt -Ov -o ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.vcf \
${exomedir}tili.qc4.ab1e03.miss10pc.HWE.vcf
 
 
###############################		trying GEMMA for association testing with univariate LMMs		###########################################################################
plink --vcf ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.vcf --make-bed --out ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss

python ${exomedir}annotate_case_ctrl.py --ctrl tili.qc4.ab1e03.miss10pc.HWE.ctrls.fam \
--case tili.qc4.ab1e03.miss10pc.HWE.cases.fam -f tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.fam

mv tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.fam.updated tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.fam
/software/team152/./gemma.linux -bfile ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss -maf 0.005 -gk 2 -o tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss #generate standardised relatedness matrix
##note that I have NOT excluded the MHC from the relatedness matrix.
/software/team152/./gemma.linux -bfile ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss -k output/tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.sXX.txt \
-eigen -o tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss

/software/team152/./gemma.linux -bfile ${exomedir}tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss \
-d ${exomedir}output/tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.eigenD.txt \
-u ${exomedir}output/tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.eigenU.txt \
-hwe 0.00000001 -maf 0.005 -lmm 4 -o tili.qc4.ab1e03.miss10pc.HWE.nodiffmiss.gemma_output


#Re-check SNPTEST2,3 PCs as covariates (3 PCs with p<0.05 in Tracy Widom test, which I haven't regenerated after latest variant QC)
plink --vcf ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.vcf --recode oxford --out ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot
#now add covariates and case/ctrl info to the .sample file
python prepare_sample_for_snptest.py -e TILI.qc3.30.evec  -i tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc2.fam -s tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample
#now run snptest using the updated sample file
snptest -data  ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.gen ${exomedir}tili.poly.biallelic.v5.gq30.miss10pc.HWEpass.qc3.annot.sample.updated \
-o ${exomedir}output/tili.qc3_output.txt -frequentist 1  -cov_names pc1 pc2 pc3 -pheno phenotype -method score 



