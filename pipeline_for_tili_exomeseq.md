# Pipeline for exome sequencing data filtering and analysis for the TILI cohort, using bcftools and vcftools

#Remove variants which don't have PASS. Keep only biallelic variants. Also keep only the samples in the list.
#At the moment we still need the IDs of 46 samples from Gareth/the Exeter guys.

bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/vcffilter1poly.log /software/team152/bcftools/bin/./bcftools view -S /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/list_of_samples_to_keep.cases.txt -O b -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili_cases.biallelic.bcf --threads 4 -m2 -M2 -f PASS /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/v31_ahmad.vcf.gz
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/vcffilter2poly.log /software/team152/bcftools/bin/./bcftools view -S /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/list_of_samples_to_keep.ctrls.txt --force-samples -O b -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili_ctrls.biallelic.bcf --threads 4  -m2 -M2 -f PASS /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/v31_ahmad.vcf.gz
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/vcffilter3poly.log /software/team152/bcftools/bin/./bcftools view -S /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/list_of_samples_to_keep.both.txt --force-samples -O b -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.bcf --threads 4 -c1:minor -m2 -M2 -f PASS /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/v31_ahmad.vcf.gz
/software/team152/bcftools/bin/./bcftools stats tili_cases.biallelic.bcf
#something is wrong and bftools aborts (glibc exception)
#to have a look at the variants in the file:
bcftools query -f'%CHROM %POS %ID %REF %ALT\n'  /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili_ctrls.biallelic.bcf > /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili_ctrls.biallelic.variants
bcftools query -f'%CHROM %POS %ID %REF %ALT\n'  /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili_ctrls.biallelic.poly.bcf > /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili_ctrls.biallelic.poly.variants
bcftools query -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER\n'  /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/v31_ahmad.vcf.gz > /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/v31_ahmad.vcf.variants
bcftools query -f %CHROM %POS %ID %REF %ALT %QUAL %FILTER %INFO\n  maf.0.01.cr90_gq20_dp8.postqc_analysis_interval_oxford_eu.all.chrs.norelout1.oxf_int.baits.uk10k.exac3.annots.vep75.correct.vqsr999.both.vcf.gz


/software/team152/bcftools/bin/./bcftools stats -S tili.poly.biallelic.bcf

test.qc1e.ped
#filter variants for missingness
vcftools --bcf tili.poly.biallelic.bcf --max-missing 0.95 --recode-bcf --recode-INFO-all --out tili.poly.biallelic.nomiss
# After filtering, kept 523313 out of a possible 559130 Sites
#get stats of site missingness to confirm
vcftools --bcf tili.poly.biallelic.nomiss.recode.bcf  --missing-site --out tili.poly.biallelic.nomiss.recode
#replace missing IDs '.' with a name containing the chromosome and position.
/software/team152/bcftools/bin/./bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' tili.poly.biallelic.nomiss.recode.bcf -o tili.poly.biallelic.nomiss.recode.rs1.bcf
#now extract the controls, to then get HWE stats for controls
bsub -q normal -J vcf_filter -G team152 -n4 -R "span[hosts=1] select[mem>2500] rusage[mem=2500]" -M2500 -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/vcfextract_ctrls.log /software/team152/bcftools/bin/./bcftools view -S /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/list_of_samples_to_keep.ctrls.txt --force-samples -O b -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.ctrls.bcf --threads 4  -m2 -M2 -f PASS /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.bcf

#get HWE stats
vcftools --bcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.ctrls.bcf  --hardy --out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.ctrls.bcf
#to remove a subset we need their IDs, and need to make sure they have rsID or otherwise something in that column.
#for now I have used the chr_pos_ref_alt format but I'd rather get proper dbSNP rsIDs.
/software/team152/bcftools/bin/./bcftools view /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.bcf|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}'|grep -v '#'> tili.poly.biallelic.nomiss.recode.rs1.variants
python create_snplist_to_extract.py #this outputs a file with rsIDs, which should be fine for vcftools, though it does not work in practice
cat tili.poly.biallelic.nomiss.recode.rs1.ctrls.bcf.hwe|awk '{if($6>1e-08){print $1"\t"$2}}'> snps_to_keep_after_hwe_check.txt #this is for bcftools -R option, though also doesn't work atm

#option 1: Using vcftools (doesn't work atm)
vcftools --bcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.bcf \
--exclude /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.outofHWE.txt --recode-bcf \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.bcf
#option 2: Using bcftools (doesn't work as it is, but works once we bgzip and idenx the bcf file -see below)
/software/team152/bcftools/bin/./bcftools view -R /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/snps_to_keep_after_hwe_check.txt \
/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.bcf \
-o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.bcf 

#bcftools complains that the bcf file wasn't zipped with bgzip. To overcome this, I use the -Ob and index commands of bcftools as follows:
/software/team152/bcftools/bin/./bcftools view  /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.nomiss.recode.rs1.bcf -I -Ob -o antegamisou.bcf  #i) bgzip
bcftools index antegamisou.bcf #ii) indexing
# Having finished these two steps, the command for removing the SNPs that are out of HWE can now be run:
# (Note that bcftools *does not* like a header in the file with the list of variants/regions, so make sure there is none.)
#bsub -q normal -J vcf_filter -G team152 -n1 -R "span[hosts=1] select[mem>4500] rusage[mem=4500]" -M4500 -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/vcf_hwefilt.log /software/team152/bcftools/bin/./bcftools view -R /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/snps_to_keep_after_hwe_check.txt -c1:minor -m2 -M2 -f PASS -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1b.bcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/antegamisou 

#03/Oct: I am also removing duplicate entries because they seem to give trouble downstream with multiallelic SNPs (I used R, handle_duplicates_in_tili_exomeseq.R)
/software/team152/bcftools/bin/./bcftools concat --rm-dup all -a -o antegamisou_rmvdup.bcf  antegamisou.bcf
/software/team152/bcftools/bin/./bcftools view -R /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/snps_to_keep_after_hwe_check_nodup.txt -c1:minor -m2 -M2  -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.bcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/antegamisou.bcf 
#the following should work (Petr's suggestion):
/software/team152/bcftools/bin/./bcftools view -T /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/snps_to_keep_after_hwe_check_nodup.txt -c1:minor -m2 -M2  -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.bcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/antegamisou.bcf 
/software/team152/bcftools/bin/./bcftools view -T /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/snps_to_keep_after_hwe_check_nodup_updated.txt -c1:minor -m2 -M2  -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.bcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/antegamisou.bcf 

#the following was just fiddling around for alternatives before the -T solution.
#bcftools view tili.poly.biallelic.qc1.bcf | grep -v '#'| cut -f3 | awk '{a[$0]++; if(a[$0]==2) print; if (a[$0]>=2) print}'
#bcftools view tili.poly.biallelic.qc1.bcf | awk '!uniq[$3]++' > tili.poly.biallelic.qc1.unique.bcf

#list of the remaining variants:
/software/team152/bcftools/bin/./bcftools view antegamisou.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > antegamisou.variants
/software/team152/bcftools/bin/./bcftools view antegamisou_rmvdup.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > antegamisou_rmvdup.variants
/software/team152/bcftools/bin/./bcftools view tili.poly.biallelic.qc1.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.poly.biallelic.qc1.variants
/software/team152/bcftools/bin/./bcftools view tili.poly.biallelic.qc1.unique.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.poly.biallelic.qc1.unique.bcf.variants
/software/team152/bcftools/bin/./bcftools view tili.poly.biallelic.qc1.T.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.poly.biallelic.qc1.T.variants
/software/team152/bcftools/bin/./bcftools view tili.poly.biallelic.qc1.T.nodup.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.poly.biallelic.qc1.T.nodup.variants

##Epitelous nomizw to tili.poly.biallelic.qc1.T.nodup.bcf einai to arxeio pou mporoume na paroume downstream.

#Once I have a file in bed format, I can use bcftools annotate to do that
#example: bcftools annotate -c CHROM,FROM,TO,ID -a my_ids.bed.gz   -o output.vcf  input.vcf.gz
#the following works(??), but make sure that the -a file is tab-delimited, bgzipped and tabixed
tabix -s1 -b2 -e2 chr${i}.1KG_b37.alleles #see tabixme.sh
bcftools annotate -c CHROM,POS,ID,REF,ALT -a /lustre/scratch115/projects/ibdgwas/aux_files/chr20.1KG_b37.alleles.gz   -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.bcf   /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.bcf
#we actually need to do the above across all chromosomes, so i have done so in a small script:
sh /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/update_rsIDs.sh

/software/team152/bcftools/bin/./bcftools view tili.poly.biallelic.qc1.T.nodup.rsID.final.bcf|grep -v '#'|awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' > tili.poly.biallelic.qc1.T.nodup.rsID.final.variants
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

#trying to make ensembl-vep work
## -example below works
./vep -i examples/homo_sapiens_GRCh37.vcf --cache --force_overwrite --merged --dir /software/team152/ensembl-vep/cache/ --port 3337

##this still creates a separate file, this time with a ton of information (Sift, Polyphen, CADD, AFs). Again it worked.
/software/team152/ensembl-vep/./vep -i examples/homo_sapiens_GRCh37.vcf --cache --force_overwrite --merged --dir /software/team152/ensembl-vep/cache/ --offline --canonical  \
  --stats_file  /software/team152/ensembl-vep/cache/vep.stats --output_file /software/team152/ensembl-vep/cache/homo_sapiens_GRCh37.vcf.annot  --sift p --poly p --af_1kg --af_gnomad\
 --regulatory --use_given_ref --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

# now for my actual file - make a full version as above, and also see if we can make an updated vcf automatically with a small number of fields.
tab-delimited
##first convert bcf file to vcf, for input to vep
/software/team152/bcftools/bin/./bcftools view -O v -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.bcf

##now run vep 
/software/team152/ensembl-vep/./vep -i /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf  --offline  \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  /software/team152/ensembl-vep/cache/vep.stats --use_given_ref \
--output_file /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf.annot --force_overwrite \
--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

#running vep on the farm (worked fine)
bsub -q normal -J veprun -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/vep.log \
 /software/team152/ensembl-vep/./vep -i /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf  --offline  \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --stats_file  /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vep.stats --use_given_ref \
--output_file /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf.annot --force_overwrite \
--sift p --poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz

DATA_DIR="/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/"
# the above produced a separate file, now I will work on incorporating some of these info into the vcf
##note that the fields that are added to the vcf are specified explicitly. 
##note that the *.annot file can be checked for more information
bsub -q normal -J veprun2 -G team152 -n1 -R "span[hosts=1] select[mem>5500] rusage[mem=5500]" -M5500 -o /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/vep2.log \
/software/team152/ensembl-vep/./vep -i /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf --cache --offline \
--merged --cache --dir /software/team152/ensembl-vep/cache/ --use_given_ref --vcf \
--output_file /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --force_overwrite --per_gene \
--sift b -poly p --af_1kg --af_gnomad --canonical --regulatory  --plugin CADD,/software/team152/ensembl-vep/cache/cadd_datasets/1000G_phase3.tsv.gz --symbol \
--fields Consequence,SYMBOL,CANONICAL,SIFT,PolyPhen,CADD_PHRED,EUR_AF,gnomAD_AMR_AF  

## At this stage we have a vcf file with basic site-filtering (missingness, HWE deviations) and annotation info. Now we can do some sample filtering too

#Start by looking into individuals with high missingness and remove any with missingness >5% 
vcftools --vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --missing-indv \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot
## 11/10/2017 Actually none had missingness above 2% 

#now check heterozygosity , relatedness, --indv-freq-burden 
vcftools --vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --het \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot

vcftools --vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --het \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot

vcftools --vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --indv-freq-burden \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot 

vcftools --vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --relatedness \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot 
 
vcftools --vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --relatedness2 \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot 

#relatedness results are actually different to the ones we had seen in the genotype data. Think whether it's worth checking only at thinned, common sites 
vcftools --vcf /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --relatedness --maf 0.05 \
--out /lustre/scratch115/projects/crohns/exome/TIH/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.common


##The list of related individuals doesn't overlap with that from the genotype data. Double check with Plink, after pruning.
DATADIR="/lustre/scratch115/projects/crohns/exome/TIH/"
plink --vcf ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --exclude ${DATADIR}aux_and_intermediate_files/high_ld.txt \
--indep-pairwise 1000 50 0.2 --maf 0.01 --out ${DATADIR}/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot
##then get IBS estimates
plink --vcf ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf  --exclude ${DATADIR}/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.prune.out --genome \
--out ${DATADIR}/exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.forIBS


##based on these, remove samples with high pairwise ibs ** for now just from the merged dataset. Then make a list of the ones to remove due to being PCA outliers, and THEN remove them
#from the respective datasets.
plink --allow-no-sex --bfile TILI_final_ctrls.qc3 --bmerge  TILI_final_cases.qc3 --exclude TILI_set_for_IBS.prune.out --make-bed --out TILI_merged.qc3.pruned
plink --allow-no-sex --bfile TILI_merged.qc3.pruned   --recode --remove high_ibs_TILI_qc3.samples.toremove.txt --out TILI_merged.qc3.pruned.filtered

#having a look to see why there are discrepancies between exomeseq and genotype data in terms of relatedness calculations.
#The ones that look like siblings in the genotype data are 109038428,1019749979
bcftools view -s 1109038428,1019749979 tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf|awk '{split($8,arr,";");split($10,arrA,":");split($11,arrB,":");{print $1,$2,$3,$4,$5,arr[1],arrA[1],arrB[1]}}'|more
##add a random sample in to see concordance
bcftools view -s 1109038428,1019749979,1100261686 tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf|awk '{split($8,arr,";");split($10,arrA,":");split($11,arrB,":");split($12,arrC,":");{print $1,$2,$3,$4,$5,arr[1],arrA[1],arrB[1],arrC[1]}}'|more

plink --bfile  TILI_merged.qc3.forQC --keep tmp.keep --recode vcf-iid --out exomeseq/check_dups --keep-allele-order #checking the two individuals that seem to be the same person

#there doesn't seem to be high concordance at all...
#check IBS based on these ~4750 SNPs alone.

plink --bfile  ${DATADIR}TILI_merged.qc3.forQC --allow-no-sex --extract ${DATADIR}exomeseq/snps_in_common_between_exomeseq_and_genotype_data.txt --genome --out ${DATADIR}exomeseq/genodata.shared.snps
plink --vcf ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --allow-no-sex --extract ${DATADIR}exomeseq/snps_in_common_between_exomeseq_and_genotype_data.txt --genome --out ${DATADIR}exomeseq/exomedata.shared.snps

#To assess whether an ID mismatch is likely, I will compare the genotypes of the two pairs of samples wth an IBS of ~1
#1008416606,1019747727,1019749979,1109038428
bcftools view -s 1008416606,1019747727,1019749979,1109038428 tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf -o exome.shared.highIBS
plink --bfile  ${DATADIR}TILI_merged.qc3.forQC --allow-no-sex  --keep ${DATADIR}exomeseq/tmp.keep2 --recode vcf-iid --out ${DATADIR}exomeseq/geno.shared.highIBS --keep-allele-order #checking the two individuals that seem to be the same person
cat ${DATADIR}exomeseq/exome.shared.highIBS.vcf |awk '{split($8,arr,";");split($10,arrA,":");split($11,arrB,":");split($12,arrC,":");split($13,arrD,":");{print $1,$2,$3,$4,$5,arr[1],arrA[1],arrB[1],arrC[1],arrD[1]}}' >exome.test

#From the comparison of the two files for the sites in common where the allele order was the same, it seems that there was very high agreement between the 
# samples  1008416606,1019747727 from the exome data and the samples 1019749979,1109038428 from the genotype data.

#from Gareths file:
IBD785	1008416606
IBD1693	1019747727
IBD1098	1019749979
IBD2704	1109038428

### Now that the mis -alignment looks very likely I will merge the exomeseq and genotype files for the 3276 sites in common and in the same order, then run IBS estimation
### to identify corresponding samples.
cat ${DATADIR}/exomeseq/variant_comparison_between_genotype_and_exomeseq_data.3.txt |awk '{print $1}' > ${DATADIR}/exomeseq/vars_in_common_and_same_order.txt
plink --vcf ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf --allow-no-sex --keep-allele-order --extract ${DATADIR}exomeseq/vars_in_common_and_same_order.txt --make-bed --out ${DATADIR}exomeseq/exomedata.shared.snps
cat ${DATADIR}exomeseq/exomedata.shared.snps.fam |awk '{print "ex"$1,"ex"$2,$3,$4,$5,$6}' > ${DATADIR}exomeseq/exomedata.shared.snps.fam.2
cat ${DATADIR}exomeseq/genodata.shared.snps.fam |awk '{print "geno"$1,"geno"$2,$3,$4,$5,$6}' > ${DATADIR}exomeseq/genodata.shared.snps.fam.2
#now prepare vcf files and then merge
plink --bfile ${DATADIR}exomeseq/genodata.shared.snps --keep-allele-order --recode vcf-iid --out ${DATADIR}exomeseq/genodata.shared.snps
plink --bfile ${DATADIR}exomeseq/exomedata.shared.snps --keep-allele-order --recode vcf-iid --out ${DATADIR}exomeseq/exomedata.shared.snps
#bgzip and tabix the files then run bcftools to merge:
bcftools merge  ${DATADIR}exomeseq/exomedata.shared.snps.vcf.gz  ${DATADIR}exomeseq/genodata.shared.snps.vcf.gz -O v -o ${DATADIR}exomeseq/merged.shared.snps.vcf --threads 3 --force-samples
#now everything is ready to check IBS.
plink --vcf ${DATADIR}exomeseq/merged.shared.snps.vcf --allow-no-sex --const-fid --genome --out ${DATADIR}exomeseq/merged.shared.snps
##Use file  ${DATADIR}exomeseq/merged.shared.snps to match exome and genotype sample IDs

###################################
## Gene-based tests using EPACTS ##
###################################

#For file preparation look at prepare_group_file_for_EPACTS.sh

#now prepare a pseudo ped file with 
 head -n500  ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf|grep QUAL|sed 's/\t/\n/g' >  ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf.samples
#with some basic R I prepared the pseudoped file: dummy_pheno.ped
EPACTS_DIR="/software/team152/EPACTS/" 
CHR=21
${EPACTS_DIR}bin/epacts group --vcf ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf.gz \
  --groupf ${DATADIR}exomeseq/tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf.${CHR}.annot.fc.genes.epacts --out ${DATADIR}exomeseq/test.gene.${CHR}.skato \
  --ped ${DATADIR}exomeseq/dummy_pheno.ped --max-maf 0.05 \
   --pheno DISEASE --test skat --skat-o --run 2


  