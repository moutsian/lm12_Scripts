##Prepare group file for EPACTS based on ensembl-vep annotation

##Input
annofile="tili.poly.biallelic.qc1.T.nodup.rsID.final.vcf.annot"
vcffile="tili.poly.biallelic.qc1.T.nodup.rsID.annot.vcf"
chrom=20

# ** Note: If annotation list above is modified, the code will need to be edited **
annotations_list=("5_prime_UTR_variant"\
 "3_prime_UTR_variant" "coding_sequence_variant" "downstream_gene_variant" "frameshift_variant" "incomplete_terminal_codon_variant" \
 "inframe_deletion" "inframe_insertion" "intergenic_variant" "intron_variant" "mature_miRNA_variant"\
 "missense_variant" "non_coding_transcript_exon_variant" "NMD_transcript_variant" "protein_altering_variant" "regulatory_region_variant"\
 "splice_acceptor_variant" "splice_donor_variant" "splice_region_variant" "start_lost" "start_retained_variant"\
 "stop_gained" "stop_lost" "stop_retained_variant" "synonymous_variant" "TF_binding_site_variant" "upstream_gene_variant")


#The functional coding we used for IBDseq paper included the following:
# o	Functional coding (‘fc’) includes the following annotations from Ensembl’s variant effect predictor (VeP): 
#frameshift_variant, stop_gained, initiator_codon_variant, splice_donor_variant, splice_acceptor_variant, missense_variant, stop_lost, inframe_deletion, and inframe_insertion.

#below I've prepared the list for this annotation -note that we could modify and that "initiator_codon_variant" seems to be missing
fc="${annotations_list[4]}\\|${annotations_list[21]}\\|${annotations_list[22]}\\|${annotations_list[16]}\\|${annotations_list[17]}\\|${annotations_list[18]}\\|${annotations_list[11]}\\|${annotations_list[7]}"

#save a file for further processing 
grep ${fc} ${annofile} |awk -v CHR=${chrom} '{split($2,arr,":");{if(arr[1]==CHR){print $0}}}'> ${annofile}.${chrom}.annot.fc.tmp

## Below is an example line from the EPACTS accepted format. 
## ARFRP1  20:62331850_A/C 20:62331871_C/T 20:62332005_C/T
##Note that the variants are ordered (should be already fine in our case) and that the two alleles (REF/ALT) are required.
## I will have to get these from the vcf file
grep -v '#' ${vcffile} |awk -v CHR=${chrom} '{if($1==CHR){print $1,$2,$3,$4,$5}}' > ${vcffile}.${chrom}.variants 

#now load this into prepare_group_file_for_EPACTS_Rpart.R to produce the group file

#note that EPACTS needs a pseudo ped file to use for phenotype info