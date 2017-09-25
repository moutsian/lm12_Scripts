##

LOWFREQ=1;
trait="uc"
utrait=toupper(trait);
ltrait=tolower(trait);
#a HWE threshold of 1e-07 has been used. Note that what is reported in the HWE column is -log10P for HWE.
GWAS1out=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS1.outofHWE.controls.txt",head=T)
GWAS2out=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS2.outofHWE.controls.txt",head=T)
GWAS3out=read.table("/lustre/scratch115/teams/anderson/ibd_conditional/snp_qc/snp_stats.GWAS3.outofHWE.controls.txt",head=T)
ourdata=read.table("~lm12/IBD_conditional/step2.NOVEL.from_all_papers.with_ichip_info.txt",head=T)
