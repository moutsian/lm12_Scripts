#!/usr/bin/python
DIR = "/lustre/scratch115/projects/crohns/exome/TIH/exomeseq/"
FILE = "tili.poly.biallelic.nomiss.recode.rs1.variants" #contains the correct IDs
VARFILE = DIR+FILE
CHROMS=[]
POSITIONS=[]
IDs=[]
POS_DICT={}
with open(VARFILE,'r') as f:
	for line in f.readlines():
		li=line.lstrip()
		if not li.startswith("#"):
			CHR, POS, ID, REF, ALT, SCORE, QUAL, INFO = li.split()
			CHROMS.append(CHR)
			POSITIONS.append(POS)
			IDs.append(ID)
			POS_DICT[CHR+POS]=ID
#print IDs[:10]
f.close()

HWEFILENAME="tili.poly.biallelic.nomiss.recode.rs1.ctrls.bcf.hwe"
HWEFILE=DIR+HWEFILENAME
OUTFILE=DIR+"tili.poly.biallelic.nomiss.recode.rs1.outofHWE.txt"
OUTOFHWE=[]
with open(HWEFILE,'r') as f:
	f.readline() #skip header
	for line in f.readlines():
		pHWE=line.split()[5]
		Chr=line.split()[0]
		Pos=line.split()[1]
		Chr_Pos=Chr+Pos
		if float(pHWE)<1e-08:
			if Chr_Pos in POS_DICT:
			#idx=POSITIONS.index(line.split()[1])
			#if CHROMS[idx]==line.split()[0]:
			#OUTOFHWE.append(IDs[idx])
				OUTOFHWE.append(POS_DICT[Chr_Pos])
f.close()
#now save OUTOFHWE list to file
with open(OUTFILE,'w') as outf:
	for item in OUTOFHWE:
		outf.write("%s\n" % item)
outf.close()
