#!/usr/bin/python
DIR = "/lustre/scratch115/projects/crohns/exome/TIH/"
FILE = "high_ibs_TILI_qc1d.txt" #contains the pairs of samples with high IBS
VARFILE = DIR+FILE
FID1=[]
IID1=[]
FID2=[]
IID2=[]
with open(VARFILE,'r') as f:
	f.readline() #skip header
	for line in f.readlines():
		li=line.lstrip()
		if not li.startswith("#"):
			fid1, iid1, fid2, iid2, rest = li.split('\t',5)
			FID1.append(fid1)
			IID1.append(iid1)
			FID2.append(fid2)
			IID2.append(iid2)
#print FID1[:]
f.close()

#get missingness per sample in
MISSFILENAME="TILI_merged.qc1d.imiss"
MISSFILE=DIR+MISSFILENAME
OUTFILE=DIR+"high_ibs_TILI_qc1d.samples.toremove.txt"
MISSDATA={}
with open(MISSFILE,'r') as f:
	f.readline() #skip header
	for line in f.readlines():
		fid=line.split()[0]
		miss=line.split()[5]
		if (fid in FID1) or (fid in FID2) :
			MISSDATA[fid]=miss
f.close()

#now pick a sample to remove from each high IBS pair based on missingness
with open (OUTFILE,'w') as outf:
	for i, item in enumerate(FID1):
		if MISSDATA[item] > MISSDATA[FID2[i]] :
			outf.write("{}\t{}".format(item, IID1[i]))
		else:
			outf.write("{}\t{}".format(FID2[i], IID2[i]))
		outf.write("\n");		 
outf.close()
