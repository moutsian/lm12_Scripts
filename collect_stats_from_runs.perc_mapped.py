#!/usr/bin/python
import os, string
from fnmatch import fnmatch
import sys
run_to_check = sys.argv[1]
print "Getting stats for run "+run_to_check
root = '/lustre/scratch115/projects/paxgene/SalmonOutput_newruns/'
name_pattern = "meta_info.json"
path_pattern = "*sorted*"
runfiles=[]
filecheck="*"+run_to_check+"*"
for path, subdirs, files in os.walk(root):
    for name in files:
	if fnmatch(path, path_pattern):
        	if fnmatch(name, name_pattern):
			fullpath=os.path.join(path,name)
			if fnmatch(fullpath,filecheck):
				runfiles.append(fullpath)

out_data=[]
for myfile in runfiles:
	with open(myfile,'r') as f:
		for line in f.readlines():
			li=line.lstrip()
			if  li.startswith('"percent_mapped'):
				num_assigned_fragments= li.split()[1]
				run=myfile.split('/')[7].split('_')[0]
				lane=myfile.split('/')[7].split('_')[1]
				sample=myfile.split('/')[7].split('_')[2]
				out_data.append ("{}\t{}\t{}\t{}".format(run, lane, sample, num_assigned_fragments[:-1]))			
f.close()
#now save output file
OUTFILE=root+"run" + run_to_check + ".sorted.percent_mapped_fragments.txt"
with open(OUTFILE,'w') as outf:
	for item in out_data:
		outf.write("{}\n".format(item,"\n"))
outf.close()
