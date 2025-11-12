# will read in a directory of phylip alignments, then run PAML codeml on each alignment and put the results file for each gene into a new directory.

# Directories must be input using the full path name!

# adapted from code given by Claire Dubin


import os
from sys import argv
import subprocess as sp


def Write_ctl_file(file,model):

	gene = file.split(".")[0]

	#makes a temp control file every time to avoid editing main control file
	template_ctl_file = open('/bigrock_home/tkang/paml2022/branchsite/ctl/'+species+'/paml_'+model+'.ctl')

	wf = open('/bigrock_home/tkang/paml2022/branchsite/ctl/'+species+'/tempctl/temp_ctl_'+model+'.ctl','w')

	for line in template_ctl_file:
		if "seqfile" in line:
			wf.writelines("\tseqfile = "+indir+file+"\n")
		elif "outfile" in line:
			wf.writelines("\toutfile = "+outdir+gene+".paml"+model+"out\n")
		else:
			wf.writelines(line)

	wf.close()

def Run_codeml(model):

	cmd = ['codeml', '/bigrock_home/tkang/paml2022/branchsite/ctl/'+species+'/tempctl/temp_ctl_'+model+'.ctl']
	sp.call(cmd)

if __name__ == "__main__":
	indir = argv[1]
	outdir = argv[2]
	species = argv[3]

	all_in_files = os.listdir(indir)

	for file in all_in_files:

		Write_ctl_file(file,"alt")
		Write_ctl_file(file,"nul")
		Run_codeml("alt")
		Run_codeml("nul")

