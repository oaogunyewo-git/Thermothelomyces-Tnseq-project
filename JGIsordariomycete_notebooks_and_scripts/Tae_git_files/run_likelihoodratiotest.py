#load required packages
from sys import argv
import os
from scipy.stats import chi2
import statsmodels.stats.multitest as smm
import pandas as pd
from Bio.Phylo.PAML import codeml
import datetime

#################################
#this script takes the directory that your codeml output files are found in, and outputs a single time stamped file with the LRT results
#input is: python /directory/run_likelihoodratiotest.py /output_directory/ /codeml_results/directory/ species
#################################

#function to read in codeml files
def parse_paml(filename):
	invalid = 0
	fileID = filename.split("/")[-1]
	gene = fileID.split(".")[0]
	model = fileID.split("paml")[-1].split("out")[0]
	#initiate lnL and np as None, in case they are invalid
	lnL = None
	np = None
	try:
		results = codeml.read(filename)

		#in each model, filter for unreliable omega values
		if model == "alt":
			if sum([x != 999.0 for x in results["NSsites"][0]["parameters"]["omega"]]) == 2:
				lnL = results["NSsites"][0]["lnL"]
				np = len(results["NSsites"][0]["parameters"]["parameter list"].split(" "))
		elif model == "nul":
			if results["NSsites"][0]["parameters"]["omega"] != 999.0:
				lnL = results["NSsites"][0]["lnL"]
				np = len(results["NSsites"][0]["parameters"]["parameter list"].split(" "))

	except: #sometimes PAML couldn't run the alt or nul, so I have this check
		invalid += 1

	#print invalid/2, "invalid PAML files."
	return gene, model, lnL, np

#run likelihood ratio test
def run_lrtest(dat):
	stamp = str(datetime.datetime.now())
	stamp2 = stamp.replace(" ","_")
	final = {}
	wf = open(outdir+stamp2+"_codeml.likelihoodratios","w")
	wf.writelines("Genes\tlstat\tddof\n")
	for i in dat:
		if len(dat[i]) != 2:
			continue
		if i not in final:
			final[i] = []
		lstat = 2*(dat[i]["alt"][0] - dat[i]["nul"][0])
		ddof = dat[i]["alt"][1] - dat[i]["nul"][1]
		final[i] = [lstat,ddof]
		wf.writelines(i+"\t"+str(lstat)+"\t"+str(ddof)+"\n")
	wf.close()

	gene_list = []
	pval_list = []
	
	for i in final:
		chisq = chi2.sf(final[i][0],final[i][1])
		gene_list.append(i)
		pval_list.append(chisq)
	
	#multiple testing correction
	fdrbh_output = smm.multipletests(pval_list,method="fdr_bh")
	adj_pvals = fdrbh_output[1].tolist()
	zipped_pvals = zip(pval_list, adj_pvals)
	zipped_gene_info = zip(gene_list, zipped_pvals)
	unpacked_gene_info = [(gene, pval, adjusted_pval) for (gene, (pval, adjusted_pval)) in zipped_gene_info]
	pval_df = pd.DataFrame(unpacked_gene_info, columns=['gene', 'pval', 'adjusted_pval'])

	pval_df.sort_values(by='pval', inplace=True)
	pval_df.to_csv(outdir+stamp2+'_codeml.likelihoodratio_chisqtest', sep='\t', index=False)

#start script here
if __name__ == "__main__":

	outdir = argv[1]
	paml_dir = argv[2]
	species = argv[3]
	paml_files = os.listdir(paml_dir)
	print len(paml_files)/2, "genes with PAML output for", species

	all_dat = {}
	for paml_file in paml_files:
		Gene,Model,LNL,NP = parse_paml(paml_dir+paml_file)
		if LNL == None or NP == None:
			continue
		if Gene not in all_dat:
			all_dat[Gene] = {}
		all_dat[Gene][Model] = [LNL,NP]

	print len(all_dat), "genes with real dN/dS values for", species

	run_lrtest(all_dat)
	print "Done for", species, "!"
