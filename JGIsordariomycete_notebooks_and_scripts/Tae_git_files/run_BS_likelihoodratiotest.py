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
	#check for valid lnL values
	lnL = None
	try:
		results = codeml.read(filename)
		lnL = results["NSsites"][2]["lnL"]

	except:
		invalid += 1

	#print invalid/2, "invalid PAML files."
	return gene, model, lnL

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
		lstat = 2*(dat[i]["alt"] - dat[i]["nul"])
		ddof = 1
		final[i] = [lstat,ddof]
		wf.writelines(i+"\t"+str(lstat)+"\t"+str(ddof)+"\n")
	wf.close()

	gene_list = []
	pval_list = []
	
	for i in final:
		chisq = chi2.sf(final[i][0],final[i][1])
		gene_list.append(i)
		pval_list.append(chisq/2)

	#multiple testing correction
	fdrbh_output = smm.multipletests(pval_list,method="fdr_bh")
	adj_pvals = fdrbh_output[1].tolist()
	zipped_pvals = zip(pval_list, adj_pvals)
	zipped_gene_info = zip(gene_list, zipped_pvals)
	unpacked_gene_info = [(gene, pval, adjusted_pval) for (gene, (pval, adjusted_pval)) in zipped_gene_info]
	pval_df = pd.DataFrame(unpacked_gene_info, columns=['gene', 'pval', 'adjusted_pval'])

	pval_df.sort_values(by='pval', inplace=True)
	pval_df.to_csv(outdir+stamp2+'_codeml.likelihoodratio_chisqtest', sep='\t', index=False)

if __name__ == "__main__":

	outdir = argv[1]
	paml_dir = argv[2]
	species = argv[3]
	paml_files = os.listdir(paml_dir)
	print len(paml_files)/2, "genes with PAML output for", species

	all_dat = {}
	for paml_file in paml_files:
		Gene,Model,LNL = parse_paml(paml_dir+paml_file)
		if LNL == None:
			continue
		if Gene not in all_dat:
			all_dat[Gene] = {}
		all_dat[Gene][Model] = LNL

	print len(all_dat), "genes with real dN/dS values for", species

	run_lrtest(all_dat)
	print "Done for", species, "!"

