from sys import argv
from Bio import AlignIO
import os
import numpy as np
import pickle

#######################################################################################################
#converts MUSCLE alignment output into phylip format for PAML
#adapted from code given by Claire Dubin
#######################################################################################################

#read in FASTA files, removes those with a premature stop codon and too many gaps in the alignment
def read_files(filename):

	fileID = filename.split("/")[-1].split(".")[0]
	geneID = fileID.split("_")[0]
	stops = ["TAA","TAG","TGA"]

	seqdat = AlignIO.read(filename,"fasta")

	tempdict = {}
	for seq in seqdat:
		tempseq = seq.seq
		if float(tempseq.count("-"))/len(tempseq) > 0.1:
			continue
		stop = False
		for i in range(0,len(tempseq),3):
			if tempseq[i:i+3] in stops:
				stop = True
				break
		if stop == False:
			tempdict[seq.id] = str(seq.seq)

	return geneID, tempdict

#convert to phylip format
def make_phylip(geneID,genedat):
	for dat in genedat:
		for j in dat:
			if j.startswith("EN"):
				mmid = j
				mmseq = dat[j]
				mmlen = len(mmseq)
			elif j.startswith("MGP_S"):
				msid = j
				msseq = dat[j]
				mslen = len(msseq)
			elif j.startswith("MGP_C"):
				mcid = j
				mcseq = dat[j]
				mclen = len(mcseq)

		if mmlen == mslen == mclen:
			#matched += 1
			align_len = mmlen

			#fileID = [i,".".join(idat.keys())]
			fileID = "_".join([geneID,mmid,msid,mcid])

			wf = open(outdir+fileID+".phylip","w")
			wf.writelines("3 "+str(align_len)+"\n")

			wf.writelines("mm\n")
			for l in range(0,len(mmseq),60):
				wf.writelines(mmseq[l:l+60]+"\n")

			wf.writelines("ms\n")
			for l in range(0,len(msseq),60):
				wf.writelines(msseq[l:l+60]+"\n")

			wf.writelines("mc\n")
			for l in range(0,len(mcseq),60):
				wf.writelines(mcseq[l:l+60]+"\n")

			wf.close()

if __name__ == "__main__":
	outdir = argv[1]
	muscle_files = argv[2:]

	all_dat = {}
	all_ct = 0
	added = 0

	for muscle_file in muscle_files:
		all_ct += 1
		gene_ID, temp_dict = read_files(muscle_file)
		if len(temp_dict) == 3: #make sure all three species alignments are acceptable
			if gene_ID not in all_dat:
				all_dat[gene_ID] = []
			all_dat[gene_ID].append(temp_dict)
			added += 1

	print all_ct, "total transcript alignments,", added, "with less than 10% gaps and no in frame stop codons,", float(100*added)/all_ct, "% of total."

	for gene in all_dat:
		make_phylip(gene,all_dat[gene])

	print "Done!"

