import sys
from sys import argv
from Bio import SeqIO

#######################################################################################################
#takes in the CDS sequences from gffread, organizes them into a single FASTA file for all species
#just looking at sequences from Mus musculus, M. spretus, and M. caroli
#######################################################################################################

#parse file of homologous genes from ensembl
def parse_homs(filename,sep="\t"):
	homsdict = {}
	with open(filename) as f:
		f.readline()
		for line in f:
			linedata = line.strip().split(sep)
			if linedata[5] == "ortholog_one2one":
				homsdict[linedata[0]] = linedata[1]
	return homsdict

#parse CDS data from gffread
def parse_cds(filename):
	cds_seqs = SeqIO.parse(filename,"fasta")
	cdsdict = {}

	stops = ["TAA","TAG","TGA"]

	for cds in cds_seqs:
		cdsid = cds.id
		cdsseq = str(cds.seq)
		if cdsseq[-3:] in stops:
			cdsseq = cdsseq[:-3]
		cdsdict[cdsid] = str(cdsseq)

	return cdsdict

#parse transcript annotations from ensembl
def parse_transcripts(filename,homsdict,sep="\t"):
	tdict = {}
	with open(filename) as f:
		f.readline()
		for line in f:
			linedata = line.strip().split(sep)
			if len(linedata) != 8:
				continue
			if linedata[0] in homsdict:
				geneID = homsdict[linedata[0]]
				if geneID not in tdict:
					tdict[geneID] = []
				ms_trans = linedata[5].replace("_P","_T")
				mc_trans = linedata[7].replace("_P","_T")
				tdict[geneID].append([linedata[2],ms_trans,mc_trans])
	return tdict

#organize CDS data
def get_seqs(tdict,mmcds,mscds,mccds):
	totalct = 0
	same = 0
	for i in tdict:
		for j in tdict[i]:
			mmID,msID,mcID = j

			if mmID in mmcds and msID in mscds and mcID in mccds:
				mmseq = mmcds[mmID]
				mmcds.pop(mmID)

				msseq = mscds[msID]
				mscds.pop(msID)

				mcseq = mccds[mcID]
				mccds.pop(mcID)

				totalct += 1

				wf = open(outdir+i+"_"+mmID+"_"+msID+"_"+mcID+".fa","w")
				wf.writelines(">"+mmID+"\n"+mmseq+"\n>"+msID+"\n"+msseq+"\n>"+mcID+"\n"+mcseq+"\n")
				wf.close()

	print len(tdict), "genes with homology and stable transcript IDs."
	print totalct, "homologous transcripts."


if __name__ == "__main__":

	outdir = argv[1]
	homs_file = argv[2]
	transcripts_file = argv[3]
	cds_files = argv[4:]

	homs_dict = parse_homs(homs_file)
	t_dict = parse_transcripts(transcripts_file,homs_dict)

	for cds_file in cds_files:
		fileID = cds_file.split("/")[-1]
		if fileID.startswith("GRC"):
			mm_cds = parse_cds(cds_file)
			print len(mm_cds), "CDS seqs in M. mus GTF."
		elif fileID.startswith("SPR"):
			ms_cds = parse_cds(cds_file)
			print len(ms_cds), "CDS seqs in M. spr GTF."
		elif fileID.startswith("CAR"):
			mc_cds = parse_cds(cds_file)
			print len(mc_cds), "CDS seqs in M. car GTF."

	get_seqs(t_dict,mm_cds,ms_cds,mc_cds)

