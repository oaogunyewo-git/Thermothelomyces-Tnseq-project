import sys
import csv
import pandas as pd

# store the list of all single-copy orthogroup names
single_copy_ortho = []
with open('singlecopyorthogroups.txt') as proj2file:
	proj2reader = csv.reader(proj2file, delimiter =' ')
	for line in proj2reader:
		single_copy_ortho.append(line[0]) # the append here implies that you are asking the program to read the orthogroups and and stick the corresponding name to it

# MAKE A LIST OF SPECIES NAMES YOU ACTUALLY ARE INTERESTED IN
MySpeciesList = ['apoper-1', 'canare-1', 'chacoc-1','chathe-1', 'corsep-1', 'crahot-1', 'crathe-1', 'dicfun-1', 'gelunc-1', 'mycthe-1', 'neudis4-', 'parapp-1', 'parhyr-1', 'podans-0', 'podapp-1', 'podcur-1', 'thehet-1', 'theter-1', 'theher-1', 'thethe-1', 'triant-1']	# MAKE A LIST OF SPECIES NAMES YOU ACTUALLY ARE INTERESTED IN

with open('Orthogroups.txt') as proj3file: # the ''   here makes this a string
        Proj3reader = csv.reader(proj3file, delimiter =' ')
        for line in Proj3reader:
	 	my_orthogroup = line[0]
		ret = my_orthogroup.split(":")				#[:] was used because the orthogroups were separated by :
		my_orthogroup_no_colon = ret[0]
		if my_orthogroup_no_colon in single_copy_ortho:
               		genes_in_my_orthogroup = line[1:]
			orthogroup_file_name = my_orthogroup_no_colon + ".fasta"
			for gene in genes_in_my_orthogroup:
				my_header = ">" + gene
				all_the_elements_name = gene.split("_") 	#_ was used for splitting because the specie and gene names were separated by [_]
				sp_name = all_the_elements_name[0]		# space cannot be included in a variable
				if(sp_name in MySpeciesList):  			# PUT IN AN IF-IN STATEMENT TO ASK IF SP_NAME IS AMONG THOSE OF INTEREST
					my_sp_filename = 'individuals_fastas_CDS_nucl/' + sp_name + '.fasta' 	#new string of letters which is the sp_name.fasta; rem to provide the directory for it to be read
					found_my_gene = 0 							#Read the dataset and nothing found
					with open(my_sp_filename) as Proj5file: 				# I already have a string here so there is no need to keep the name in quote
						for line_seq in Proj5file:
							if(my_header in line_seq): 				# if this is the header line for my gene
								found_my_gene = 1 				#here, it recongizes that the first gene ID has been found 
								line_seq_data = line_seq
							else:
								if(found_my_gene == 1):					# if we are past the header line for my gene 
									if(">" in line_seq):	# if we are at a new header line for the next gene, imoplies it has already passed it (data shoukd be collected before this) 
										break  					#The break is the end of the loop
									else:
										line_seq_data = line_seq_data + line_seq
#							if(line_seq == my_header):	# == used as an evaluation statement to ask if a defined statement is true or not
					gene_to_orthogroup_file = line_seq_data
					orthogroup_file = open(orthogroup_file_name, "a")	#'a' is for appending to file multiple times
					orthogroup_file.write(line_seq_data)
					orthogroup_file.close()
				
