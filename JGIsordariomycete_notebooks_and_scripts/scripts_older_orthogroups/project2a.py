import sys
import csv
import pandas as pd

# store the list of all single-copy orthogroup names
single_copy_ortho = []
with open('singlecopyorthogroups.txt') as proj2file:
	proj2reader = csv.reader(proj2file, delimiter =' ')
	for line in proj2reader:
		single_copy_ortho.append(line[0])

with open('Projorthogroups.txt') as proj1file:
        Proj1reader = csv.reader(proj1file, delimiter =' ')
        for line in Proj1reader:
		my_orthogroup = line[0]
		ret = my_orthogroup.split(":")
		my_orthogroup_no_colon = ret[0]
		if(my_orthogroup_no_colon in single_copy_ortho):
			print(my_orthogroup_no_colon)
#                mylist = line[1:]
