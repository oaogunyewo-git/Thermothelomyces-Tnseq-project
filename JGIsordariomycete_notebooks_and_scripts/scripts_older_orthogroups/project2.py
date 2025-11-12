import sys
import csv
import pandas as pd

with open('singlecopyorthogroups.txt') as proj2file:
	proj2reader = csv.reader(proj2file, delimiter =' ')
	for line in Proj2reader:
		print(line[0])
