import sys
import csv
import pandas as pd
with open('Projorthogroups.txt') as proj1file:
	Proj1reader = csv.reader(proj1file, delimiter =' ')
	for line in Proj1reader:
		print(line[0])
