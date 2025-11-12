import sys
import csv
import pandas as pd
with open('complete_orthogroups.txt') as proj3file:
	Proj3reader = csv.reader(proj3file, delimiter =' ')
	for line in Proj3reader:
		print(line [0])

