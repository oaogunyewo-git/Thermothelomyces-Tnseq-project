import pandas as pd
import os
import subprocess
#import matplotlib.pyplot as plt
#makes viewing pandas tables better
#pd.set_option('display.max_colwidth', 0)
#change this to where your FASTQs are
#How big is each FASTQ file?
FASTQ_directory = '/usr2/people/shollyt22/shollyt22/TnSeq_Results_may_2023/TnSeq_M003593'


fastq_files = ['OORB001_S1_L001_R1_001.fastq.gz', 'OORB001_S1_L001_R2_001.fastq.gz']

for f in fastq_files:

    size = !du -sh {f}
    size = size[0].split('\t')[0]
    
	print('{}: {}'.format(f, size))

