

import os
import sys
import subprocess as sp


my_genelist = ['OG0003058.fasta', 'OG0003060.fasta']
for myfasta_name in my_genelist:	

	os.system('rm firstline.ctl')
	os.system('rm lastline.ctl')
	# strip out OG name from mygene
	ret = myfasta_name.split('.fasta')	
	myorthogroup_name = ret[0]
	myseqfile = 'seqfile = ' + '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/JGI_files_PAML/2_Orthogroups_trial/new_fasta_and_alignment_file/' + myorthogroup_name + 'msf' + '\n'
	my_alt_outputfile = 'outfile = ' + '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/JGI_files_PAML/2_Orthogroups_trial/output_files' + myorthogroup_name + '.altout' + '\n'
	my_nul_outputfile = 'outfile = ' + '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/JGI_files_PAML/2_Orthogroups_trial/output_files' + myorthogroup_name + '.nulout' + '\n'	#\n used to return
	f = open('firstline.ctl','a')
	f.write(myseqfile)
	f.close()
	f = open('lastline.ctl', 'a')
	f.write(my_alt_outputfile)
	f.close()
	os.system('cat firstline.ctl lastline.ctl default_revised_alt.ctl > final.ctl')
	os.system('codeml final.ctl')
