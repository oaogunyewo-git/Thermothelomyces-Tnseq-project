#####  pipeline for filtering raxml-ng gene trees with treeshrink
#
# this script is aimed to to take all the gene trees obtained with the pipeline A
# it groups them ina single file, and launches the TreeShrink program on this file
#
# if at least one individual is to remove from a gene tree,
# this script launches the pipeline A for the concerned orthogroup,
# giving the information of individuals to remove
# else, it puts the raxml-ng gene tree in a directory for being used (later) in the pipeline C
# 
#
# by Lucas Bonometti
#
#
#####  what does the pipeline do ?
#
# 1) first it creates its own environment, building the directories it needs if they do not yet exist
#
# 2) then it takes all the gene trees  as  workingdirectory/raxmlng/orthogroup/orthogroup.raxml.support
#    and it puts them in the same directory (tree2shrink)
#    and it concatenates them in a single multi-tree file (2shrink.trees) in a specific order
#
# 3) it launches TreeShrink with the following parameters
#         -b 20        (an individual is not removed if it reduces the tree size by less than 20 %)
#         -q 0.02      (an individual is removed from a gene tree if its value on this tree
#                       is on the 2% highest of the individual on the dataset)
#
# 4) if no individual is to be removed from a gene tree, the script puts the gene tree
#    in the directeory for the pipeline C (shrink2collapse)
#
# 5) if at least one individual is to be removed from a gene treee,
#    the script writes the concerned individual(s) in toremove_orthogroup_number.txt
#    and launches the pipeline A for the first 400 orthogroups to re-analyse
#    and it prepares the launchers for the following orthogroups (400 by 400)
#
#
##### points of warning
#
# 1) do not modifie the subdirectories in the working directory between pipelines A and B
#
# 2) this pipeline must be launched only once in an analysis
#    (launching TreeShrink on already skhrinked trees would remove more individuals than necessary)
#
#
###### libraries

import os
import subprocess
import sys
import argparse


##### creation of an argument parser
#
# those are the arguments needed to the script

parser = argparse.ArgumentParser(description="script filtering gene trees from pipeline A using TreeShrink")

parser.add_argument("-d","--work_dir",
                    help="# the path to the working directory (default : ./ )",
                    required=False, default="./")

args = parser.parse_args()



##### no home-made function is used on the script



##### there is where the script begins


### the first part is aimed to creates its own environment to the script


d=args.work_dir            # "d" is the path to the working directory
if d[-1]!="/": d=d+"/"

contend = []             # contend (with a "d") is the content in directories of the working directory
for item in os.listdir(d):
    if os.path.isdir(d+item):
        contend.append(item)


nlog=1
while "log_pipeB_"+str(nlog)+".log" in os.listdir(d):
    nlog+=1

with open(d+"log_pipeB_"+str(nlog)+".log","w") as log :
    
    # log is the log file, in which we will write allong the script
    
    log.write("** initialization done\n")
    
    drax=d+"raxmlng/"      # drax is the directory containing all the trees obtained with the pipeline A
    
    listog=[]      # listog is an ordered list of all the orthogroups with a suitable tree (OG.raxml.support)
    
    subdrax=os.listdir(drax)   # subdrax is a list of the orthogorups with a directory in drax
                               # (the gene trees are in those sub-directories)
    
    for OG in subdrax:
        if OG+".raxml.support" in os.listdir(drax+OG):
            listog.append(OG)
    
    
    
    
    
    dt2s=d+"tree2shrink"         # "dt2s" is the directory containing the trees to input in TreeShrink
    if "tree2shrink" not in contend:
        sp = subprocess.check_call(["mkdir",dt2s])
    elif len(os.listdir(dt2s))>0:
        sp = subprocess.Popen("rm "+dt2s+"/*",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    dt2s=dt2s+"/"
    
    
    
    
    ds2c=d+"shrink2collapse"    # "ds2c" is the directory containing the trees to input in the pipeline C
    if "shrink2collapse" not in contend:
        sp = subprocess.check_call(["mkdir",ds2c])
    ds2c=ds2c+"/"
    
    
    
    
    dsh=d+"shrinker"     # "dsh" is the directory in which TreeShrink is launch (with input and output data)
    if "shrinker" not in contend:
        sp = subprocess.check_call(["mkdir",dsh])
    elif "2shrink.trees" in os.listdir(dsh):
        sp = subprocess.call(["rm",dsh+"/2shrink.trees"])
    nout=1
    while "output_"+str(nout) in os.listdir(dsh):
        nout+=1
    sp = subprocess.check_call(["mkdir",dsh+"/output_"+str(nout)])
    dsh=dsh+"/"
    
    
    log.write("** environment crated\n")
    
    
    
    ##### in this part, we copy all the single gene trees in t2s
    ##### and we concatenate all the gene trees in a single file 2shrink.trees in dsh
    
    with open(dsh+"2shrink.trees","w") as OUT:
        for OG in listog:
            sp = subprocess.check_call(["cp",drax+OG+"/"+OG+".raxml.support",dt2s+OG+".nw"])
            with open(drax+OG+"/"+OG+".raxml.support",'r') as IN:
                OUT.write(IN.read())
    
    if len(os.listdir(dt2s))!=len(listog):
        log.write("\nERROR : problem as we copied the files in tree2shrink\n")
        sys.exit()
    
    log.write("** single gene trees copied in tree2shrink\n")
    log.write("** 2shrink.trees created\n")
    
    ##### then we launch treeshrink with the following parameters
    #             -b 20
    #             -q 0.02
    
    sp = subprocess.Popen("module load treeshrink/1.3.7 && module load r/3.6.3 "+
                          " && run_treeshrink.py -t "+
                                 dsh+"2shrink.trees -b 20 -q 0.02 -o "+
                                 dsh+"output_"+str(nout)+" -O shrinked -m per-species > "+
                                 dsh+"/output_"+str(nout)+"/log_treeshrink_"+str(nout)+".log",
                                 shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    
    if list(sp.communicate())[1].decode()!='':
        log.write("\nERROR as doing treeshrink\n\n")
        log.write(list(sp.communicate())[1].decode())
        sys.exit()
    
    
    log.write("** TreeShrink finished\n\n")
    
    
    ##### then we analyse the output of TreeShrink
    
    
    tolaunch=[]    # tolaunch is the list of orthogroups for which the pipeline A needs to be re-launched
    
    
    
    with open(dsh+"output_"+str(nout)+"/shrinked.txt",'r') as IN:
        # we read the output from TreeShrink containing the individuals to remove
        # 1 line = 1 orthogroup
        
        rankog=0   # rankog is the rank of the concerned orthogroup in listog
                   # rankog = 0 for line 1,  1 for line 2,  2 for line 3,  etc.
        
        
        for line in IN:
            
            OG=listog[rankog]      # OG is the concerned orthogroup
            
            
            # toremove is the list of the individuals to remove in the OG tree
            toremove=line.strip().split()
            
            
            # toremove is of length 0 if no individual is to be remove
            # so we just copy the gene tree of OG in s2c for the pipeline C
            if len(toremove)==0:
                
                log.write(OG+"\tOK\t")
                
                if OG+"_skrinked.nw" not in os.listdir(ds2c):
                    sp = subprocess.check_call(["cp",dt2s+OG+".nw",ds2c+OG+"_skrinked.nw"])
                    
                    log.write("copied\n")
                
                else :
                    
                    log.write("no-copied\n")
                    
                    
                #### /!\ warning : this script does not launch automatically the pipeline C
            
            
            
            
            # if the is at least one individual to remove from the orthogroup OG
            # we write a toremove file for the orthogroup with the individuals to remove
            # and we add the orthogroup to the tolaunch list
            else:
                
                log.write(OG+"\trerun\t")
                
                ntr=1
                while "toremove_"+OG+"_"+str(ntr)+".txt" in os.listdir(d):
                    ntr+=1
                with open(d+"toremove_"+OG+"_"+str(ntr)+".txt","w") as OUT:
                    for itr in toremove:
                        OUT.write(itr+"\n")
                
                log.write("written\n")
                
                tolaunch.append(OG)
            
            rankog+=1    # we increase rankog as we change of line
    
    log.write("\n** analyse of the TreeShrink output done\n")
    
    
    ##### finally we write the jobs for the re-run and we re-run the pipeline A for the first 400 orthogroups
    
    dicoja={1:[]}  # dicoja is a dictionnary with the OG of tolaunch by group of 400
    
    nja=1    # nja is the number of the serie of 400 orthogroups
    
    nlp=0    # nlp is the number of orthogroups already in the current serie
    
    for OG in tolaunch:      # here we implement dicoja
        nlp+=1
        if nlp>400:
            nja+=1
            dicoja[nja]=[]
            nlp=1
        dicoja[nja].append(OG)
    
    
    
    # nserie is the number of the serie of job_array we write for the pipeline A
    # (ex. nserie = 2 if it is the second time we run the pipeline A)
    nserie=1
    while "job_array_pA_"+str(nserie)+"-1.sh" in os.listdir(d):
        nserie+=1
    
    # here we write the job arrays to launch
    for nja in dicoja:
        with open(d+"job_array_pA_"+str(nserie)+"-"+str(nja)+".sh","w") as OUT:
            OUT.write("#!/bin/bash\n\n")
            for OG in dicoja[nja]:
                OUT.write("sbatch "+d+"pipe_A_*_"+OG+".sh\n")
    
    log.write("** job arrays written\n")
    
    # finally we re-run the pipelina A for the first 400 orthogroups
    if len(tolaunch)>0:
        subprocess.check_call(["sbatch","-p","long",d+"job_array_pA_"+str(nserie)+"-1.sh"])
    
    log.write("** first job array launched (400 OGs)\n")
    
    if len(dicoja)>1:
        log.write("\n!!!! warning !!!!\n "+str(len(dicoja)-1)+
                  " other job array(s) to launch manually later\n")
    
    log.write("\nEND")

























