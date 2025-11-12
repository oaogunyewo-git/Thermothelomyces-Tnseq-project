#####  pipeline for collapsing low-support branches in gene trees and launch astral
#
# this script is aimed to take a set of gene trees obtained after pipeline A and checked by pipeline B
# if a gene trees is not to be modified after TreeShrink (pipeline B),
# this script takes it directly from the directory fulfilled by the pipeline B,
# otherwise the gene tree is taken directly from the pipeline A output directory (raxml/OG/OG.raxml.support)
# then this script roots the gene trees and collapse low-support branches using the newick utilities
# finally, it launches astral on the whole gene trees dataset
#
# by Lucas Bonometti
#
#
#####  whats does this script do ?
#
# 1) 
#
#
##### points of warning
#
# 
#
#
#
#
###### libraries

import os
import sys
import subprocess
import argparse


##### creation of an argument parser
#
# those are the arguments needed to the script

parser = argparse.ArgumentParser(description="script filtering a single-copy orthogroup and "+
                                 "launch the cleaned alignment on raxml-ng")

parser.add_argument("-d","--work_dir",
                    help="# the path to the working directory (default : ./ )",
                    required=False, default="./")

parser.add_argument("-a","--astral",
                    help="# the directory where is astral",
                    required=True)

parser.add_argument("-l","--limit",
                    help="# the minimal bootstrap value under which we collapse (default : 10)",
                    required=False,default=10,type=int)

parser.add_argument("-o","--outgroup",
                    help="# the outgroup for the final tree (default : lolmin-1)",
                    required=False,default="lolmin-1")

args = parser.parse_args()



##### there begins the script


limit=args.limit        # limit is the limit for defining low-support branches

d=args.work_dir            # "d" is the path to the working directory
if d[-1]!="/": d=d+"/"

contend = []             # contend (with a "d") is the content in directories of the working directory
for item in os.listdir(d):
    if os.path.isdir(d+item):
        contend.append(item)


nlog=1                # nlog is the number of the log file
while "log_pipeC_"+str(nlog)+".log" in os.listdir(d):
    nlog+=1

with open(d+"log_pipeC_"+str(nlog)+".log","w") as log :
    
    # log is the log file, in which we will write allong the script
    
    log.write("** initialization done\n")
    
    drax=d+"raxmlng/"      # drax is the directory containing all the trees obtained with the pipeline A
    ds2c=d+"shrink2collapse/"    # "ds2c" is the directory containing the trees to input in the pipeline C
    
    dc2a=d+"collapse2astral"    # "dc2a" is the directory containing the trees for astral
    if "collapse2astral" not in contend:
        sp = subprocess.call(["mkdir",dc2a])
    elif len(os.listdir(dc2a))>0:
        sp = subprocess.Popen("rm "+dc2a+"/*",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    dc2a=dc2a+"/"
    
    doa=d+"output_astral"    # "doa" is the directory containing outputs of astral
    if "output_astral" not in contend:
        sp = subprocess.call(["mkdir",doa])
    elif len(os.listdir(doa))>0:
        sp = subprocess.Popen("rm "+doa+"/*",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    doa=doa+"/"
    
    das=args.astral                 # "das" is the directory in which there is astral
    if das[-1]!="/": das=das+"/"
    
    
    log.write("** environment building done\n")
    
    ##### now we have built the environmenet of the script
    
    listog1=[]      # listog1 is the list of orthogroups not re-run on pipeline A after pipeline B
    
    for tree in os.listdir(ds2c):
        if tree[-6:]=="ked.nw":
            listog1.append(tree.replace("_skrinked.nw",""))
        
    
    listog2=[]      # listog2 is the list of orthogroups re-run on pipeline A after pipeline B
    
    subdrax=os.listdir(drax)   # subdrax is a list of the orthogorups with a directory in drax
                               # (the gene trees are in those sub-directories)
    
    
    # here we copy all the re-run gene trees in shrink2collapse
    for OG in subdrax:
        if OG not in listog1:
            if OG+".raxml.support" in os.listdir(drax+OG):
                listog2.append(OG)
                sp = subprocess.check_call(["cp",drax+OG+"/"+OG+".raxml.support",ds2c+OG+"_skrinked.nw"])
    
    
    
    listog=listog1+listog2    # listog is the list of all orthogroups on which astral will be run
    
    log.write("** all the gene trees are in shrink2collapse\n")
    
    ##### this step is for rooting the gene trees and collapse the low-support branches
    
    for OG in listog:
        
        # first we root avery gene tree
        sp = subprocess.Popen(d+"nw_reroot "+ds2c+OG+"_skrinked.nw > "+ds2c+OG+".nw",
                              shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        
        if list(sp.communicate())[1].decode()!='':
            log.write("\n!!!!! error : error as we use nw_reroot on "+OG+" \n")
            log.write(list(sp.communicate())[1].decode())
            sys.exit()
        
        
        
        # and finally we collapse low-support branches and we put the obtained tree in collapse2astral
        sp = subprocess.Popen(d+"nw_ed "+ds2c+OG+".nw 'i & b<"+str(limit)+"' o > "+dc2a+OG+".nw",
                              shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        
        if list(sp.communicate())[1].decode()!='':
            log.write("\n!!!!! error : error as we use nw_ed on "+OG+" \n")
            log.write(list(sp.communicate())[1].decode())
            sys.exit()
    
    log.write("** rooting and collapsing done\n")
    
    if len(os.listdir(dc2a))!=len(listog):
        log.write("\nERROR : problem as we collapse the files in collapse2astral\n")
        sys.exit()
    
    
    ninp=1
    while "input"+str(ninp)+".trees" in os.listdir(das):
        ninp+=1
    
    sp = subprocess.Popen("cat "+dc2a+"* > "+das+"input"+str(ninp)+".trees",
                          shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    
    if list(sp.communicate())[1].decode()!='':
        log.write("\n!!!!! error : error as we concatene gene trees for astral\n")
        log.write(list(sp.communicate())[1].decode())
        sys.exit()
    
    sp = subprocess.Popen("module load picard/2.23.5 && java -jar "+das+"astral.5.7.5.jar -i "+
                          das+"input"+str(ninp)+".trees -o "+doa+"arbre_especes.nw --outgroup "+args.outgroup+
                          " 2> "+doa+"logastral"+str(ninp)+".log",
                          shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    
    if list(sp.communicate())[1].decode()!='':
        log.write("\n!!!!! error : error when using astral\n")
        log.write(list(sp.communicate())[1].decode())
        sys.exit()
    
    log.write("** astral done\n\n")
    log.write("YOU WIN")



































