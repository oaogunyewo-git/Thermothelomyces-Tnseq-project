#####  pipeline for filtering orthogroups alignments and create a gene tree with RAxML-NG
#
# this script is aimed to take a single-copy orthogroup (at most 1 sequence per individual),
# to filter the alignment using TranslatorX, and finally to launch RAxML-NG on this alignmnet
#
# by Lucas Bonometti
#
#
#####  which kinds of filtering ?
#
# 1) it removes the positions sustained by less than 50% of the individuals (through TranslatorX)
#
# 2) it removes "fragmented" sequences
#    i.e. sequences with more than k % of empty positions (- or N) after TranslatorX
#
# 3) it renames sequences with the name of the indivudual only
#    ex. >ind_gene_transcript   -->   >ind
#
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
from Bio import SeqIO
import subprocess
import argparse


##### creation of an argument parser
#
# those are the arguments needed to the script

parser = argparse.ArgumentParser(description="script filtering a single-copy orthogroup and "+
                                 "launch the cleaned alignment on raxml-ng")

parser.add_argument("-f","--fasta_sco",
                    help="# the path to the original orthogroup fasta file",
                    required=True)

parser.add_argument("-d","--work_dir",
                    help="# the path to the working directory (default : ./ )",
                    required=False, default="./")

parser.add_argument("-m","--minlength",
                    help="# the minimal length of the cleanali to conduct raxml-ng (default : 100)",
                    required=False,default=100,type=int)

parser.add_argument("-l","--lvl_frag",
                    help="# the minimal proportion of positions with a nucleotide in the cleanali "+
                    "to considere a sequence as non-fragmented (default : 0.40 )",
                    required=False,default=0.4,type=float)

parser.add_argument("-i","--individuals",
                    help="# file of all the individuals potentially present",
                    required=True)

args = parser.parse_args()



##### the functions used latter on the script



# 1. translatorx

# this function is used to launch translatorx as a subprocess
# it takes a fasta of sequences as entry and creates the outputs in directory
# it then removes everything but the interesting part (nt_cleanali)
# it returns 0,0 if well functionnated

def translatorx(fasta,directory,OG,wd):
    sptr = subprocess.Popen([wd+"translatorX.sif",
                            "-i",fasta,
                            "-o",directory+OG,
                            "-p","F",
                            "-g","'-b5=half'"],
                           stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    
    sptr = sptr.communicate()
    sp1,sp2=list(sptr)[0].decode(),list(sptr)[1].decode()
    
    sprm = subprocess.Popen("rm "+
                           directory+OG+".aa* "+
                           directory+OG+".html "+
                           directory+OG+".mafft.log "+
                           directory+OG+".nt1* "+
                           directory+OG+".nt2* "+
                           directory+OG+".nt3* "+
                           directory+OG+".nt_ali.fasta",
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    sp3,sp4=list(sprm.communicate())[0].decode(),list(sprm.communicate())[1].decode()
    
    return sp1,sp2,sp3,sp4



# 2. check_raxml

# this function is used to evaluate the cleanali obtained
# it takes as entry the fasta of the cleaali
# it first launches "--check" in raxml-ng
# and next launches "--parse" in raxml-ng
# if it detectes an ERROR in the cleanali file, it stops the program
# othermise it will create a binary file (.raxml.rba)
# and will return the estimated memory required and the recommanded number of cpus

def check_raxml(fasta,log,wd,OG):
    sp = subprocess.check_output("module load raxml-ng/0.9.0 && raxml-ng --check --msa "+
                                 fasta+" --model GTR+G",shell=True)    
    sp = str(sp).split("\\n")
    warning = [match for match in sp if "WARNING" in match]
    error = [match for match in sp if "ERROR" in match]
    note = [match for match in sp if "NOTE" in match]
    
    if len(error)>0:
        log.write("\nerror in the final fasta : it can't be used by RAxML-NG\n")
        for i in error:
            log.write(i+"\n")
        with open(wd+"ERROR_"+OG+".txt","w") as erfile:
            erfile.write("error in the final fasta : it can't be used by RAxML-NG\n")
            for i in error:
                erfile.write(i+"\n")
        sys.exit()
    
    if len(warning)>0:
        log.write("\nsome warnings have been provided\n")
        for i in warning:
            log.write(i+"\n")
        log.write("\n")
    
    if len(note)>0:
        log.write("\na reduced alignment has been created by raxml-ng\n")
        sp = subprocess.check_call(["rm",fasta+".raxml.reduced.phy"])
        log.write("and has been deleted since\n")
    
    sp = subprocess.check_output("module load raxml-ng/0.9.0 && raxml-ng --parse --msa "+
                                 fasta+" --model GTR+G",shell=True)
    
    sp = str(sp).split("\\n")
    error = [match for match in sp if "ERROR" in match]
    note = [match for match in sp if "NOTE: Reduced alignment" in match]
    memory = [match for match in sp if "Estimated memory requirements" in match]
    cpu = [match for match in sp if "Recommended number of threads / MPI processes" in match]
    
    if len(error)>0:
        log.write("\nerror in the final fasta : it can't be used by RAxML-NG\n")
        for i in error:
            log.write(i+"\n")
        with open(wd+"ERROR_"+OG+".txt","w") as erfile:
            erfile.write("error in the final fasta : it can't be used by RAxML-NG\n")
            for i in error:
                erfile.write(i+"\n")
        sys.exit()
    
    if len(note)>0:
        log.write("\na reduced alignment has been created by raxml-ng\n")
        sp = subprocess.check_call(["rm",fasta+".raxml.reduced.phy"])
        log.write("and has been deleted since\n")
    
    memory=memory[0].split(":")[-1].strip()
    cpu=cpu[0].split(":")[-1].strip()
    
    return memory,cpu


##### this part gives a list of the whole set of individuals

with open(args.individuals,'r') as indivi:
    setind=indivi.read().splitlines()



#####  this part is aimed to indicate on what we work (whiech fasta of which orthogroup)

OG=args.fasta_sco.split("/")[-1].split(".")[0]   # "OG" is the name of the orthogroup
fasta0=args.fasta_sco                            # "fasta0" is the name of the original fasta


##### this part is aimed to build the environment of the script

d=args.work_dir            # "d" is the path to the working directory
if d[-1]!="/": d=d+"/"

content = os.listdir(d)    # "content" is the list of elements in the working directory

n=1
while "log_"+OG+"_"+str(n)+".log" in content:
    n+=1

with open(d+"log_"+OG+"_"+str(n)+".log","w") as log :
    
    log.write("**initialisation complete\n")
    
    dr=d+"reduced_sco"         # "dr" is the directory containing the orthogroups from which 
    if "reduced_sco" not in content:
        sp = subprocess.call(["mkdir",dr])
    elif OG+"_reduced.fasta" in os.listdir(dr):
        sp = subprocess.call(["rm",dr+"/"+OG+"_reduced.fasta"])
    dr=dr+"/"
    
    dtrx=d+"translatorx"          # "dtrx" is the name of the translatorx directory for the orthogroup 
    if "translatorx" not in content:
        sp = subprocess.call(["mkdir",dtrx])
    dtrx=dtrx+"/"
    if OG not in os.listdir(dtrx):
        sp = subprocess.call(["mkdir",dtrx+OG])
    elif len(os.listdir(dtrx+OG))>0:
        sp = subprocess.Popen("rm "+dtrx+OG+"/*",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    dtrx=dtrx+OG+"/"
    
    drax=d+"raxmlng"
    if "raxmlng" not in content:
        sp = subprocess.call(["mkdir",drax])
    drax=drax+"/"
    if OG not in os.listdir(drax):
        sp = subprocess.call(["mkdir",drax+OG])
    elif len(os.listdir(drax+OG))>0:
        sp = subprocess.Popen("rm "+drax+OG+"/*",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    drax=drax+OG+"/"
    
    log.write("**environmenet built\n")
    
    ##### this part is aimed to generate a first reduced fasta based on the original fasta
    
    # we load the original fasta as a dictionnary "fa"
    fa=SeqIO.to_dict(SeqIO.parse(fasta0,"fasta"))
    
    # we change all the genes names into individuals names  (ex. john_g12_t3 -> john)
    for elet in fa:
        fa[elet].id=fa[elet].id.split("_")[0]
        fa[elet].name=fa[elet].id.split("_")[0]
        fa[elet].description=fa[elet].id.split("_")[0]
        fa[elet].seq=fa[elet].seq.upper()
    
    # we write a file containing the individuals not present in the original fasta 
    with open("completeness_"+OG+".txt","w") as comp:
        present=[]
        for gene in fa:
            present.append(gene.split("_")[0])
        for tind in setind:
            if tind not in present:
                comp.write(tind+"\n")
    
    # if theer is at least one toremove file in the working directory (from the 2nd pipeline)
    ntr=1
    fa2=fa.copy()
    while "toremove_"+OG+"_"+str(ntr)+".txt" in content:
        with open(d+"toremove_"+OG+"_"+str(ntr)+".txt",'r') as toremove:
            lr=toremove.read().splitlines()
        for rem in lr:
            for gene in fa:
                if gene.split("_")[0]==rem and gene in fa2:
                    fa2.pop(gene)
        ntr+=1
    fa=fa2.copy()
    # we remove from "fa" all the genes that have been identified as "to removed" in the 2n pipeline
    
    # if there is at least one fragmented file in the working directory (from the 2nd pipeline)
    fa2=fa.copy()
    nfr=1
    while "fragmented_"+OG+"_"+str(nfr)+".txt" in content:
        with open("fragmented_"+OG+"_"+str(nfr)+".txt",'r') as fragmented:
            lf=fragmented.read().splitlines()
        for frag in lf:
            for gene in fa:
                if gene.split("_")[0]==frag and gene in fa2:
                    fa2.pop(gene)
        nfr+=1
    fa=fa2.copy()
    # we remove from "fa" all the genes that have been identified as "to removed" in the 2n pipeline
    
    # we write the first reduced fasta from "fa"
    SeqIO.write(fa.values(),dr+OG+"_reduced.fasta","fasta")
    
    
    log.write("**first reduced fasta written\n")
    
    #####  this part is the main loop of the script that finally write a cleanali fasta file to put in raxml-ng
    
    turn = 0  # counts the number of loops done
    
    fasta_cleanali=dtrx+OG+".nt_cleanali.fasta"  # "fasta_cleanali" is the name of the cleanali file
    
    fragmented=[]    # "fragmented" is the list of the sequences to remove because too fragmented
    l0=0             # l0 is the length of "fragmented" at the begining of each loop
    l1=1             # l1 is the length of "fragmented" at the end of each loop
    
    while l0!=l1:                # we stop the loop when we no more remove individuals from the reduces fasta
        l0=len(fragmented)
        
        if len(os.listdir(dtrx))>0:
            sp = subprocess.Popen("rm "+dtrx+"*",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        
        # we launch translatorx on the reduced fasta  
        sp1,sp2,sp3,sp4 = translatorx(dr+OG+"_reduced.fasta",dtrx,OG,d)

        if sp2!='':
            if "Illegal division by zero" in sp1:
                log.write("\nwarning : translatorX wrote an empty sequence in the cleanali\n")
            else:
                log.write("\ntranslatorX returned an error\n")
            if OG+".nt_cleanali.fasta" not in os.listdir(dtrx):
                log.write("\ntranslatorX did not write the cleanali file\n")
                with open(d+"ERROR_"+OG+".txt","w") as erfile:
                    erfile.write("translatorX did not write the cleanali file")
                sys.exit()
        if sp3!='' or sp4!='':
            with open(d+"ERROR_"+OG+".txt","w") as erfile:
                erfile.write("removal of translatorx elements went wrong\n"+sp3+"\n"+sp4)
            log.write("\nremoval of translatorx elements went wrong\n"+sp3+"\n"+sp4)
            sys.exit()
        
        # "fa" is now the dictionnary of the cleanali fasta created by translatorx
        fa=SeqIO.to_dict(SeqIO.parse(fasta_cleanali,"fasta"))
        
        # we add to fragmented all the genes with more than k % of empty positions (- or N) on the alignment
        for gene in fa:
            seq=fa[gene].seq
            long=len(seq)
            empty=seq.count("-")+seq.count("N")
            if empty/long > 1-args.lvl_frag:
                fragmented.append(gene)
        
        l1=len(fragmented)
        
        
        if l0!=l1:
            
            # we remove from the reduced fasta all the fragmented gene not already removed from it
            fa=SeqIO.to_dict(SeqIO.parse(dr+OG+"_reduced.fasta","fasta"))
            for frag in fragmented:
                if frag in fa:
                    fa.pop(frag)
            
            if len(fa)==0:
                with open(d+"too-frag_cleanali_"+OG+".txt","w") as ts:
                        ts.write("removing the too fragmented genes in the cleanali of "+OG+
                                 "removed all the genes")
                log.write("\nremoving the too fragmented genes in the cleanali removed all the genes")
                sys.exit()
        
            
            sp = subprocess.call("rm "+dr+OG+"_reduced.fasta",shell=True)
            
            # we write the new reduced fasta from "fa"
            SeqIO.write(fa.values(),dr+OG+"_reduced.fasta","fasta")
        
        turn += 1
        log.write("**turn "+str(turn)+" of the loop done\n")
    
    log.write("**loop done\n")    
    
    if len(fragmented)>0:
        n=1
        while "fragmented_"+OG+"_"+str(n)+".txt" in content:
            n+=1
        with open("fragmented_"+OG+"_"+str(n)+".txt","w") as frag:
            frag.write("\n".join(fragmented))
            
    
    ##### this part is aimed to be sure that the final cleanali is long enough
    
    
    # "fa" is now the cleanali fasta
    fa=SeqIO.to_dict(SeqIO.parse(fasta_cleanali,"fasta"))
    long=len(fa[list(fa.keys())[0]].seq)   # "long" is the length of the cleaned alignment
    
    # if long is shoter than the minimal admited length, the orthogroup is removed from teh dataset
    if long < args.minlength:
        log.write("too short cleanali !!!")
        log.write("only "+str(long)+" nucleotides !!")
        with open(d+"too-short_cleanali_"+OG+".txt","w") as ts:
            ts.write("the cleanali of "+OG+"contains less than "+str(args.minlength)+
                     " nucleotides\nconsequently the orthogroup has been removed from the dataset")
        sys.exit()
    
    log.write("**cleanali length controle done\n")
    
    ##### we launch check_raxml on the clean_ali fasta  (memory is for now on not used)
    
    memory,cpu = check_raxml(fasta_cleanali,log,d,OG)
    
    log.write("**raxml-ng check and parse OK\n"+memory+" estimated\n"+cpu+" cpus recommanded\n")
    
    ##### we finally write the raxml-ng launcher and we launch it
    
    
    job=d+"job_raxml-ng_"+OG+".sh"
    
    with open(job,"w") as OUT:
        OUT.write("#!/bin/bash\n\n")
        OUT.write("module load raxml-ng/0.9.0\n")
        OUT.write("raxml-ng --all --threads "+cpu+" --msa "+fasta_cleanali+".raxml.rba --model GTR+G ")
        OUT.write("--bs-trees autoMRE{1000} --tree pars{25},rand{25} --seed 12345 --prefix "+drax+OG
                  +" --data-type DNA")
    
    subprocess.call(["sbatch","-p","long","--output","slurm-raxml"+OG+".out",
                     "--mem-per-cpu=1GB","--cpus-per-task="+cpu,job])
    
    log.write("**raxml-ng launched as recommanded\n\nEND")















