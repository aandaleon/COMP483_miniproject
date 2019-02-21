# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:55:25 2019

@author: Angela
"""

import Bio #well I can't install Biopython so
import numpy as np
import os
#import pandas as pd

#make directory and change into it
os.system("mkdir Angela_Andaleon")
os.chdir("Angela_Andaleon")

#get reference numbers for strains
strain_dict = {'HM27':'APNU00000000', 'HM46':'APNY00000000', 'HM65':'APNX00000000', 'HM69':'APNV00000000'}

#retrieve assemblies and count contigs & bp
UPEC_log = open("UPEC.log", "w")
for strain in strain_dict:
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AP/" + strain_dict[strain][2:4] + "/" + strain_dict[strain][0:4] + "01/" + strain_dict[strain][0:4] + "01.1.fsa_nt.gz") #retrieve contigs in fasta format
    os.system("gunzip " + strain_dict[strain][0:4] + "01.1.fsa_nt.gz") #idk if gzip is installed on here
    fastas = open(strain_dict[strain][0:4] + "01.1.fsa_nt").read().splitlines() #because np.loadtxt doesn't feel like working
    contigs = 0 #intitalize counts
    bp = 0
    for line in fastas: #note: I did this before Nick installed biopython
        if line.startswith(">"): #if start of new read
            contigs += 1
        else: #str of base pairs
            bp += len(line) #contigs < 1000 bp aren't allowed in SRA
    UPEC_log.write("There are " + str(contigs) + " contigs in the assembly " + strain + ".\nThere are " + str(bp) + " bp in the assembly " + strain + ".\n") #print to log

