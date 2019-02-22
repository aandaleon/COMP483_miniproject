# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:55:25 2019

@author: Angela
"""

#I silenced all the prints to screen that I didn't make because I thought they were annoying
import os

#make directory and change into it
os.system("mkdir Angela_Andaleon")
os.chdir("Angela_Andaleon")

#get reference numbers for strains (is there a way to automate this?)
strain_dict = {'HM27':('APNU00000000', 'SRR1278956'), 'HM46':('APNY00000000', 'SRR1278960'), 'HM65':('APNX00000000', 'SRR1283106'), 'HM69':('APNV00000000', 'SRR1278963')}

#retrieve assemblies and count contigs & bp
UPEC_log = open("UPEC.log", "w")
for strain in strain_dict:
    print("Beginning analyses on " + strain + ".")
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AP/" + strain_dict[strain][0][2:4] + "/" + strain_dict[strain][0][0:4] + "01/" + strain_dict[strain][0][0:4] + "01.1.fsa_nt.gz") #retrieve contigs in fasta format
    os.system("gunzip " + strain_dict[strain][0][0:4] + "01.1.fsa_nt.gz") #idk if gzip is installed on here
    fastas = open(strain_dict[strain][0][0:4] + "01.1.fsa_nt").read().splitlines() #because np.loadtxt doesn't feel like working
    contigs = 0 #intitalize counts
    bp = 0
    for line in fastas: #note: I did this before Nick installed biopython
        if line.startswith(">"): #if start of new read
            contigs += 1
        else: #str of base pairs
            bp += len(line) #contigs < 1000 bp aren't allowed in SRA
    UPEC_log.write("There are " + str(contigs) + " contigs in the assembly " + strain + ".\nThere are " + str(bp) + " bp in the assembly " + strain + ".\n") #print to log
    print("Reported number of contigs and bp.")
    
    #annotation w/ Prokka
    UPEC_log.write("Prokka command: 'prokka " + strain_dict[strain][0][0:4] + "01.1.fsa_nt --outdir prokka_" + strain + " --prefix " + strain + " --genus Coli'\n") #write prokka command to log file
    os.system("prokka " + strain_dict[strain][0][0:4] + "01.1.fsa_nt --outdir prokka_" + strain + " --prefix " + strain + " --genus Coli") #prokka is too noisy
    anno = open("prokka_" + strain + "/" + strain + ".txt").read().splitlines()
    for line in anno: #write anno results
        if line.startswith("CDS"):
            CDS = int(line.split(": ")[1]) #store for comparing to the reference genome
            UPEC_log.write("There are " + str(CDS) + " CDS in the assembly " + strain + ".\n")
        elif line.startswith("tmRNA"):
            tmRNA = int(line.split(": ")[1])
            UPEC_log.write("There are " + str(tmRNA) + " tmRNA in the assembly " + strain + ".\n")
        elif line.startswith("tRNA"):
            tRNA = int(line.split(": ")[1])
            UPEC_log.write("There are " + str(tRNA) + " tRNA in the assembly " + strain + ".\n")
    
    #write discrepancies of CDS & tRNA
    if CDS >= 4140 and tRNA >= 89: #some weird conditionals to get the "more than" and "less than" correct in the print statements
        UPEC_log.write("There are " + str(CDS - 4140) + " CDS more and " + str(tRNA - 89) + " more tRNA than the RefSeq in assembly " + strain + ".\n")
    elif CDS <= 4140 and tRNA <= 89:
        UPEC_log.write("There are " + str(4140 - CDS) + " CDS less and " + str(89 - tRNA) + " less tRNA than the RefSeq in assembly " + strain + ".\n")
    elif CDS >= 4140 and tRNA <= 89:
        UPEC_log.write("There are " + str(CDS - 4140) + " CDS more and " + str(89 - tRNA) + " less tRNA than the RefSeq in assembly " + strain + ".\n")
    elif CDS <= 4140 and tRNA >= 89:
        UPEC_log.write("There are " + str(4140 - CDS) + " CDS less and " + str(tRNA - 89) + " more tRNA than the RefSeq in assembly " + strain + ".\n")
    print("Wrote annotations and discrepancies from the reference genome.")
    
    #get transcriptome
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/" + strain_dict[strain][1][0:6] + "/" + strain_dict[strain][1] + "/" + strain_dict[strain][1] + ".sra > /dev/null 2>&1")
    os.system("fastq-dump -I --split-files " + strain_dict[strain][1] + ".sra") #get reads
    print("Retrived RNA-seq reads.")
    
    #build index, map reads and quantify expression
    os.system("bowtie2-build prokka_" + strain + "/" + strain + ".fna " + strain) #build index
    os.system("cp "+ strain_dict[strain][0][0:4] + "01.1.fsa_nt " + strain + ".fa")  #copy fasta into new name so tophat can find it
    os.system("tophat2 -o tophat_" + strain + " " + strain + " " + strain_dict[strain][1] + "_1.fastq  " + strain_dict[strain][1] + "_2.fastq") #map reads
        #quantify expression
    print("Finished writing analyses on strain " + strain + ".")

UPEC_log.close()
    
