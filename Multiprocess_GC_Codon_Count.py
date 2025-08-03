# -*- coding: utf-8 -*-
"""
Created on Tue May 13 21:59:57 2025

@author: okmaz
"""

import os 
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool 

fileloc="D:/Filtered_Bacteria_CDS/Filtered_Actinomycetota_CDS"
genome_cDNA_list=sorted(os.listdir(fileloc))

def Codon_GC_Definer(genome):
    with open(os.path.join(fileloc,genome), "r") as input_file: 
        genedict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
        ####produce library with keys being codons
        gccontent=[0]*500
        ####produce a list the elngth required 
        lengthlist=[0]*500
    for gene in genedict.values():
        for p in range(min(500, len(gene.seq) // 3)):
            GCcount= gene.seq[p*3:((p*3)+3)].count("G")+gene.seq[p*3:((p*3)+3)].count("C")
            gccontent[p]+=GCcount/3
            lengthlist[p]+=1
    
    #####Calculate average GC per position and add to grand list 
    genname1=os.path.splitext(genome)[0]
    genname2=genname1.replace("filtered_", "")
    relgc=[]
    for n in range(len(lengthlist)):
        relpos=gccontent[n]/lengthlist[n] if lengthlist[n] else 0
        relgc.append(relpos)
        #####add relgc into a list 
    
    relgc.insert(0,genname2)
    return relgc
        

specieslist=[]
if __name__ == "__main__":
    with Pool(processes=os.cpu_count()) as pool: 
         specieslist= pool.map(Codon_GC_Definer, genome_cDNA_list)

colnames=["assembly_accession"]+list(range(1,501))
metagenedf=pd.DataFrame(specieslist, columns=colnames)
group=os.path.basename(fileloc)
group2=group.replace("Filtered_", "").replace("_CDS", "")
metagenedf.to_csv(f"D:/Bacteria_Summative_GC/{group2}_all_GC.csv", index=False)
