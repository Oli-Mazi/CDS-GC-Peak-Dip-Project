
####Import libraries
from Bio import SeqIO
import os

###"direction" variable example "D:/Bacteria_Genome_CDS/Thermotogati_Genome_CDS/", folder containing CDS fasta files, with each folder named using assembly identifier
direction=

####Create list of all files (with their pathways)
genomelist=[]
for (root,_,files) in os.walk(direction, topdown=True):
    for file in files:
        if file[-4:]==".fna":
            genomelist.append(os.path.join(root,file))

####Multi-faceted Filter (cDNA to dictionary conversion)

for genome in genomelist:
    genedict={}
    ngenedict={}
    fgenedict={}
    with open(genome, "r") as input_file: 
        genedict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    
    for gene, seq in genedict.items():
        if len(seq.seq)%3==0 and seq.seq[0:3]=="ATG" and not any(nuc in seq.seq for nuc in {"Y", "R", "N"}):
            ngenedict[gene]=seq.seq
    genedict.clear()
    for gene, seq in ngenedict.items():
        if "*" not in seq.translate()[:-1] and seq.translate()[-1]=="*":
            fgenedict[gene]=seq
    
    foldname=os.path.basename(os.path.dirname(genome))
    with open(f"D:/Filtered_Bacteria_CDS/Filtered_Thermotogati_CDS/filtered_{foldname}.fna", "w") as file:
        for head, seq in fgenedict.items():
            file.write(f">{head}\n")
            file.write(f"{seq}\n")








