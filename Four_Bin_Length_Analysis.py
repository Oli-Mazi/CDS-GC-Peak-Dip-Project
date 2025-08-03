

####Extract data 
import os 
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool 

fileloc="D:/Filtered_Bacteria_CDS/Filtered_Thermotogati_CDS"
genome_cDNA_list=sorted(os.listdir(fileloc))


####Process each length set
def Codon_GC_Definer(genome):
    with open(os.path.join(fileloc,genome), "r") as input_file: 
        genedict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
        ####produce library with keys being codons
        Lgccontent=[0]*200
        Llengthlist=[0]*200
        
        Mgccontent=[0]*500
        Mlengthlist=[0]*500
        
        UPMgccontent=[0]*500
        UPMlengthlist=[0]*500
        
        Hgccontent=[0]*500
        Hlengthlist=[0]*500
        
    lowdict= {key: value for key, value in genedict.items() if len(value) <=600}
    middict= {key: value for key, value in genedict.items() if 600<len(value) <=1500 }
    upmiddict= {key: value for key, value in genedict.items() if 1500<len(value) <=2400 }   
    highdict= {key: value for key, value in genedict.items() if len(value) >2400 }


    for gene in lowdict.values():
        for p in range(min(200, len(gene.seq) // 3)):
            GCcount= gene.seq[p*3:((p*3)+3)].count("G")+gene.seq[p*3:((p*3)+3)].count("C")
            Lgccontent[p]+=GCcount/3
            Llengthlist[p]+=1

    for gene in middict.values():
        for p in range(min(500, len(gene.seq) // 3)):
            GCcount= gene.seq[p*3:((p*3)+3)].count("G")+gene.seq[p*3:((p*3)+3)].count("C")
            Mgccontent[p]+=GCcount/3
            Mlengthlist[p]+=1
            
    for gene in upmiddict.values():
        for p in range(min(500, len(gene.seq) // 3)):
            GCcount= gene.seq[p*3:((p*3)+3)].count("G")+gene.seq[p*3:((p*3)+3)].count("C")
            UPMgccontent[p]+=GCcount/3
            UPMlengthlist[p]+=1

    for gene in highdict.values():
        for p in range(min(500, len(gene.seq) // 3)):
            GCcount= gene.seq[p*3:((p*3)+3)].count("G")+gene.seq[p*3:((p*3)+3)].count("C")
            Hgccontent[p]+=GCcount/3
            Hlengthlist[p]+=1
            
    #####Calculate average GC per position and add to grand list 
    genname1=os.path.splitext(genome)[0]
    genname2=genname1.replace("filtered_", "")
    Lrelgc=[]
    for n in range(len(Llengthlist)):
        lrelpos=Lgccontent[n]/Llengthlist[n] if Llengthlist[n] else 0
        Lrelgc.append(lrelpos)
    
    Mrelgc=[]
    for n in range(len(Mlengthlist)):
        mrelpos=Mgccontent[n]/Mlengthlist[n] if Mlengthlist[n] else 0
        Mrelgc.append(mrelpos)

    UPMrelgc=[]
    for n in range(len(UPMlengthlist)):
        upmrelpos=UPMgccontent[n]/UPMlengthlist[n] if UPMlengthlist[n] else 0
        UPMrelgc.append(upmrelpos)
        
    Hrelgc=[]
    for n in range(len(Hlengthlist)):
        hrelpos=Hgccontent[n]/Hlengthlist[n] if Hlengthlist[n] else 0
        Hrelgc.append(hrelpos)
    
    Lrelgc.insert(0,genname2)
    Mrelgc.insert(0,genname2)
    UPMrelgc.insert(0,genname2)
    Hrelgc.insert(0,genname2)
    
    return Lrelgc, Mrelgc,UPMrelgc, Hrelgc
        

Lspecieslist=[]
Mspecieslist=[]
UPMspecieslist=[]
Hspecieslist=[]
if __name__ == "__main__":
    with Pool(processes=os.cpu_count()) as pool: 
         results= pool.map(Codon_GC_Definer, genome_cDNA_list)
         for low, mid,upmid, high in results:
             Lspecieslist.append(low)
             Mspecieslist.append(mid)
             UPMspecieslist.append(upmid) 
             Hspecieslist.append(high)

lowcolnames=["assembly_accession"]+list(range(1,201))
colnames=["assembly_accession"]+list(range(1,501))

Lmetagenedf=pd.DataFrame(Lspecieslist, columns=lowcolnames)
Mmetagenedf=pd.DataFrame(Mspecieslist, columns=colnames)
UPMmetagenedf=pd.DataFrame(UPMspecieslist, columns=colnames)
Hmetagenedf=pd.DataFrame(Hspecieslist, columns=colnames)

group=os.path.basename(fileloc)
group2=group.replace("Filtered_", "").replace("_CDS", "")

Lmetagenedf.to_csv(f"D:/Bacteria_Length_GC/{group2}_short_GC.csv", index=False)
Mmetagenedf.to_csv(f"D:/Bacteria_Length_GC/{group2}_medium_GC.csv", index=False)
UPMmetagenedf.to_csv(f"D:/Bacteria_Length_GC/{group2}_mediumupper_GC.csv", index=False)
Hmetagenedf.to_csv(f"D:/Bacteria_Length_GC/{group2}_long_GC.csv", index=False)