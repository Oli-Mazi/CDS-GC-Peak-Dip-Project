
from Bio import Entrez
import pandas as pd
import time

####Metadata pathway 
metadata=pd.read_csv("")

###Enter email
Entrez.email=""

def lineage_finder(txd):
    time.sleep(1)
    stream=Entrez.efetch(id=txd, db="taxonomy", retmode="xml")
    record=Entrez.read(stream)
    stream.close()
    lineage=record[0]["Lineage"]
    return lineage

metadata["Lineage"]=metadata["taxid"].apply(lineage_finder)


