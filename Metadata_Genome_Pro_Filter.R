####Read CSV file
metageno=read.csv(file.choose())

####Filter out Scaffolds and Contigs
metageno2=metageno[metageno$assembly_level %in% c("Complete Genome","Chromosome"),]

#####Genus division
library(dplyr)
library(stringr)

metageno3=metageno2 %>%mutate(genus = word(organism_name, 1))
metageno3$genus=str_replace_all(metageno3$genus, "\\[|\\]|'", "")


rm(metageno, metageno2)

####Filter out Chromosome-level assemblies if Complete-Level exists
genuslist=unique(metageno3$genus)
metageno4=metageno3[0,]
for (gen in genuslist){
  tempdata=metageno3[metageno3$genus==gen,]
  if (any(tempdata$assembly_level=="Complete Genome")){
    metageno4=rbind(metageno4,tempdata[tempdata$assembly_level=="Complete Genome",])
  } else{
    metageno4=rbind(metageno4,tempdata)
  }
}


####Select Variants with smaller Ungap v Gap difference 

metageno4$dsize=metageno4$genome_size-metageno4$genome_size_ungapped
genuslist=unique(metageno4$genus)
metageno5=metageno4[0,]

for (gentype in genuslist){
  tempdata=metageno4[metageno4$genus==gentype,]
  metageno5=rbind(metageno5, tempdata[tempdata$dsize== min(tempdata$dsize),])
}

rm(metageno3, gen, gentype, genuslist, tempdata, metageno4)

####Remove any further redundant genera via Gene number
genuslist=unique(metageno5$genus)
metageno6=metageno5[0,]

metageno6 <- metageno5 %>%group_by(genus) %>% filter(protein_coding_gene_count == max(protein_coding_gene_count))


####Select the one with the larger genome 
genuslist=unique(metageno6$genus)
metageno7=metageno6[0,]

for (gentype in genuslist){
  tempdata=metageno6[metageno6$genus==gentype,]
  metageno7=rbind(metageno7, tempdata[tempdata$genome_size== max(tempdata$genome_size),])
}

Achecking=as.data.frame(table(metageno7$genus))
Achecking2=as.data.frame(table(Achecking$Freq))

####Remove redundant 
metageno8 <- metageno7 %>%group_by(species_taxid) %>%slice_sample(n = 1) %>%ungroup()

Achecking=as.data.frame(table(metageno8$genus))
Achecking2=as.data.frame(table(Achecking$Freq))

####Remove Genus and dsize
metageno8$dsize <- NULL
metageno8$genus <- NULL

####Save CSV
write.csv(metageno8,file='D:/Assembly_Metadata/bacteria_onepergenus.csv', row.names=FALSE)
