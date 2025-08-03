####Read CSV file
metageno=read.csv(file.choose())

####Filter out Scaffolds and Contigs
metageno2=metageno[metageno$assembly_level %in% c("Complete Genome","Chromosome"),]

#####Genus division
library(dplyr)
library(stringr)

metageno3=metageno2 %>%mutate(genus = word(organism_name, 1))
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

check=as.data.frame(table(metageno6$genus))
check2=as.data.frame(table(check$Freq))
####Remove Genus and dsize
metageno6$dsize <- NULL
metageno6$genus <- NULL

####Save CSV
write.csv(metageno6,file='D:/Assembly_Metadata/Invert_onepergenus.csv', row.names=FALSE)
