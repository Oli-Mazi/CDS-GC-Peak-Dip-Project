####Load packages
library(biomartr)
library(stringr)
library(dplyr)
library(Biostrings)
library(parallel)
library(ggplot2)
library(ggsci)
library(caper)
library(ggtree)
library(rotl)
library(ape)
library(moments)

####Load pathway

metadata=read.csv(file.choose())

dropoff="D:/Bacteria_Genomes_Full/Thermotogati_Genomes_Full"


####Loop through list of accessions and download genomes

for (genome in metadata$assembly_accession){
  tryCatch({
  getGenome(db = "refseq", organism = genome, skip_bacteria =FALSE, path=file.path(dropoff, genome), gunzip=TRUE)
  }, error = function(e) {print(paste("skipped:", genome))})
  
}



####Create CSV containing all the data######

####Load list of all genome names####

genfiles =list.files("D:/Eukaryote_Genomes_Full/Fungi_Genome_Full", pattern = paste0("_genomic_refseq", ".*\\.fna$"), full.names = TRUE)



####Calc Length and GC content per genome chromosomes####

finaldata=data.frame(assembly_accession = character(),chr = character(),Length = numeric(),GC_content = numeric(),stringsAsFactors = FALSE)


for (pathway in genfiles) {
  accession=sub(".*(GCF_\\d+\\.\\d+).*", "\\1", pathway)
  tryCatch({
    genome=readDNAStringSet(pathway)
    
    genomedata=data.frame(
      chr = names(genome),
      Length = width(genome),
      GC_content = rowSums(letterFrequency(genome, c("G", "C"), as.prob = TRUE)),
      stringsAsFactors = FALSE
    )
    genomedata=genomedata %>%filter(!grepl("scaffold|NotPlaced|unplaced|plastid|contig|Unknown|chrUn|mitochondrion|plasmid|MIT|API|chloroplast|apicoplast|chromosome Z|chromosome W|chromosome Y", chr))
    rm(genome)
    gc()
    
    if (!any(grepl("chromosome V|chromosome I|chromosome II", genomedata$chr))) {
      genomedata <- genomedata %>% filter(!grepl("chromosome X", chr))
    }
    
  
    
    result=data.frame(assembly_accession = accession,chr = genomedata$chr, Length = genomedata$Length,GC_content = genomedata$GC_content,stringsAsFactors = FALSE)
    
    
    
    finaldata=bind_rows(finaldata, result)
    print(accession)
  }, error = function(e) {print(paste("skipped:", accession))})
}



####Write it ####

write.csv(finaldata, "D:/Eukaryote_Genomes_Full/Chr_GC_Data/fungi_chr_data.csv", row.names = FALSE)




####Linear Regression analysis####



metafiles =list.files("D:/Eukaryote_Genomes_Full/Chr_GC_Data", pattern = paste0("_chr_data", ".*\\.csv$"), full.names = TRUE)
chrdata= bind_rows(lapply(metafiles, read.csv))


chrdata$Length=chrdata$Length/1000000
chrdata$GC_content=chrdata$GC_content*100
assemblylist=unique(chrdata$assembly_accession)


linearoutcome=data.frame("assembly_accession"=as.character(), "intercept"=as.numeric(), "slope"=as.numeric(), "p_value"=as.numeric(), "chr_number"=as.numeric(), "p_division"=as.numeric(), "p_binary"=as.numeric(), "Chr_Len_mean"=as.numeric(), "Chr_Len_SD"=as.numeric(), "Chr_Len_Skew"=as.numeric())

for (assembly in assemblylist){
  species=chrdata %>% filter(grepl(assembly,assembly_accession ))
  
  if (nrow(species)>2){
    model=lm(GC_content~Length,data=species)
    tempdata=data.frame("assembly_accession"=assembly, "intercept"=model$coefficients[1], "slope"=model$coefficients[2], "p_value"=summary(model)$coefficients["Length","Pr(>|t|)"], "chr_number"=nrow(species), "p_division"=NA, "p_binary"=NA, "Chr_Len_mean"=mean(species$Length),"Chr_Len_SD"=sd(species$Length), "Chr_Len_Skew"=skewness(species$Length))
    
    
    linearoutcome=rbind(linearoutcome,tempdata)
  }
  

}




linearoutcome$p_division=ifelse(linearoutcome$p_value > 0.05,"P-value>0.05","P-value<=0.05")
linearoutcome$p_binary=ifelse(linearoutcome$p_value > 0.05,1,0)



##################Chromosome Size analysis

metafiles =list.files("D:/Eukaryote_Genomes_Full/Chr_GC_Data", pattern = paste0("_chr_data", ".*\\.csv$"), full.names = TRUE)
chrdata= bind_rows(lapply(metafiles, read.csv))


chrdata$Length=chrdata$Length/1000000
assemblylist=unique(chrdata$assembly_accession)


meansize=data.frame("assembly_accession"=as.character(),"Chr_Len_mean"=as.numeric())

for (assembly in assemblylist){
  species=chrdata %>% filter(grepl(assembly,assembly_accession ))
    tempdata=data.frame("assembly_accession"=assembly, "Chr_Len_mean"=mean(species$Length))
    
    
    meansize=rbind(meansize,tempdata)
}


metafiles =list.files("D:/Assembly_Metadata/Eukaryote_Metadata", pattern = paste0("_onepergenus", ".*\\.csv$"), full.names = TRUE)
data= bind_rows(lapply(metafiles, read.csv))

identifiers=tnrs_match_names(data$organism_name)
identifiers=identifiers[!is.na(identifiers$ott_id),]


identifiers$presence=is_in_tree(identifiers$ott_id)

identifiers=identifiers[identifiers$presence==TRUE,]


tree=tol_induced_subtree(ott_ids=identifiers$ott_id)
plot(tree, direction="upwards")
tree$tip.label <- gsub("_ott[0-9]+", "", tree$tip.label)

tree=compute.brlen(tree)  

gccordata=data.frame("Level"=character(),"GC_type"=character(),"Slope"=numeric(), "intercept"=numeric(), "p_value_slope"=numeric(), "p_value_intercept"=numeric())

metadata=merge(data, meansize, by="assembly_accession")


levlist=unique(metadata$Level)

gc_list=list("_all_GC", "_GC1", "_GC2", "_GC3", "_GC4")


for (studied_gc in gc_list){
  
  
  genefiles =list.files("D:/Eukaryote_Summative_GC", pattern = paste0(studied_gc, ".*\\.csv$"), full.names = TRUE)
  metagene= bind_rows(lapply(genefiles, read.csv))
  
  metagene$means=rowMeans(metagene[,2:41], na.rm = TRUE)-rowMeans(metagene[, 401:501], na.rm = TRUE)
  
  

  colnames(meansize) = c("assembly_accession", "chrmetric")
  
  metadata= merge(metadata, meansize, by = "assembly_accession")
  
  #
  
  means_data=data.frame(assembly_accession = metagene$assembly_accession, means = metagene$means)
  
  collective_data=merge(means_data,metadata , by = "assembly_accession", all = TRUE)
  collective_data$organism_name=tolower(collective_data$organism_name)
  species_names=data.frame(organism_name=identifiers$search_string, species=identifiers$unique_name)
  
  
  for (lev in levlist){
    collective_data_2nd=subset(collective_data, Level ==lev)
    
    prefinal_data=merge(species_names, collective_data_2nd, by = "organism_name", all = TRUE)
    final_data=data.frame(species=prefinal_data$species, codon40mean=prefinal_data$means, chrmetric=prefinal_data$chrmetric)
    
    final_data$species=gsub(" ", "_", final_data$species)
    
    
    final_data2=subset(final_data, !is.na(final_data$species))
    
    final_data2 <- final_data2[!duplicated(final_data2$species), ]
    
    compdata=comparative.data(tree, final_data2, species,vcv=TRUE, na.omit=TRUE, force.root = TRUE)
    
    model=pgls(codon40mean ~ chrmetric, data = compdata)
    summary(model)
    
    gcname=gsub("_all_|_", "", studied_gc)
    
    tempdata=data.frame("Level"=lev,"GC_type"=gcname,"Slope"=model$model$coef[2], "intercept"=model$model$coef[1], "p_value_slope"=summary(model)$coefficients["chrmetric","Pr(>|t|)"], "p_value_intercept"=summary(model)$coefficients["(Intercept)","Pr(>|t|)"])
    
    gccordata=rbind(gccordata,tempdata)
    

  }
  
  
}


write.csv(gccordata, "C:/Users/mean_size_gcdiff_indiv.csv", row.names = FALSE)


