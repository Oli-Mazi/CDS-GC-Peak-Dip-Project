####Load libraries
library(rotl)
library(ape)

#####Load Dataset
metadata=read.csv(file.choose())

####Convert organism names to a dataframe containing standardised names and identifiers
identifiers=tnrs_match_names(metadata$organism_name)

####Remove any organisms which do not have an ott_id (aka an identifier); (This means certain organism from the dataset will not be in the phylogeny)

identifiers=identifiers[!is.na(identifiers$ott_id),]

#####As some organisms may not be in the tree, filter out the missing organisms (Again, certain organisms wont appear in the phylogeny)
identifiers$presence=is_in_tree(identifiers$ott_id)
identifiers=identifiers[identifiers$presence==TRUE,]

###Make the tree
tree=tol_induced_subtree(ott_ids=identifiers$ott_id)

####Tree modifies labels by adding the ott_id, thus requires cleaning (also important for analyses)
tree$tip.label <- gsub("_ott[0-9]+", "", tree$tip.label)

###Plot the tree 
plot(tree, direction="upwards")

####If branch length required (e.g. in PGLS analyses), below adds branch lengths, although all branch lengths will be equal from root to tip.
tree=compute.brlen(tree)


