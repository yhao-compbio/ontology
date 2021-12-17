# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script processes raw ontology data to construct DTox hierarchy 


# Functions 
source("src/functions.R");


## 0. Input arguments
Args			<- commandArgs(T);
pathway_source		<- Args[1];		# source of pathways 
root_pathway_file	<- Args[2];		# input file that contains processed root pathways 
pathway2pathway_file	<- Args[3];		# input file that contains parent/child relationships between pathways  
gene2pathway_file	<- Args[4];		# input file that contains gene annotations of pathways 
output_file		<- Args[5];		# output folder 
min_pathway_size	<- as.numeric(Args[6]);	# minimal size of pathway to be included in hierarchy  
gene2reaction_file	<- Args[7];		# optional input file that contains gene annotations of chemical reactions (only applies when source of pathway is 'Reactome')  
structure2target_file	<- Args[8];		# optional input file that contains mapping between structure features to target genes derived from faeture selection (only applies when input feature type is structural)
feature_file		<- Args[9];		# input file that contains features of DTox model (input layer of hierarchy) 

## 1. Process pathway annotation and relation datasets for DTox hierarchy construction 
# read in root pathways, pathway relation file, and pathway gene annotation file  
root_pathway <- readLines(root_pathway_file);
pathway_relation <- read.delim(file = pathway2pathway_file, sep = "\t", header = F);
pathway_genes <- read.delim(file = gene2pathway_file, sep = "\t", header = F);
# remove duplicate gene annotations   
if(pathway_source == "Reactome"){
	human_id <- which(pathway_genes$V6 %in% "Homo sapiens")
	pathway_genes <- pathway_genes[human_id, c(1,2)]
	pathway_genes$V1 <- sapply(pathway_genes$V1, function(pgv) strsplit(pgv, "-")[[1]][[1]])
	pathway_genes <- unique(pathway_genes)
}
# remove pathways that do not pass minimal size threhold, and are not descendants of root pathways    
filter_results1 <- filter.pathways(root_pathway, pathway_relation, pathway_genes, min_pathway_size);
pathway_relation1 <- filter_results1$relation;
annotation_mat1 <- filter_results1$annotation; 

## 2. Process reaction annotation data for DTox hierarchy construction  
if(gene2reaction_file != "NA"){
	# read in reaction gene annotatin file  
	reaction_genes <- read.delim(file = gene2reaction_file, sep = "\t", header = F)
	# only keep reactions in human 
	human_id <- which(reaction_genes$V6 %in% "Homo sapiens")
	reaction_genes <- reaction_genes[human_id, c(1,2)]
	# remove duplicate gene annotations  
	reaction_genes$V1 <- sapply(reaction_genes$V1, function(rgv) strsplit(rgv, "-")[[1]][[1]])
	reaction_genes <- unique(reaction_genes)	
	# map reactions to lowest-level pathways, then remove pathways that are not ancestors of mapped lowest-level pathways
	filter_results2 <- filter.reactions(reaction_genes, pathway_relation1, annotation_mat1)
	lowest_level2 <- filter_results2$lowest
	pathway_relation2 <- filter_results2$relation
	annotation_mat2 <- filter_results2$annotation
	# add specified parameters to output file name 
	output_file <- paste(output_file, "ps", min_pathway_size, "re_1", sep = "_")
}
if(gene2reaction_file == "NA"){
	# obtain pathways at the lowest level of hierarchy (no children pathways) 
	lowest_level2 <- unique(setdiff(pathway_relation1[,2], pathway_relation1[,1]))	
	pathway_relation2 <- pathway_relation1
	annotation_mat2 <- annotation_mat1
	# add specified parameters to output file name 
	output_file <- paste(output_file, "ps", min_pathway_size, "re_0", sep = "_")
}

## 3. Process structure feature connections for DTox hierarchy construction  
# read in input features 
input_features <- readLines(feature_file);
# build a list in which each element contains genes/reactions/pathways in a layer of DTox hierarchy 
knowledge_hierarchy <- NULL;
knowledge_hierarchy[[1]] <- input_features;
if(structure2target_file != "NA"){
	# read in relation matrix bewteen structure features and target genes of interest 
	structure_target <- read.delim(file = structure2target_file, sep = "\t", header = F)
	# filter pathways by target genes of interest 
	filter_results3 <- filter.targets(structure_target[,1], lowest_level2, pathway_relation2, annotation_mat2)
	hidden_target3 <- filter_results3$inter_target 
	knowledge_hierarchy[[2]] <- hidden_target3
	lowest_level3 <- filter_results3$inter_lowest 	
	pathway_relation3 <- filter_results3$relation
	annotation_mat3 <- filter_results3$annotation
	# builds connections between target genes of interest, lowest-level pathways, and higher level pathways 
	knowledge_hierarchy_mat <- get.knowledge.hierarchy(hidden_target3, lowest_level3, annotation_mat3, pathway_relation3)
	# add connections between structure features and target genes of interest, in order to form DTox hierarchy
	ht_id <- which(structure_target[,1] %in% hidden_target3)
	knowledge_hierarchy_mat <- rbind(structure_target[ht_id,], knowledge_hierarchy_mat)
	# add specified parameters to output file name 
	output_file <- paste(output_file, "st_1", sep = "_")
}
if(structure2target_file == "NA"){
	# filter pathways by target genes of input features 
	filter_results3 <- filter.targets(input_features, lowest_level2, pathway_relation2, annotation_mat2)
	hidden_target3 <- filter_results3$inter_target
	lowest_level3 <- filter_results3$inter_lowest
	pathway_relation3 <- filter_results3$relation
	annotation_mat3 <- filter_results3$annotation
	# builds connections between target genes of interest, lowest-level pathways, and higher level pathways to form DTox hierarchy 
	knowledge_hierarchy_mat <- get.knowledge.hierarchy(hidden_target3, lowest_level3, annotation_mat3, pathway_relation3)
	# add specified parameters to output file name 
	output_file <- paste(output_file, "st_0", sep = "_")
}

## 4. Construct and sort DTox hierarchy using the processed datasets 
# use child/parent relationships to sort pathways in the DTox hierarchy by layer
sort_results <- sort.nodes.by.hierarchy(lowest_level3, pathway_relation3, knowledge_hierarchy);
# convert entities in knowledge matrix to their corresponding node numbers in DTox hierarchy  
convert_results <- convert.node.hierarchy.to.numbers(knowledge_hierarchy_mat, sort_results$layer, sort_results$node);
# obtain pathway size and root pathway name information 
info_results <- get.node.information(sort_results$node, annotation_mat3, root_pathway);

## 5. Write DTox hierarchy files to output files 
# node layer file 
write.table(sort_results$layer, file = paste(output_file, "_layer.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T);
# node number file 
write.table(sort_results$node, file = paste(output_file, "_node.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T);
# node connection file by structure/gene/pathway node name
write.table(convert_results$origin, file = paste(output_file, "_knowledge_by_name.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T);
# node connection file by node number 
write.table(convert_results$number, file = paste(output_file, "_knowledge_by_node.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T);
# pathway node size file  
write.table(info_results$size, file = paste(output_file, "_node_size.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T);
# root node name file 
write.table(info_results$root, file = paste(output_file, "_root.tsv", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T);
