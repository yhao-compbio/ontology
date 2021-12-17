# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script shuffles the hierarchy of input DTox model while preserving the distribution of connections between node layers 


# Functions 
source("src/functions.R");


## 0. Input arguments
Args		<- commandArgs(T);
knowledge_file	<- Args[1];	# input file that contains 
layer_file	<- Args[2];	# input file that contains 
out_file	<- Args[3];	# output file of shuffled connections of DTox hierarchy

## 1. Obtain layer information of DTox 
# read in layer number matrix, build vector that maps node number and layer number
layer_df <- read.delim(file = layer_file, header = T, sep = "\t"); 
node_layer <- as.character(layer_df$layer_number);
names(node_layer) <- layer_df$node;
# obtain the nodes of each layer 
layer_node_list <- group.vector.by.categories(layer_df$layer_number, layer_df$node);

## 2. Shuffle the DTox hierarchy while preserving the distribution of connections between node layers
# read in knowledge matrix 
knowledge_df <- read.delim(file = knowledge_file, header = T, sep = "\t");
# set seed number 
set.seed(0);
# iterate by each parent node from layer 1+ (those with connections to children nodes)  
shuffle_child <- sapply(knowledge_df$children_node, function(kdcn){
	# obtain children nodes of current node 
	kdcn_child <- strsplit(kdcn, ",")[[1]];
	# obtain layer summary of the children nodes 
	kdcn_layer <- node_layer[kdcn_child];
	kdcn_table <- table(kdcn_layer);
	# iterate by the distinct layer of the children nodes 
	kdcn_sample <- mapply(function(nkt, kt){
		# sample the same number of nodes from the current layer  
		nkt_sample <- sample(layer_node_list[[nkt]], kt);
		nkt_sample <- paste(sort(nkt_sample), collapse = ",");
		return(nkt_sample); 
	}, names(kdcn_table), kdcn_table);
	# aggregate sampled nodes from different layers 
	kdcn_sample <- paste(kdcn_sample, collapse = ",");
	return(kdcn_sample);	
});
# store shuffled children nodes of each parent node in data frame form, output 
shuffle_df <- knowledge_df;
shuffle_df$children_node <- shuffle_child;
write.table(shuffle_df, file = out_file, sep = "\t", row.names = F, col.names = T, quote = F);
