# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script generates root pathway files and number codes for constructing DTox hierarchy 


## Functions 
source("src/functions.R");


## 0. Input arguments
input_root_file	<- "downloads/reactome/root_pathways.txt";	# input file that contains root pathway names and ids 
output_folder	<- "data/reactome/root/";			# output folder that stores processed root files for building DTox hierarchy 
out_map_file	<- "data/reactome/root_file_map.tsv";		# output file that provides mapping between root pathway names and number code for DTox hierarchy

## 1. Generate all possible combinations of root pathways, write to files that named by number code for DTox hierarchy 
# read in input root pathay file  
input_root_df <- read.delim(file = input_root_file, header = T, sep = "\t");
# generate all possible combinations of root pathways, as well as a number code for each root combination 
root_combo <- generate.all.element.combinations(input_root_df$root_id);
names(root_combo) <- sapply(names(root_combo), function(nrc) paste("rt_", nrc, sep = ""));
# write processed root pathway IDs to output folder, name files after number code 
root_combo_write <- mapply(function(rc, nrc){
	out_wd <- paste(output_folder, nrc, ".txt", sep = "");
	writeLines(rc, out_wd);
	return(1);
}, root_combo, names(root_combo));

## 2. Generate mapping between root pathway names and number code for DTox hierarchy 
# obtain IDs of each root pathways
root_names <- input_root_df$root_name;
names(root_names) <- input_root_df$root_id;
# obtain names of each root pathway combination 
root_combo_names <- sapply(root_combo, function(rc) paste(root_names[rc], collapse = ","));
# obtain IDs of each root pathway combination
root_combo_ids <- sapply(root_combo, function(rc) paste(rc, collapse = ","));
# build mapping data frame and output 
root_map_df <- data.frame(names(root_combo), root_combo_ids, root_combo_names);
colnames(root_map_df) <- c("file_id", "root_ids", "root_names");
write.table(root_map_df, file = out_map_file, sep = "\t", row.names = F, col.names = T, quote = F); 


