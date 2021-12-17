# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script generates shell scripts that processes reactome pathway data to construct DTox hierarchy.


## functions 
source("src/functions.R");


## 0. Input arguments
pathway_source		<- "Reactome";
root_file_loc		<- "data/reactome/root/";
pathway2pathway_file	<- "downloads/reactome/ReactomePathwaysRelation.txt";
gene2pathway_file	<- "downloads/reactome/UniProt2Reactome_All_Levels.txt";
output_loc		<- "data/reactome/hierarchy/";
gene2reaction_file	<- "downloads/reactome/UniProt2ReactomeReactions.txt";
structure2target_file	<- "https://raw.githubusercontent.com/yhao-compbio/TTox/master/data/compound_target_0.25_binary_feature_select_implementation/fingerprint_maccs_analysis/fingerprint_maccs_select_features_mc_0.85_target_structure.tsv";
structure_file		<- "data/feature/maccs_fingerprints.txt";
target_file		<- "data/feature/target_profile_from_maccs_fingerprints.txt";
job_name		<- "get_knowledge_hierarchy_reactome"
N_cores			<- 4;

## 1. Obtain names of input and output files related to root pathways 
# list all input files of root pathways of interest
root_files <- list.files(root_file_loc);
root_input <- sapply(root_files, function(rf) paste(root_file_loc, rf, sep = ""));
# generate name for each output files of root pathways
root_name <- sapply(root_files, function(rf) strsplit(rf, ".txt")[[1]][[1]]); 
root_output <- sapply(root_name, function(rn) paste(output_loc, rn, sep = ""));

## 2. Generate parts of commands
part <- NULL;
# parts that include pathway source, root pathways of interest, pathway relation file, pathway annotation file, output root file 
part[[1]] <- mapply(function(ri, ro){
	rio <- paste("Rscript", "src/get_knowledge_hierarchy.R", pathway_source, ri, pathway2pathway_file, gene2pathway_file, ro, sep = " ");
	return(rio);
}, root_input, root_output);
# parts that include minimal size threshold of pathways  
part[[2]] <- c(5, 10, 20);
# parts that include reaction annotation file 
part[[3]] <- c(gene2reaction_file, "NA");
# parts that include structure-target connectoin file, and feature list file
part[[4]] <- c(paste(structure2target_file, structure_file, sep = " "), paste("NA", target_file, sep = " "));

## 3. Generate commands for jobs 
commands <- generate.all.possible.hyperparameter.combinations(part);
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
