# !/usr/bin/env python
## created by Yun Hao @MooreLab 2021
## This script generates shell scripts that shuffles the hierarchy of optimal DTox model for tox21 datasets


## functions 
source("src/functions.R");


## 0. Input arguments
perf_file		<- "/home/yunhao1/project/DTox/data/compound_target_probability_tox21_implementation/compound_target_probability_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv";
reaction_structure	<- "re_0_st_0"; 
input_folder		<- "data/reactome/hierarchy/"
output_folder		<- "data/reactome/hierarchy_shuffle/";
N_cores 		<- 1;
job_name		<- "shuffle_hierarchy_reactome"

## 1. Obtain hyperparameters of input DTox hierarchy for datasets
perf_df <- read.delim(file = perf_file, header = T, sep = "\t");
perf_hp <- sapply(perf_df$hyperparameter_setting, function(pdhs) strsplit(pdhs, "_xs_", fixed = T)[[1]][[1]]);

## 2. Generate commands for jobs 
# iterate by hyperparameters for each dataset
commands <- mapply(function(ph){
	# node connection file by node number 
	hp_knowledge <- paste(input_folder, ph, "_st_0_knowledge_by_node.tsv", sep = ""); 
	# node layer file
	hp_layer <- paste(input_folder, ph, "_st_0_layer.tsv", sep = "");
	# output file 
	hp_op <- paste(output_folder, ph, "_st_0_knowledge_by_node_shuffle.tsv", sep = "");
	# command  
	hp_command <- paste("Rscript", "src/shuffle_hierarchy.R", hp_knowledge, hp_layer, hp_op, sep = " ");
	return(hp_command);
}, perf_hp);
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
