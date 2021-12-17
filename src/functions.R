# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2021
## This script contains basic functions to process pathway ontology data


## This function generates all possible combinations of input elements, as well as a number code for each element combination (first number represent length, i.e. number of elements, second number represent order in the corresponding length group, e.g.: 11, 12, 13, 21, 22, 23, 33 for combo of 3). 
generate.all.element.combinations <- function(ele_vec){
	## 0. Input arguments 
		# ele_vec: vector that contains elements 
	
	## 1. Generate list containing all possible combinations of elements   
	# obtain the number of elements  
	ev_len <- length(ele_vec);
	# iterate by all possible lengths of element combinations 
	ele_combo_list <- lapply(1:ev_len, function(el){
		el_combo <- combn(ev_len, el);
		el_combo_list <- lapply(1:ncol(el_combo), function(nec) ele_vec[el_combo[,nec]])
		return(el_combo_list);
	});
	ele_combo_list1 <- unlist(ele_combo_list, recursive = F); 	

	## 2. Generate a number code for each element combination
	# obtain the number of combinations at each possible length, and the number of digit of first number 
	ecl_len <- sapply(ele_combo_list, length);
	n_all_char <- nchar(max(ecl_len));
	# obtain the first number of each combination, representing length
	ecl_combo_len <- mapply(function(lel,el) rep(lel,el), 1:length(ecl_len), ecl_len);
	ecl_combo_len <- unlist(ecl_combo_len);
	# obtain the second number of each combination, representing order  
	ecl_order <- lapply(ecl_len, function(el) 1:el);
	ecl_order <- unlist(ecl_order);
	# put digits of two numbers together
	ecl_id_char <- mapply(function(ecl, eo){
		n_zero <- n_all_char - nchar(eo);
		zero_char <- paste(rep('0', n_zero), collapse = "");
		id_char <- paste(ecl, zero_char, eo, sep = "");
		return(id_char);
	}, ecl_combo_len, ecl_order);
	names(ele_combo_list1) <- ecl_id_char;

	return(ele_combo_list1);
}


## This function filters pathways by size and relationship with root pathways  
filter.pathways <- function(roots, relation_mat, annotation_mat, min_size){
	## 0. Input arguments 
		# roots: root pathways  
		# relation_mat: matrix that contains parent/child relationships between pathways (1st column: parent; 2nd column: child) 
		# annotation_mat: matrix that contains gene annotation of pathways (1st column: pathway; 2nd column: gene)
		# min_size: minimal size of pathway to be included in hierarchy  

	## 1. Remove pathways that do not satisfy minimal size threshold 
	# obtain names of pathways that pass minimal size threshold  
	path_size <- table(annotation_mat[,2]);
	path1 <- names(path_size)[which(path_size >= min_size)];
	# only keep the pathways that pass thershold in the relation matrix
	rela_id1 <- which(relation_mat[,1] %in% path1);
	rela_id2 <- which(relation_mat[,2] %in% path1);		
	rela_id <- intersect(rela_id1, rela_id2);
	relation_mat <- relation_mat[rela_id, ];

	## 2. Remove pathways that are not descendants of root pathways 
	# start with root pathways
	parents <- roots;
	children <- logical(0);
	# add children of current pathways to current pathways, until no more can be added
	while (length(parents) != length(children)){
		parents <- unique(union(parents, children))
		child_id <- which(relation_mat[,1] %in% parents)
		children <- union(parents, unique(relation_mat[child_id, 2]))
	}
	# only keep descendants of root pathways in the relation matrix
	relation_mat <- relation_mat[child_id, ];
	# only keep descendants of root pathways in the annotation matrix
	annot_id <- which(annotation_mat[,2] %in% children);
	annotation_mat <- annotation_mat[annot_id, ];

	return(ls = list(relation = relation_mat, annotation = annotation_mat));
}


## This function groups a vector by categories of its elements.
group.vector.by.categories <- function(cate, vec){	
	## 0. Input arguments:
		# cate: category of vectors
		# vec: vector
	
	## 1. Sort vector by order of categorie
	vec_od <- order(cate);
	cate <- cate[vec_od]
	vec <- vec[vec_od];

	## 2. Group elements of same category together
        # obtain unique categorie
	cate_table <- table(cate);
	# obtain lower bound index/upper bound index of each unique category
	lower_ids <- cumsum(cate_table);
	upper_ids <- c(0, lower_ids[-length(lower_ids)]) + 1;
	# return list of vectors 
	vec_list <- mapply(function(li, ui) vec[li:ui], lower_ids, upper_ids, SIMPLIFY=F);
        names(vec_list) <- names(cate_table);

	return(vec_list);
}


## This function performs child-parent mapping between reactions and pathways, then filters reactions by the mapping 
filter.reactions <- function(react_mat, relation_mat, annotation_mat){
	## 0. Input arguments: 
		# react_mat: matrix that contains gene annotation of reactions
		# relation_mat: matrix that contains parent/child relationships between pathways
		# annotation_mat: matrix that contains gene annotation of pathways
	
	## 1. Obtain gene annotations of lowest-level pathways and reactions 
	# obtain pathways at the lowest level of hierarchy (no children pathways)
	lowest_pathways <- unique(setdiff(relation_mat[,2], relation_mat[,1]));
	# obtain annotated genes of each lowest-level pathway
	lp_id <- which(annotation_mat[,2] %in% lowest_pathways);
	lp_gene_list <- group.vector.by.categories(annotation_mat[lp_id,2], annotation_mat[lp_id,1]);
	# obtain annotated genes of each reaction 
	react_gene_list <- group.vector.by.categories(react_mat[,2], react_mat[,1])
	
	## 2. Perform child-parent relationship mapping between reactions and lowest-level pathways 
	# iterate by reaction  
	react_lp_list <- mapply(function(nrgl, rgl){
		# find the intersection between current reaction and each lowest-level pathway 
		lp_inter <- sapply(lp_gene_list, function(lgl) length(intersect(lgl, rgl)));
		# assign reaction as child of a lowest-level pathway if the reaction is a subset of the lowest-level pathway
		p_id <- which(lp_inter == length(rgl));
        	if(length(p_id) == 0)	p_mat <- matrix(NA, 0, 0)
		else{
			# output child/parent relationship in matrix form
			p_mat <- matrix(NA, length(p_id), 2)
			p_mat[,1] <- names(lp_gene_list)[p_id]
			p_mat[,2] <- nrgl
		}
		return(p_mat);
	}, names(react_gene_list), react_gene_list, SIMPLIFY = F);
	# remove reactions that are mapped to more than one lowest-level pathway 
	react_lp_len <- sapply(react_lp_list, nrow);
	react_lp_list <- react_lp_list[react_lp_len == 1];
	react_lp_mat <- do.call(rbind, react_lp_list);

	## 3. Filter the relation and annotation matrix after mapping
	# start with lowest-level pathways that are mapped with reactions  
	children <- unique(react_lp_mat[,1]);
	parents <- logical(0);
	# add parents of current pathways to current pathways, until no more can be added  
	while (length(children) != length(parents)){
		children <- unique(union(children, parents))
		c_id <- which(relation_mat[,2] %in% children)
		parents <- union(children, unique(relation_mat[c_id,1]))
	}
	# only keep ancestors of mapped lowest-level pathways in the relation matrix 
	relation_mat <- relation_mat[c_id,];	
	combine_relation_mat <- rbind(react_lp_mat, relation_mat);
	# only keep ancestors of mapped lowest-level pathways in the pathway annotation matrix 
	a_filter_id1  <- which(annotation_mat[,2] %in% children);
	annotation_mat <- annotation_mat[a_filter_id1, ]; 
	# only keep mapped reactions in the reaction annotation matrix
	lowest <- unique(react_lp_mat[,2]);
	a_filter_id2 <- which(react_mat[,2] %in% lowest);
	react_mat <- react_mat[a_filter_id2, ];
	# combine filtered pathway and reaction annotation matrix
	combine_annotation_mat <- rbind(react_mat, annotation_mat);	
	
	return(ls = list(lowest = react_lp_mat[,2], relation = combine_relation_mat, annotation = combine_annotation_mat));		
}


## This function filter pathways by target genes of interest
filter.targets <- function(all_targets, lowest, relation_mat, annotation_mat){
	## 0. Input arguments
		# all_targets: vector that contains target genes of interest
		# lowest: vector that contains lowest-level pathways
		# relation_mat: matrix that contains parent/child relationships between pathways
		# annotation_mat: matrix that contains gene annotation of pathways
	
	## 1. Identify targets for DTox hierarchy 
	# obtain annotated genes of lowest-level pathways 
	lowest_id <- which(annotation_mat[,2] %in% lowest);
	lowest_mat <- annotation_mat[lowest_id, ];
	# find overlaps between targets of interest and target genes annotated in lowest-level pathways 
	inter_target <- intersect(lowest_mat[,1], all_targets);

	## 2. Filter the relation and annotation matrix by identified targets
	# identify lowest-level pathways that are annotated with identified targets of DTox hierarchy 
	it_id <- which(lowest_mat[,1] %in% inter_target);
	inter_lowest <- unique(lowest_mat[it_id,2]);
	# start with lowest-level pathways that are annotated with identified targets of DTox hierarchy 
	children <- inter_lowest;
	parents <- logical(0);
	# add parents of current pathways to current pathways, until no more can be added   
	while (length(children) != length(parents)){
		children <- unique(union(children, parents))
		c_id <- which(relation_mat[,2] %in% children)
		parents <- union(children, unique(relation_mat[c_id,1]))
	}
	# only keep ancestors of target lowest-level pathways in the relation matrix
	relation_mat <- relation_mat[c_id, ]
	# only keep ancestors of target lowest-level pathways in the annotation matrix
	annot_id <- which(annotation_mat[,2] %in% children);
	annotation_mat <- annotation_mat[annot_id, ];
	
	return(ls = list(inter_target = inter_target, inter_lowest = inter_lowest, relation = relation_mat, annotation = annotation_mat));
}


## This function builds connections between target genes of interest, lowest-level pathways, and higher level pathways to form DTox hierarchy 
get.knowledge.hierarchy <- function(inter_target, lowest, annotation_mat, relation_mat){
	## 0. Input arguement 
		# inter_target: vector that contains target genes included in the DTox hierarchy  
		# lowest: vector that contains lowest-level pathways
		# relation_mat: matrix that contains parent/child relationships between pathways
		# annotation_mat: matrix that contains gene annotation of pathways 

	## 1. Build connections between target genes of interest and lowest-level pathways
	# identify annotated target genes of interest for each lowest-level pathway
	it_id <- which(annotation_mat[,1] %in% inter_target);
	l_id <- which(annotation_mat[,2] %in% lowest);
	itl_id <- intersect(it_id, l_id);	
	lowest_target <- group.vector.by.categories(annotation_mat[itl_id,2], annotation_mat[itl_id,1]);
	# separate each annotated target gene of interest by ',', store in data frame form  
	lowest_target_vec <- sapply(lowest_target, function(lt) paste(lt, collapse = ","));
	lowest_target_mat <- data.frame(names(lowest_target), lowest_target_vec);

	## 2. Build connections between higher level pathways
	# identify children pathways for each higher level pathway  
	path_children <- group.vector.by.categories(relation_mat[,1], relation_mat[,2]);
	# separate each children pathway by ',', store in data frame form
	path_children_vec <- sapply(path_children, function(pc) paste(pc, collapse = ","));
	path_children_mat <- data.frame(names(path_children), path_children_vec);
	
	## 3. Combine two connection data frames to form DTox hierarchy
	colnames(lowest_target_mat) <- colnames(path_children_mat) <- c("V1", "V2")
	knowledge_mat <- rbind(lowest_target_mat, path_children_mat);

	return(knowledge_mat);
}


## This function uses child/parent relationships to sort pathways in the DTox hierarchy by layer
sort.nodes.by.hierarchy <- function(lowest, relation_mat, path_layer){
	## 0. Input arguments 
		# lowest: vector that contains lowest-level pathways included in DTox hierarchy
		# relation_mat: matrix that contains parent/child relationships between pathways
		# path_layer: list that contains nodes in the DTox hierarchy that have already been sorted (structure layer if available, target gene layer)
	
	## 1. Sort pathways by the order of layer that no child appears in a layer after its parent 
	# assign lowest-level pathways to the next available layer 
	i <- length(path_layer) + 1;
	path_layer[[i]]	<- lowest;
	# start with lowest-level pathways 
	sorted_paths <- lowest;
	all_paths <- unique(c(relation_mat[,1], relation_mat[,2]));
	path_children <- group.vector.by.categories(relation_mat[,1], relation_mat[,2]);
	# add qualified parents of current layer pathways to the next available layer, until no more can be added
	while(length(sorted_paths) != length(all_paths)){
		# identify all parents of current layer pathways
		c_id <- which(relation_mat[,2] %in% path_layer[[i]])
		parent_paths <- unique(relation_mat[c_id, 1])
		# check if all children of identified parent pathways have been sorted in the previous layers
		parent_path_id <- sapply(path_children[parent_paths], function(pcpp){
			pcpp_diff <- setdiff(pcpp, sorted_paths)
			if(length(pcpp_diff) == 0)	return(1)
			else	return(0)
		})
		# only assign the parent pathways whose children have been all sorted to the next available layer  
		i <- i + 1
		path_layer[[i]] <- parent_paths[parent_path_id == 1]
		# update pathways that have been sorted 
		sorted_paths <- union(sorted_paths, path_layer[[i]])
	}
	
	## 2. Number the layer and node in the sorted order
	# Number each layer in the sorted order, and store the layer number of each node in data frame (layer number starts with 0)
	layer_number <- mapply(function(pl, lpl){
		rep(lpl-1, length(pl))
	}, path_layer, 1:length(path_layer), SIMPLIFY = F);
	layer_number_vec <- unlist(layer_number);
	layer_mat <- data.frame(1:length(layer_number_vec)-1, layer_number_vec);
	colnames(layer_mat) <- c("node", "layer_number");
	# Number each node in the order of layer number, store node number in data frame (node number starts with 0) 
	node_name_vec <- unlist(path_layer)
	node_mat <- data.frame(1:length(node_name_vec)-1, node_name_vec);
	colnames(node_mat) <- c("node", "node_name");

	return(ls = list(hierarchy = path_layer, layer = layer_mat, node = node_mat));	
}


## This function converts entities in knowledge matrix to their corresponding node numbers in DTox hierarchy  
convert.node.hierarchy.to.numbers <- function(knowledge_mat, layer_mat, node_mat){
	## 0. Input arguments 
		# knowledge_mat: matrix that contains the structure features (if availabel), target genes, or pathways that each entity is connected to in DTox hierarchy
		# layer_mat: matrix that contains layer number of each node in DTox hierarchy  
		# node_mat: matrix that contains number of each node in DTox hierarchy  


	## 1. Sort pathways in knowledge matrix by node the order in DTox hierarchy  
	layer_id <- which(layer_mat$layer_number > 0);
	km_id <- sapply(layer_id, function(li) which(knowledge_mat[,1] %in% node_mat$node_name[[li]]));
	knowledge_mat1 <- knowledge_mat[km_id, ];
	colnames(knowledge_mat1) <- c("node_name", "children_node_name");

	## 2. Map entities in knowledge matrix to their node numbers
	# build a vector that contains mapping between node number and name  
	node_name_map <- node_mat$node;
	names(node_name_map) <- node_mat$node_name;
	# map target genes/children pathways of each pathway to their node numbers 
	knowledge_number <- sapply(knowledge_mat1[,2], function(km2){
		km2_s <- strsplit(km2, ",")[[1]];
		kms_number <- sort(node_name_map[km2_s]);
		kms_np <- paste(kms_number, collapse = ",");
		return(kms_np);
	});
	# store mapped node numbers in data frame form 
	knowledge_number_df <- data.frame(layer_mat$node[layer_id], knowledge_number);
	colnames(knowledge_number_df) <- c("node", "children_node");

	return(ls = list(origin = knowledge_mat1, number = knowledge_number_df));
} 


## This function obtains information about pathway (size) and root nodes (name) in DTox hierarchy 
get.node.information <- function(node_mat, annotation_mat, root){
	## 0. Input arguments 
		# node_mat: matrix that contains number of each node in DTox hierarchy   
		# annotation_mat: matrix that contains gene annotation of pathways  
		# root: vector that contains root pathways 
	
	## 1. Obtain sizes of pathways in DTox hierarchy
	# compute the number of annotated genes for each pathway 
	all_path_genes <- group.vector.by.categories(annotation_mat[,2], annotation_mat[,1]);
	all_path_len <- sapply(all_path_genes, length);
	# obtain the number of annotated genes for pathways included in DTox hierarchy
	node_size <- rep(1, nrow(node_mat));
	names(node_size) <- node_mat$node_name;
	node_size[names(all_path_len)] <- all_path_len;
	# store pathway size of nodes in data frame form
	node_size_df <- data.frame(node_mat$node, node_size);
	colnames(node_size_df) <- c("node", "size");
	
	## 2. Obtain names of root pathways in DTox hierarchy
	root_id <- which(node_mat$node_name %in% root);
	root_df <- node_mat[root_id, ];
	colnames(root_df) <- c("root", "root_name");
	
	return(ls = list(size = node_size_df, root = root_df));
}


## This function generates all possible combinations for a list of hyperparameters 
generate.all.possible.hyperparameter.combinations <- function(hp_list){
	## 0. Input arguments 
		# hp_list: list of hyperparameters   

	## 1. Generate all possible hyperparameter combinations
	# iterate by hyperparameters 
	hp_current <- hp_list[[1]];
	for (i in 2:length(hp_list)){
		# next hyperparameter
		hp_next <- hp_list[[i]];
		# combine two hyperparmeters  
		hp_current_vec <- rep(hp_current, each = length(hp_next));
		hp_next_vec <- rep(hp_next, times = length(hp_current));
		hp_current <- mapply(function(hcv, hnv) paste(hcv, hnv, sep = " "), hp_current_vec, hp_next_vec)
	}

	return(hp_current)
}


## This function generates executable shell scripts that will run input commands 
generate.parallel.bash.files <- function(all_commands, N_group, job_name, folder){
	## 0. Input arguments 
		# all_commands: a vector of commands that are to be run
		# N_group: number of groups that commands will be split into  
		# job_name: name of job 
		# folder: name of folder where shell scripts will be written 

	## 1. Assign the indices of commands to each group 
	# obtain the number of commands in each group
	N_group_member <- ceiling(length(all_commands) / N_group);
	# obtain upper bound of index for each group 
	upper_bound <- 1:N_group * N_group_member;
	ub_id <- min(which(upper_bound >= length(all_commands)));
	upper_bound <- upper_bound[1:ub_id];
	upper_bound[[ub_id]] <- length(all_commands);
	# obtain lower bound of index for each group 
	lower_bound <- 0:(ub_id-1) * N_group_member + 1;
	# assign commands to each group (lower bound - upper bound)
	command_list <- mapply(function(lb, ub) all_commands[lb:ub], lower_bound, upper_bound, SIMPLIFY = F);

	## 2. write commands into executable shell scripts
	# name executable shell scripts by "job_nameX.sh" (X is the index of script)
	setwd(folder);
	c_file_name <- sapply(1:ub_id, function(gn) paste(job_name, gn, ".sh", sep = ""));
	# iterate by script
	write_sub_files <- mapply(function(cl, cfn){
		# write commands into the script
		writeLines(cl, cfn);
		# make the script executable
		system(paste("chmod", "775", cfn, sep=" "));
		return(1);
	},command_list,c_file_name,SIMPLIFY=F);

	## 3. write an executable shell script that runs all the scripts above  
	final_command <- sapply(c_file_name,function(cfn) paste("./", folder, cfn," &", sep = ""));
	# name executable shell scripts by "job_name.sh" 
	final_file <- paste(job_name, ".sh", sep = "");
	writeLines(final_command, final_file);
	# make the script executable
	system(paste("chmod","775", final_file, sep = " "));

	return("DONE");
}


