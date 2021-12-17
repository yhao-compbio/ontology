# This folder contains source code used by the repository.

## R/python scripts 

+ [`generate_root_file.R`](generate_root_file.R) generates root pathway files and number codes for constructing DTox hierarchy.

+ [`get_knowledge_hierarchy.R`](get_knowledge_hierarchy.R) processes raw ontology data to construct DTox hierarchy. [`run/run_get_knowledge_hierarchy_reactome.R`](run/run_get_knowledge_hierarchy_reactome.R) generates shell scripts that processes reactome pathway data to construct DTox hierarchy.

+ [`shuffle_hierarchy.R`](shuffle_hierarchy.R) shuffles the hierarchy of input DTox network while preserving the distribution of connections between node layers. [`run/run_shuffle_hierarchy.R`](run/run_shuffle_hierarchy.R) generates shell scripts that shuffles the hierarchy of optimal DTox model for tox21 datasets.

+ [`functions.R`](functions.R) contains basic functions to process pathway ontology data.

## Executable shell scripts

+ [`run/get_knowledge_hierarchy_reactome.sh`](run/get_knowledge_hierarchy_reactome.sh) runs [`get_knowledge_hierarchy.R`](get_knowledge_hierarchy.R) on reactome pathway data and 15 different root pathway combinations of interest.

+ [`run/shuffle_hierarchy_reactome.sh`](run/shuffle_hierarchy_reactome.sh) runs [`shuffle_hierarchy.R`](shuffle_hierarchy.R) on optimal DTox model for tox21 datasets.
