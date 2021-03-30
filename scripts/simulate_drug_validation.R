   #############################################################################
#### Simulate the effect of the list of drug  and compare with the patients  ######
   #############################################################################

library(pacman)
devtools::install_github("thomasp85/patchwork")

pacman::p_load("here", "hipathia", "utils", "stringr", "genefilter", "tidyr","data.table","limma", "ggplot2", "patchwork")

## Load normalized expression data
logcpm <- readRDS(file = here("rds", "logcpm_validating_cohorts.rds"))
metadata_init <- fread(file = here("data", "feedback_guillermo_Marzo2021", "metadata_reformatted.tsv"))

## Subset only disease samples
logcpm_dis <- logcpm[, which(colnames(logcpm) %in% metadata_init$samples[metadata_init$diagnosis != "CTRL"])]

trans_data <- translate_data(logcpm_dis, "hsa")
saveRDS(trans_data, file = here("rds", "trans_data_onlydiseases_validationCohort.rds"))

metadata_init <- metadata_init[metadata_init$diagnosis != "CTRL",]

metadata_init$type <- "controls"

## Loading Pathways (all and only physiological)

path_list <- read.table(file = here("data", "physiological_paths.tsv"), sep = "\t") #physiological_pathways
pathways <- load_pathways("hsa")
pathways_phy <- load_pathways("hsa", pathways_list = path_list$V2 )

subpathways.list <- lapply(pathways_phy$pathigraphs, function(x){names(x$effector.subgraphs)})

subpathways_phy <- data.frame(unlist(subpathways.list), stringsAsFactors = F)

colnames(subpathways_phy) <- "hipathia"


#### DO KO , CALCULATE HIPATHIA AND COMPARE THE log2 FC and plot it ###############

## List of desired drugs to test
drugs_totest <- read.delim(file = here("data", "feedback_guillermo_Marzo2021", "drugs_for_testing.tsv"), header = F, sep = " ")
colnames(drugs_totest) <- c ( "drug_name", "drug_type")

drugbank <- read.delim(file = "/home/m3m/INFO_PROYECTO/drugbank/drugbank_drug-bindings_v5.1.8.tsv", header = T)
drugbank <- drugbank[drugbank$drug_binding == "target_gene",]

## Check how many drug names are found in the drugbank 5.1.8 db 
sum(drugs_totest$drug_name %in% drugbank$drug_name) ## 18/18 ALL!! ^^ 

## Get gene targets from the desired drugs

drug_targetsDF <- drugbank[drugbank$drug_name %in% drugs_totest$drug_name, c(2,4,5,8,9,10,12,13,14)]
drug_targetsDF <- drug_targetsDF[which(!is.na(drug_targetsDF$entrez_id)), ]

targets <- unique(drug_targetsDF$entrez_id)    

length(targets) ## 30 gene targets    

sum(drugs_totest$drug_name %in% drug_targetsDF$drug_name) ### Check if all the drugs have a targ    et with and entrez ID


## Load Hipathia and chech which targets are in the pathways object

sum(targets %in% pathways$all.genes) ## 18/31 en pathways_phy y 20/30 en pathways

drug_targetsDF$in_Hipathia <- drug_targetsDF$entrez_id %in% pathways$all.genes

drug_targetsDF <- drug_targetsDF[drug_targetsDF$in_Hipathia == TRUE , ]    

sum(drugs_totest$drug_name %in% drug_targetsDF$drug_name) ## 12/18 can be tested in Hipathia

## Drugs to test which are in Hipathia .. Export for future reports
drugs_totest[drugs_totest$drug_name %in% drug_targetsDF$drug_name, ]

# write.table(drug_targetsDF, file = here("results", "drugstotest_inHipathia.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")


### 2. USE MADE FUNCTIONS FOR THE KO -- SIMULATING THE DRUG ####

KO_drogas <- lapply(drugs_totest$drug_name[-1], function(x){simulate_drug_validation(drug = x, folder = "on_2dataset", expression_matrix = trans_data)})


