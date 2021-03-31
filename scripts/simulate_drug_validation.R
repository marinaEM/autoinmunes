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

KO_drogas <- lapply(drugs_totest$drug_name, function(x){simulate_drug_validation(drug = x, folder = "on_2dataset", expression_matrix = trans_data)})

names(KO_drogas) <- drugs_totest$drug_name

KO_drogas <- KO_drogas[!sapply(KO_drogas, is.null)]

saveRDS(KO_drogas, file = here("rds", "KO_drugs", "on_2dataset", "FC_df_12drugs.rds"))


### RE DO the plots for those drugs where the line is not well positioned

### Sifalimumab
        # drug <- "Sifalimumab"
        # folder <- "on_2dataset/reDO"
        # mediana <- 0.0005
        # FC_df <- KO_drogas$Sifalimumab

### Tozilizumab

        drug <- "Tocilizumab"
        folder <- "on_2dataset/reDO"
        mediana <- 0.00007
        FC_df <- KO_drogas$Tocilizumab

        
        FCs <- ggplot(data = FC_df,
                      mapping = aes(x = ranked_patients,
                                    y = logFC,
                                    color = Cluster)) +
          ggtitle(paste0('Log FC of the circuits activity absolute values on the ranked patients after simulation of ',drug ,' drug effect' )) +
          geom_point(size = 3) +
          geom_hline(yintercept = mediana, linetype = "dashed") +
          theme_minimal()+
          theme(axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                legend.text = element_text(size = 18),
                legend.title =  element_text(size = 24),
                plot.title = element_text(size = 20, face = "bold"))
        
        FCs_noCluster2 <- ggplot(data = FC_df[!FC_df$Cluster== "2",],
                                 mapping = aes(x = ranked_patients,
                                               y = logFC,
                                               color = Cluster)) +
          ggtitle(paste0('Log FC of the circuits activity absolute values on the ranked patients after simulation of ',drug ,' drug effect without Cluster 2' )) +
          geom_point(size = 3) +
          geom_hline(yintercept = mediana, linetype = "dashed") +
          theme_minimal()+
          theme(axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                legend.text = element_text(size = 18),
                legend.title =  element_text(size = 24),
                plot.title = element_text(size = 20, face = "bold"))
        
        
        
        ## With CLUSTER2 
        bars_df_res <- data.frame(table(FC_df$Cluster[FC_df$logFC >= mediana]))
        
        bars_df_NOres <- data.frame(table(FC_df$Cluster[FC_df$logFC < mediana]))
        
        bars_responders <- ggplot(data = bars_df_res, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
          ggtitle("High response patients")+
          geom_bar(stat="identity") +
          theme_minimal() +
          theme(axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                legend.title = element_blank(),
                legend.position = "none") +
          coord_flip()
        
        
        bars_nonresponders <- ggplot(data = bars_df_NOres, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
          ggtitle("Low response patients")+
          geom_bar(stat="identity") +
          theme_minimal() +
          theme(axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                legend.title = element_blank(),
                legend.position = "none") +
          coord_flip()
        
        ## NO CLUSTER2 
        FCdf_noC2 <- FC_df[FC_df$Cluster!= "2",] 
        
        bars_df_res_noCluster2 <- data.frame(table(FCdf_noC2$Cluster[FCdf_noC2$logFC >= mediana]))
        bars_df_NOres_noCluster2 <- data.frame(table(FCdf_noC2$Cluster[FCdf_noC2$logFC < mediana]))
        
        bars_responders_noC2 <- ggplot(data = bars_df_res_noCluster2, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
          ggtitle("High response patients")+
          geom_bar(stat="identity") +
          theme_minimal() +
          theme(axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                legend.title = element_blank(),
                legend.position = "none") +
          coord_flip()
        
        
        bars_nonresponders_noC2 <- ggplot(data = bars_df_NOres_noCluster2, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
          ggtitle("Low response patients")+
          geom_bar(stat="identity") +
          theme_minimal() +
          theme(axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                legend.title = element_blank(),
                legend.position = "none") +
          coord_flip()
        
        
        png(filename = here("results", "KO_drugs", folder, paste0("FCplot_validation_",drug,"_KO_clusters.png")), width = 5000, height = 3000, res = 200)
        print((FCs | bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1))) 
        dev.off()
        
        png(filename = here("results", "KO_drugs", folder,paste0("FCplot_validation",drug,"_KO_withoutCluster2.png")), width = 5000, height = 3000, res = 200)
        print((FCs_noCluster2 | bars_responders_noC2/plot_spacer()+ bars_nonresponders_noC2) + plot_layout(widths = c(3, 1))) #+ plot_layout(guides = 'collect')  
        dev.off()

