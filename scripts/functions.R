      #########################################################################
################# FUNCTIONS FOR THE AUTOINMUNES PROJECT ############################
      ######################################################################

#################################################################
###### Function to apply to each disease vs control comparison 
########     (previously dividing the matrixes)  16.03.2021
##################################################################  
      
      get_diff_expressed_pathways_diseases <- function(matrix_metadata, pathvals){
        
        disease <- names(table(matrix_metadata$diagnosis))[names(table(matrix_metadata$diagnosis)) != "CTRL"]
        message(paste0("performing limma and extracting significant circuits for... ", disease))
        
        ## Select the columns corresponding to each disease
        samples <- matrix_metadata$samples
        idx <- which(colnames(pathvals) %in% samples)
        disease_path_vals <- pathvals[, idx]
        
        ##  Disease vs control limma
        sample_group_limma <- matrix_metadata[, c(1,2,4,5,12,13,14)]
        design <- model.matrix( ~ 0  + type + age + gender + pool + RIN , data = sample_group_limma ) #,+ antimalarials + immunosuppresants + biologicals + steroids + systemic.antibiotics ,  data = sample_group_limma)
        colnames(design) <- gsub(pattern = "typeDISEASE", replacement ="disease", x = colnames(design))
        colnames(design) <- gsub(pattern = "typeCONTROLS", replacement ="C", x = colnames(design))
        rownames(design) <- sample_group_limma$samples
        cont.matrix <- makeContrasts( diseasevsC=disease-C, levels=design)
        
        fit <- lmFit(disease_path_vals , design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit3 <- eBayes(fit2)
        tableLimma <- topTable(fit3, number = rownames(path_vals), adjust.method="bonferroni", sort.by="p")
        top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
        top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
        
        anot_all <- data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
                               logFC = top_LimmaALL$logFC,
                               P.val = top_LimmaALL$P.Value,
                               FDR.Pval = top_LimmaALL$adj.P.Val,
                               UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
                               GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
                               stringsAsFactors = F)
        
        # write.xlsx(anot_all, file = here( "results", paste0(disease, "vsC_afterlimma_annot.xlsx")),row.names = F)  
        
        pathways_dereg <- str_split(anot_all$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.)
        
        message(length(pathways_dereg))
        message(length(anot_all$Circuit))
        
        write.xlsx(pathways_dereg, file = here( "results", paste0(disease, "pathways_dereg.xlsx")),row.names = F, col.names = F) 
        
        return(anot_all)
      }
      
      
#################################################################
###### Function to apply to each cluster vs control comparison 
########     (previously dividing the matrixes)  16.03.2021
##################################################################  
      
      get_diff_expressed_pathways_cluster <- function(matrix_metadata, pathvals){
        
        cluster <- names(table(matrix_metadata$cluster))[names(table(matrix_metadata$cluster)) != "CONTROLS"]
        message(paste0("performing limma and extracting significant circuits for... ", cluster))
        
        ## Select the columns corresponding to each cluster
        samples <- matrix_metadata$samples
        idx <- which(colnames(pathvals) %in% samples)
        cluster_path_vals <- pathvals[, idx]
        
        ##  Disease vs control limma
        sample_group_limma <- matrix_metadata[, c(1,2,4,5,12,13)]
        design <- model.matrix( ~ 0  + cluster + age + gender + pool + RIN, data = sample_group_limma ) #, + systemic.antibiotics ,  data = sample_group_limma)
        colnames(design) <- gsub(pattern = "clusterCLUSTER[0-9]", replacement ="disease", x = colnames(design))
        colnames(design) <- gsub(pattern = "clusterCONTROLS", replacement ="C", x = colnames(design))
        rownames(design) <- sample_group_limma$samples
        cont.matrix <- makeContrasts( diseasevsC=disease-C, levels=design)
        
        fit <- lmFit(cluster_path_vals , design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit3 <- eBayes(fit2)
        tableLimma <- topTable(fit3, number = rownames(path_vals), adjust.method="bonferroni", sort.by="p")
        top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
        top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
        
        anot_all <- data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
                               logFC = top_LimmaALL$logFC,
                               P.val = top_LimmaALL$P.Value,
                               FDR.Pval = top_LimmaALL$adj.P.Val,
                               UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
                               GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
                               stringsAsFactors = F)
        
         write.xlsx(anot_all, file = here( "results", paste0(cluster, "vsC_afterlimma_annot.xlsx")),row.names = F)  
        
        pathways_dereg <- str_split(anot_all$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.)
        
        message(length(pathways_dereg))
        message(length(anot_all$Circuit))
        
        write.xlsx(pathways_dereg, file = here( "results", paste0(cluster, "pathways_dereg.xlsx")),row.names = F, col.names = F) 
        
        return(anot_all)
      }

############################################################################################
#### Function to simulate the effect of a drug by KO of single targets, previously 
#######  load trans_data, pathways and metadata without controls - PLOT of log2 FC 16.03.2021
###############################################

simulate_drug_targetonebyone <- function(target, expression_matrix = trans_data, metadata_matrix_controls = metadata_init, 
                                         drug_table = drug_targetsDF, pathways_db = pathways, folder = NA, PublicID = PublicID){
  if( is.na(folder)) {
    folder = ""
     }else{
    message(paste0("Folder where the results will be saved... ", folder))
  }
  
  
  ## Check that the target is in the pathways
  
  if(target %in% pathways_db$all.genes){
    
    ## Define drug name 
    drug <- drug_table$drug_name[match(target, drug_table$entrez_id)]
    
    target_name <- drug_table$gene_name[match(target, drug_table$entrez_id)]
    
    message(paste0("Simulating...", drug, " with KO of ", target_name," in expression dataset"))
    
    metadata_KO <- metadata_matrix_controls
    metadata_KO$PublicID <- paste0(metadata_KO$PublicID,"_" ,target_name)
    metadata_KO$type <- "treated"
    
    trans_data_KO <- expression_matrix
    trans_data_KO[target, ] <- 0.0001
    colnames(trans_data_KO) <- paste0(colnames(trans_data), "_", target_name)
    
    ## Normalize expression data all together (controls + repeated expression matrix with KO)
    trans_data_KO_C <- cbind(trans_data_KO, expression_matrix)
    exp_data_KO_C <- normalize_data(trans_data_KO_C)
    
    message("Calculating Hipathia")
    ## Using Hipathia to compute the signal (subpathways) 
    results_KO <- hipathia(exp_data_KO_C, pathways_db, decompose = FALSE, verbose=FALSE)
    saveRDS(results_KO, file = here("rds","KO_drugs", folder, paste0(drug,"_KO_", target_name ,"_results.rds")))
    
    message("saving results of KO")
    
    path_vals_KO <- get_paths_data(results_KO, matrix = TRUE)
    # path_vals_phy_KO <- path_vals_KO[which(rownames(path_vals_KO) %in% subpathways_phy$hipathia),]
    # anot_cir <- data.frame(Circuit = get_path_names(pathways, rownames(path_vals_phy_KO)),
    #                        UniprotKB = get_pathways_annotations(rownames(path_vals_phy_KO), pathways, dbannot= "uniprot" ,collapse = T),
    #                        GO = get_pathways_annotations(rownames(path_vals_phy_KO), pathways, dbannot= "GO" ,collapse = T), 
    #                        stringsAsFactors = F)
    
    
    ## LIMMA FOR ALL SAMPLES 
    message("Building up Limma models")
    metadata <- rbind(metadata_KO, metadata_matrix_controls)
    target_matrix_KO <- data.frame(samples = colnames(path_vals_KO), 
                                   cluster = metadata$DISEASE[match(colnames(path_vals_KO), metadata$PublicID)],
                                   diagnosis = metadata$Diagnosis[match(colnames(path_vals_KO), metadata$PublicID)],
                                   age = metadata$Age[match(colnames(path_vals_KO), metadata$PublicID)],
                                   gender = metadata$Gender[match(colnames(path_vals_KO), metadata$PublicID)],
                                   antimalarials = metadata$Antimalarials[match(colnames(path_vals_KO), metadata$PublicID)],
                                   immunosuppresants = metadata$Immunosuppresants[match(colnames(path_vals_KO), metadata$PublicID)],
                                   biologicals = metadata$Biologicals[match(colnames(path_vals_KO), metadata$PublicID)],
                                   steroids = metadata$Steroids[match(colnames(path_vals_KO), metadata$PublicID)],
                                   systemic.antibiotics = metadata$Systemic.Antibiotics[match(colnames(path_vals_KO), metadata$PublicID)],
                                   public.center = metadata$PublicCenter[match(colnames(path_vals_KO), metadata$PublicID)],
                                   pool = gsub(pattern = "pool", replacement = "", metadata$Pool[match(colnames(path_vals_KO), metadata$PublicID)]),
                                   RIN = metadata$RIN[match(colnames(path_vals_KO), metadata$PublicID)],
                                   type = metadata$type[match(colnames(path_vals_KO), metadata$PublicID)],
                                   stringsAsFactors = F)
    
    sample_group_limma <- target_matrix_KO[, c(1,2,4,5,6,7,8,9,10,12,13,14)] # Note that for Limma the Control has to be first in the levels
    design <- model.matrix( ~ 0  + type + age + gender + antimalarials + immunosuppresants + biologicals 
                            + steroids + systemic.antibiotics + pool + RIN ,  data = sample_group_limma)
    colnames(design) <- gsub(pattern = "typecontrols", "C", colnames(design)) %>% gsub("typetreated", "treated", . )
    
    rownames(design) <- sample_group_limma$samples
    cont.matrix <- makeContrasts( treatedvsC=treated-C, levels=design)
    fit <- lmFit(path_vals_KO , design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)
    
    tableLimma <- topTable(fit3, number = rownames(path_vals_KO), adjust.method ="bonferroni", sort.by="p")
    top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    
    message("Saving top deregulated circuits after Bonferroni corrections with limma")
    
    if(dim(top_LimmaALL)>0){
      
      anot_all <- data.frame(Circuit = get_path_names(pathways_db, rownames(top_LimmaALL)),
                             logFC = top_LimmaALL$logFC,
                             P.val = top_LimmaALL$P.Value,
                             FDR.Pval = top_LimmaALL$adj.P.Val,
                             UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways_db, dbannot= "uniprot" ,collapse = T),
                             GO = get_pathways_annotations(rownames(top_LimmaALL), pathways_db, dbannot= "GO" ,collapse = T), 
                             stringsAsFactors = F)
      
      write.table(anot_all, file = here("results", "KO_drugs", folder, paste0(drug, "_KO_deregCir_afterLimmaAnnot.tsv")), sep = "\t", quote = F,
                  row.names = F, col.names = T)
      
      drug_FC <- apply(path_vals_KO, 2, mean)
      
      FC <- log2(drug_FC[1:918] /drug_FC[919:length(drug_FC)])
      
      names(FC) <- gsub(paste0("_",target_name), "", names(FC))
      
      FC_df <- data.frame(patients = names(FC),
                          logFC = abs(FC),
                          ranked_patients = rank(FC),
                          Diagnosis = metadata_init$Diagnosis[match(names(FC), metadata_init$PublicID)] ,
                          Cluster =  metadata_init$DISEASE[match(names(FC), metadata_init$PublicID)] ,
                          stringsAsFactors = F) %>% .[order(.$ranked_patients, decreasing = F),]
      message("Plotting the LogFC with and without CLUSTER 2")  
      ## PLOT IT
      
      ## Calculate the median of the logFCs to divide the reponsive / non responsive groups
      mediana <- median(FC_df$logFC)
      
      FCs <- ggplot(data = FC_df,
                    mapping = aes(x = ranked_patients,
                                  y = logFC,
                                  color = Cluster)) +
        ggtitle(paste0('Log FC of the circuits activity absolute values on the ranked patients after simulation of ',drug ,' drug effect on ', target_name )) +
        geom_point(size = 3) +
        geom_hline(yintercept = mediana, linetype = "dashed") +
        theme_minimal()+
        theme(axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text = element_text(size = 18),
              legend.title =  element_text(size = 24),
              plot.title = element_text(size = 22, face = "bold"))
      
      FCs_noCluster2 <- ggplot(data = FC_df[!FC_df$Cluster== "CLUSTER2",],
                               mapping = aes(x = ranked_patients,
                                             y = logFC,
                                             color = Cluster)) +
        ggtitle(paste0('Log FC of the circuits activity absolute values on the ranked patients after simulation of ',drug ,' drug effect on ', target_name, "without cluster2" )) +
        geom_point(size = 3) +
        geom_hline(yintercept = mediana, linetype = "dashed") +
        theme_minimal()+
        theme(axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text = element_text(size = 18),
              legend.title =  element_text(size = 24),
              plot.title = element_text(size = 22, face = "bold"))
      
      
      
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
      FCdf_noC2 <- FC_df[FC_df$Cluster!= "CLUSTER2",] 
      
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
      
      
      png(filename = here("results", "KO_drugs", folder, paste0("FCplot",drug,"_KO_", target_name,"clusters.png")), width = 5000, height = 3000, res = 200)
      (FCs | bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1)) #+ plot_layout(guides = 'collect')  
      dev.off()
      
      png(filename = here("results", "KO_drugs", folder, paste0("FCplot",drug,"_KO_", target_name,"clusters_withoutCluster2.png")), width = 5000, height = 3000, res = 200)
      (FCs_noCluster2 | bars_responders_noC2/plot_spacer()+ bars_nonresponders_noC2) + plot_layout(widths = c(3, 1)) #+ plot_layout(guides = 'collect')  
      dev.off()
      
      # ggsave(last_plot(), filename = here("results", "KO_drugs", paste0("FCplot",drug,"_KO_", target_name,"clusters.png")), dpi = 400 )
      
      
    }else{
      
      message("No significant circuits detected after bonferroni p-value adjustment")
    }
    
    
  }else{
    
    message('The target of that drug could not be found in Hipathia')  
    
  }
  
  message("Saving all the plots and tables.. DONE!!")
  
  return(FC_df)
  
}


############################################################################################
#### Function to simulate the effect of a drug by KO of all its gene targets, previously 
#######  load trans_data, pathways and metadata without controls - PLOT of log2 FC 16.03.2021
###############################################

simulate_drug <- function(drug, expression_matrix = trans_data, metadata_matrix_controls = metadata_init, 
                          drug_table = drug_targetsDF, pathways_db = pathways, folder = NA){
  if(is.na(folder)){
    folder = ""
  }else{
    message(paste0("Folder where the results will be saved... ", folder))
  }
  
  possible_targets <- unique(drug_table$entrez_id[drug_table$drug_name == drug])
  
  ## Check that the target is in the pathways
  
  if(any(possible_targets %in% pathways_db$all.genes)){
    
    idx <- possible_targets %in% pathways_db$all.genes
    
    targets_name <- drug_table$gene_name[match(possible_targets[idx], drug_table$entrez_id)]
    
    message(paste0("Simulating...", drug, " with KO of ", paste(targets_name, collapse=",")," in expression dataset"))
    
    metadata_KO <- metadata_matrix_controls
    metadata_KO$PublicID <- paste0(metadata_KO$PublicID,"_" ,drug)
    metadata_KO$type <- "treated"
    
    trans_data_KO <- expression_matrix
    
    trans_data_KO[which(rownames(trans_data_KO) %in% possible_targets[idx]), ] <- 0.0001
    colnames(trans_data_KO) <- paste0(colnames(trans_data), "_", drug) 
    
    
    ## Normalize expression data all together (controls + repeated expression matrix with KO)
    trans_data_KO_C <- cbind(trans_data_KO, expression_matrix)
    exp_data_KO_C <- normalize_data(trans_data_KO_C)
    
    message("Calculating Hipathia")
    ## Using Hipathia to compute the signal (subpathways) 
    results_KO <- hipathia(exp_data_KO_C, pathways_db, decompose = FALSE, verbose=FALSE)
    saveRDS(results_KO, file = here("rds","KO_drugs", folder, paste0(drug,"_KO_", paste(targets_name, collapse="_") ,"_results.rds")))
    
    message("saving results of KO")
    
    path_vals_KO <- get_paths_data(results_KO, matrix = TRUE)
    # path_vals_phy_KO <- path_vals_KO[which(rownames(path_vals_KO) %in% subpathways_phy$hipathia),]
    # anot_cir <- data.frame(Circuit = get_path_names(pathways, rownames(path_vals_phy_KO)),
    #                        UniprotKB = get_pathways_annotations(rownames(path_vals_phy_KO), pathways, dbannot= "uniprot" ,collapse = T),
    #                        GO = get_pathways_annotations(rownames(path_vals_phy_KO), pathways, dbannot= "GO" ,collapse = T), 
    #                        stringsAsFactors = F)
    
    
    ## LIMMA FOR ALL SAMPLES 
    message("Building up Limma models")
    metadata <- rbind(metadata_KO, metadata_matrix_controls)
    target_matrix_KO <- data.frame(samples = colnames(path_vals_KO), 
                                   cluster = metadata$DISEASE[match(colnames(path_vals_KO), metadata$PublicID)],
                                   diagnosis = metadata$Diagnosis[match(colnames(path_vals_KO), metadata$PublicID)],
                                   age = metadata$Age[match(colnames(path_vals_KO), metadata$PublicID)],
                                   gender = metadata$Gender[match(colnames(path_vals_KO), metadata$PublicID)],
                                   antimalarials = metadata$Antimalarials[match(colnames(path_vals_KO), metadata$PublicID)],
                                   immunosuppresants = metadata$Immunosuppresants[match(colnames(path_vals_KO), metadata$PublicID)],
                                   biologicals = metadata$Biologicals[match(colnames(path_vals_KO), metadata$PublicID)],
                                   steroids = metadata$Steroids[match(colnames(path_vals_KO), metadata$PublicID)],
                                   systemic.antibiotics = metadata$Systemic.Antibiotics[match(colnames(path_vals_KO), metadata$PublicID)],
                                   center = metadata$PublicCenter[match(colnames(path_vals_KO), metadata$PublicID)],
                                   pool = gsub(pattern = "pool", replacement = "", metadata$Pool[match(colnames(path_vals_KO), metadata$PublicID)]),
                                   RIN = metadata$RIN[match(colnames(path_vals_KO), metadata$PublicID)],
                                   type = metadata$type[match(colnames(path_vals_KO), metadata$PublicID)],
                                   stringsAsFactors = F)
    
    sample_group_limma <- target_matrix_KO[, c(1,2,4,5,6,7,8,9,10,12,13,14)] # Note that for Limma the Control has to be first in the levels
    design <- model.matrix( ~ 0  + type + age + gender + antimalarials + immunosuppresants + biologicals 
                            + steroids + systemic.antibiotics + pool + RIN ,  data = sample_group_limma)
    colnames(design) <- gsub(pattern = "typecontrols", "C", colnames(design)) %>% gsub("typetreated", "treated", . )
    
    rownames(design) <- sample_group_limma$samples
    cont.matrix <- makeContrasts( treatedvsC=treated-C, levels=design)
    fit <- lmFit(path_vals_KO , design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)
    
    tableLimma <- topTable(fit3, number = rownames(path_vals_KO), adjust.method ="bonferroni", sort.by="p")
    top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    
    message("Saving top deregulated circuits after Bonferroni corrections with limma")
    
    if(dim(top_LimmaALL)[1]>0){
      
      anot_all <- data.frame(Circuit = get_path_names(pathways_db, rownames(top_LimmaALL)),
                             logFC = top_LimmaALL$logFC,
                             P.val = top_LimmaALL$P.Value,
                             FDR.Pval = top_LimmaALL$adj.P.Val,
                             UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways_db, dbannot= "uniprot" ,collapse = T),
                             GO = get_pathways_annotations(rownames(top_LimmaALL), pathways_db, dbannot= "GO" ,collapse = T), 
                             stringsAsFactors = F)
      
      write.table(anot_all, file = here("results", "KO_drugs", folder, paste0(drug, "_KO_deregCir_afterLimmaAnnot.tsv")), sep = "\t", quote = F,
                  row.names = F, col.names = T)
      
    }else{
      
      message("No significant circuits detected after bonferroni p-value adjustment")
      
    }
    
    
    drug_FC <- apply(path_vals_KO, 2, mean)
    
    FC <- log2(drug_FC[1:918] /drug_FC[919:length(drug_FC)])
    
    names(FC) <- gsub(paste0("_",drug), "", names(FC))
    
    FC_df <- data.frame(patients = names(FC),
                        logFC = abs(FC),
                        ranked_patients = rank(FC),
                        Diagnosis = metadata_init$Diagnosis[match(names(FC), metadata_init$PublicID)] ,
                        Cluster =  metadata_init$DISEASE[match(names(FC), metadata_init$PublicID)] ,
                        stringsAsFactors = F) %>% .[order(.$ranked_patients, decreasing = F),]
    
    message("Plotting the LogFC ")  
    ## PLOT IT
    
    ## Calculate the median of the logFCs to divide the reponsive / non responsive groups
    mediana <- median(FC_df$logFC)
    
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
    
    FCs_noCluster2 <- ggplot(data = FC_df[!FC_df$Cluster== "CLUSTER2",],
                             mapping = aes(x = ranked_patients,
                                           y = logFC,
                                           color = Cluster)) +
      ggtitle(paste0('Log FC of the circuits activity absolute values on the ranked patients after simulation of ',drug ,' drug effect without cluster2' )) +
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
    FCdf_noC2 <- FC_df[FC_df$Cluster!= "CLUSTER2",] 
    
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
    
    
    png(filename = here("results", "KO_drugs", folder, paste0("FCplot",drug,"_KO_clusters.png")), width = 5000, height = 3000, res = 200)
    print((FCs | bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1))) 
    dev.off()
    
    png(filename = here("results", "KO_drugs", folder,paste0("FCplot",drug,"_KO_withoutCluster2.png")), width = 5000, height = 3000, res = 200)
    print((FCs_noCluster2 | bars_responders_noC2/plot_spacer()+ bars_nonresponders_noC2) + plot_layout(widths = c(3, 1))) #+ plot_layout(guides = 'collect')  
    dev.off()
    
    return(FC_df)
    
    message("Saving all the plots and tables.. DONE!!")
    
  }else{
    
    message(paste0('The target of ', drug, ' could not be found in Hipathia'))
    
  }
  
}



############################################################################################
#### Function to simulate the effect of a drug by KO of all its gene targets, previously 
#####  load trans_data, pathways and metadata without controls - PLOT of log2 FC
#####  for the new dataset validation cohort 16.03.2021
###############################################
### Edited for the euclidean distance instead of the logFC adding pathvals corrections by circuits length 18.04.2021


simulate_drug_validation_cluster <- function(drug, expression_matrix = trans_data, metadata_matrix_controls = metadata_init, 
                          drug_table = drug_targetsDF, pathways_db = pathways , folder = NA, folder_res , subset_paths = pathways_phy ){
  if(is.na(folder)){
    folder = ""
  }else{
    message(paste0("Folder where the results will be saved... ", folder, " ..located in ./rds and ./results"))
  }
  
  possible_targets <- unique(drug_table$entrez_id[drug_table$drug_name == drug])
  
  ## Check that the target is in the pathways
  
  if(any(possible_targets %in% pathways_db$all.genes)){
    
    idx <- possible_targets %in% pathways_db$all.genes
    
    targets_name <- drug_table$gene_name[match(possible_targets[idx], drug_table$entrez_id)]
    
    message(paste0("Simulating...", drug, " with KO of ", paste(targets_name, collapse=",")," in expression dataset"))
    
    metadata_KO <- metadata_matrix_controls
    metadata_KO$samples <- paste0(metadata_KO$samples,"_" ,drug)
    metadata_KO$type <- "treated"
    
    if(!file.exists(file = here("rds","KO_drugs", folder_res, paste0(drug,"_KO_", paste(targets_name, collapse="_") ,"_validation_results.rds")))){
      
      trans_data_KO <- expression_matrix

      trans_data_KO[which(rownames(trans_data_KO) %in% possible_targets[idx]), ] <- 0.0001
      colnames(trans_data_KO) <- paste0(colnames(trans_data), "_", drug)


      ## Normalize expression data all together (controls + repeated expression matrix with KO)
      trans_data_KO_C <- cbind(trans_data_KO, expression_matrix)
      exp_data_KO_C <- normalize_data(trans_data_KO_C)

      message("Calculating Hipathia")
      ## Using Hipathia to compute the signal (subpathways)
      results_KO <- hipathia(exp_data_KO_C, pathways_db, decompose = FALSE, verbose=FALSE)
      saveRDS(results_KO, file = here("rds","KO_drugs", folder_res, paste0(drug,"_KO_", paste(targets_name, collapse="_") ,"_validation_results.rds")))
      
        } else {
      
        results_KO <- readRDS(file = here("rds","KO_drugs", folder_res, paste0(drug,"_KO_", paste(targets_name, collapse="_") ,"_validation_results.rds")))
        
        # message("saving results of KO")
        message("reading older results of KO.. DONE")
        
        }
    
    path_vals_KO <- get_paths_data(results_KO, matrix = TRUE)
    path_vals_KO <- normalize_paths(path_vals_KO, pathways_db)
    
    ### SUBSET ONLY PHYSIOLOGICAL PATHS

    pathways_phy <- subset_paths
    
    subpathways.list <- lapply(pathways_phy$pathigraphs, function(x){names(x$effector.subgraphs)})
    
    subpathways_phy <- data.frame(unlist(subpathways.list), stringsAsFactors = F)
    
    colnames(subpathways_phy) <- "hipathia"
    
    path_vals_phy_KO <- path_vals_KO[which(rownames(path_vals_KO) %in% subpathways_phy$hipathia),]
  
    
    ## LIMMA FOR ALL SAMPLES 
    message("Building up Limma models")
    metadata <- rbind(metadata_KO, metadata_matrix_controls)
    target_matrix_KO <- data.frame(samples = colnames(path_vals_phy_KO), 
                                   cluster = metadata$cluster[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   diagnosis = metadata$diagnosis[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   age = metadata$age[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   gender = metadata$gender[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   center = metadata$center[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   pool = gsub(pattern = "pool", replacement = "", metadata$pool[match(colnames(path_vals_phy_KO), metadata$samples)]),
                                   RIN = metadata$RIN[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   type = metadata$type[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   stringsAsFactors = F)
    
    sample_group_limma <- target_matrix_KO[, c(1,4,5,6,7,8,9)] # Note that for Limma the Control has to be first in the levels
    design <- model.matrix( ~ 0  + type + age + gender + center + pool + RIN ,  data = sample_group_limma)
    colnames(design) <- gsub(pattern = "typecontrols", "C", colnames(design)) %>% gsub("typetreated", "treated", . )
    
    rownames(design) <- sample_group_limma$samples
    cont.matrix <- makeContrasts( treatedvsC=treated-C, levels=design)
    fit <- lmFit(path_vals_phy_KO , design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)
    
    tableLimma <- topTable(fit3, number = rownames(path_vals_phy_KO), adjust.method ="bonferroni", sort.by="p")
    top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    
    message("Saving top deregulated circuits after Bonferroni corrections with limma")
    
    if(dim(top_LimmaALL)[1]>0){
      
      anot_all <- data.frame(Circuit = get_path_names(subset_paths, rownames(top_LimmaALL)),
                             logFC = top_LimmaALL$logFC,
                             P.val = top_LimmaALL$P.Value,
                             Bonferroni.adj.Pval = top_LimmaALL$adj.P.Val,
                             UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), subset_paths, dbannot= "uniprot" ,collapse = T),
                             GO = get_pathways_annotations(rownames(top_LimmaALL), subset_paths, dbannot= "GO" ,collapse = T), 
                             stringsAsFactors = F)
      
      write.table(anot_all, file = here("results", "KO_drugs", folder, paste0(drug, "_KO_deregCir_validation_afterLimmaAnnot_cluster.tsv")), sep = "\t", quote = F,
                  row.names = F, col.names = T)
      
    }else{
      
      message("No significant circuits detected after bonferroni p-value adjustment")
      
    }
    
    
    # drug_FC <- apply(path_vals_KO, 2, mean)
    # FC <- log2(drug_FC[1:1441] /drug_FC[1442:length(drug_FC)])
     
    # WE WILL CALCULATE THE EUCLIDEAN DISTANCE
    
    disease_pvals <- path_vals_phy_KO[,1:1441]
    controls_pvals <- path_vals_phy_KO[,1442:dim(path_vals_phy_KO)[2]]
    
    
    FC <- (disease_pvals - controls_pvals)**2 %>% apply(., 2, sum) %>% sqrt(.)

    names(FC) <- gsub(paste0("_",drug), "", names(FC))
    
    FC_df <- data.frame(patients = names(FC),
                        Euclidean_Distance = FC,
                        Ranked_Patients = rank(FC),
                        Diagnosis = metadata_init$diagnosis[match(names(FC), metadata_init$samples)] ,
                        Cluster =  as.character(metadata_init$cluster[match(names(FC), metadata_init$samples)]) ,
                        stringsAsFactors = F) %>% .[order(.$Ranked_Patients, decreasing = F),]
    
    message("Plotting the Euclidean distances ")  
    ## PLOT IT
    
    ## Calculate the median of the Euclidean distance to divide the reponsive / non responsive groups
    mediana <- median(FC_df$Euclidean_Distance)
    
    FCs <- ggplot(data = FC_df,
                  mapping = aes(x = Ranked_Patients,
                                y = Euclidean_Distance,
                                color = Cluster)) +
      ggtitle(paste0('Euclidean distance of the circuits activity absolute values \non the ranked patients after simulation of ',drug ,' drug effect' )) +
      geom_point(size = 3) +
      geom_hline(yintercept = mediana, linetype = "dashed") +
      theme_minimal()+
      theme(axis.title.x = element_text(size = 28),
            axis.title.y = element_text(size = 28),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.text = element_text(size = 28),
            legend.title =  element_text(size = 32),
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
    
    FCs_noCluster2 <- ggplot(data = FC_df[!FC_df$Cluster== "2",],
                             mapping = aes(x = Ranked_Patients,
                                           y = Euclidean_Distance,
                                           color = Cluster)) +
      ggtitle(paste0('Euclidean distance of the circuits activity absolute values \non the ranked patients after simulation of ',drug ,' drug effect without Cluster 2' )) +
      geom_point(size = 3) +
      geom_hline(yintercept = mediana, linetype = "dashed") +
      theme_minimal()+
      theme(axis.title.x = element_text(size = 28),
            axis.title.y = element_text(size = 28),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.text = element_text(size = 28),
            legend.title =  element_text(size = 32),
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
    
    ## With CLUSTER2 
    bars_df_res <- data.frame(table(FC_df$Cluster[FC_df$Euclidean_Distance >= mediana]))
    
    bars_df_NOres <- data.frame(table(FC_df$Cluster[FC_df$Euclidean_Distance < mediana]))
    
    bars_responders <- ggplot(data = bars_df_res, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
      ggtitle("High response patients")+
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            axis.text.y = element_text(face="bold",size= 26),
            axis.text.x = element_text(face="bold",size= 20),
            legend.position = "none",
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5))+
      scale_y_continuous(limits = c(0, max(bars_df_NOres$Freq, bars_df_res$Freq))) +
      coord_flip()
    
    
    bars_nonresponders <- ggplot(data = bars_df_NOres, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
      ggtitle("Low response patients")+
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            axis.text.y = element_text(face="bold",size= 26),
            axis.text.x = element_text(face="bold",size= 20),
            legend.position = "none",
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5)) +
      scale_y_continuous(limits = c(0, max(bars_df_NOres$Freq, bars_df_res$Freq))) +
      coord_flip()
    
    
    ## NO CLUSTER2 
    FCdf_noC2 <- FC_df[FC_df$Cluster!= "2",] 
    
    bars_df_res_noCluster2 <- data.frame(table(FCdf_noC2$Cluster[FCdf_noC2$Euclidean_Distance >= mediana]))
    bars_df_NOres_noCluster2 <- data.frame(table(FCdf_noC2$Cluster[FCdf_noC2$Euclidean_Distance < mediana]))
    
    bars_responders_noC2 <- ggplot(data = bars_df_res_noCluster2, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
      ggtitle("High response patients")+
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            axis.text.y = element_text(face="bold",size= 26),
            axis.text.x = element_text(face="bold",size= 20),
            legend.position = "none",
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5))+
      scale_y_continuous(limits = c(0, max(bars_df_NOres_noCluster2$Freq, bars_df_res_noCluster2$Freq)))+
      coord_flip()
    
    
    bars_nonresponders_noC2 <- ggplot(data = bars_df_NOres_noCluster2, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
      ggtitle("Low response patients")+
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            axis.text.y = element_text(face="bold",size= 26),
            axis.text.x = element_text(face="bold",size= 20),
            legend.position = "none",
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5)) +
      scale_y_continuous(limits = c(0, max(bars_df_NOres_noCluster2$Freq, bars_df_res_noCluster2$Freq)))+
      coord_flip()
    
    
    png(filename = here("results", "KO_drugs", folder, paste0("EuDistplot_validation",drug,"_KO_clusters.png")), width = 5000, height = 3000, res = 200)
    print((FCs | bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1))) 
    dev.off()
    
    png(filename = here("results", "KO_drugs", folder,paste0("EuDistplot_validation",drug,"_KO_withoutCluster2.png")), width = 5000, height = 3000, res = 200)
    print((FCs_noCluster2 | bars_responders_noC2/plot_spacer()+ bars_nonresponders_noC2) + plot_layout(widths = c(3, 1))) #+ plot_layout(guides = 'collect')  
    dev.off()
    
    return(FC_df)
    
    message("Saving all the plots and tables.. DONE!!")
    
  }else{
    
    message(paste0('The target of ', drug, ' could not be found in Hipathia'))
    
  }
  
}

############################################################################################
#### Function to  simulated the effect of a drug by KO of all its gene targets,
####  same as function simulate_drug_validation(), and
#### plot the logFC painting the diagnosis, instead of the  CLUSTER. 13.04.2021
###############################################

simulate_drug_validation_diagnosis <- function(drug, expression_matrix = trans_data, metadata_matrix_controls = metadata_init, 
                                               drug_table = drug_targetsDF, pathways_db = pathways, folder = NA, folder_res = NA, subset_paths = pathways_phy){
  if(is.na(folder)){
    folder = ""
  }else{
    message(paste0("Folder where the results will be saved... ", folder, " ..located in ./rds and ./results"))
  }
  
  metadata_init <- metadata_matrix_controls
  
  possible_targets <- unique(drug_table$entrez_id[drug_table$drug_name == drug])
  
  ## Check that the target is in the pathways
  
  if(any(possible_targets %in% pathways_db$all.genes)){
    
    idx <- possible_targets %in% pathways_db$all.genes
    
    targets_name <- drug_table$gene_name[match(possible_targets[idx], drug_table$entrez_id)]
    
    message(paste0("Simulating...", drug, " with KO of ", paste(targets_name, collapse=",")," in expression dataset"))
    
    metadata_KO <- metadata_matrix_controls
    metadata_KO$samples <- paste0(metadata_KO$samples,"_" ,drug)
    metadata_KO$type <- "treated"
    
    if(!file.exists(file = here("rds","KO_drugs", folder_res, paste0(drug,"_KO_", paste(targets_name, collapse="_") ,"_validation_results.rds")))){
      
      trans_data_KO <- expression_matrix
      
      trans_data_KO[which(rownames(trans_data_KO) %in% possible_targets[idx]), ] <- 0.0001
      colnames(trans_data_KO) <- paste0(colnames(trans_data), "_", drug)
      
      
      ## Normalize expression data all together (controls + repeated expression matrix with KO)
      trans_data_KO_C <- cbind(trans_data_KO, expression_matrix)
      exp_data_KO_C <- normalize_data(trans_data_KO_C)
      
      message("Calculating Hipathia")
      ## Using Hipathia to compute the signal (subpathways)
      results_KO <- hipathia(exp_data_KO_C, pathways_db, decompose = FALSE, verbose=FALSE)
      saveRDS(results_KO, file = here("rds","KO_drugs", folder_res, paste0(drug,"_KO_", paste(targets_name, collapse="_") ,"_validation_results.rds")))
      
    } else {
      
      results_KO <- readRDS(file = here("rds","KO_drugs", folder_res, paste0(drug,"_KO_", paste(targets_name, collapse="_") ,"_validation_results.rds")))
      
      # message("saving results of KO")
      message("reading older results of KO.. DONE")
      
    }
    
    path_vals_KO <- get_paths_data(results_KO, matrix = TRUE)
    path_vals_KO <- normalize_paths(path_vals_KO, pathways_db)
    
    ### SUBSET ONLY PHYSIOLOGICAL PATHS
    
    pathways_phy <- subset_paths
    
    subpathways.list <- lapply(pathways_phy$pathigraphs, function(x){names(x$effector.subgraphs)})
    
    subpathways_phy <- data.frame(unlist(subpathways.list), stringsAsFactors = F)
    
    colnames(subpathways_phy) <- "hipathia"
    
    path_vals_phy_KO <- path_vals_KO[which(rownames(path_vals_KO) %in% subpathways_phy$hipathia),]
    
    ## LIMMA FOR ALL SAMPLES 
    message("Building up Limma models")
    metadata <- rbind(metadata_KO, metadata_matrix_controls)
    target_matrix_KO <- data.frame(samples = colnames(path_vals_phy_KO), 
                                   cluster = metadata$cluster[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   diagnosis = metadata$diagnosis[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   age = metadata$age[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   gender = metadata$gender[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   center = metadata$center[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   pool = gsub(pattern = "pool", replacement = "", metadata$pool[match(colnames(path_vals_phy_KO), metadata$samples)]),
                                   RIN = metadata$RIN[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   type = metadata$type[match(colnames(path_vals_phy_KO), metadata$samples)],
                                   stringsAsFactors = F)
    
    sample_group_limma <- target_matrix_KO[, c(1,4,5,6,7,8,9)] # Note that for Limma the Control has to be first in the levels
    design <- model.matrix( ~ 0  + type + age + gender + center + pool + RIN ,  data = sample_group_limma)
    colnames(design) <- gsub(pattern = "typecontrols", "C", colnames(design)) %>% gsub("typetreated", "treated", . )
    
    rownames(design) <- sample_group_limma$samples
    cont.matrix <- makeContrasts( treatedvsC=treated-C, levels=design)
    fit <- lmFit(path_vals_phy_KO , design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)
    
    tableLimma <- topTable(fit3, number = rownames(path_vals_phy_KO), adjust.method ="bonferroni", sort.by="p")
    top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    
    message("Saving top deregulated circuits after Bonferroni corrections with limma")
    
    if(dim(top_LimmaALL)[1]>0){
      
      anot_all <- data.frame(Circuit = get_path_names(subset_paths, rownames(top_LimmaALL)),
                             logFC = top_LimmaALL$logFC,
                             P.val = top_LimmaALL$P.Value,
                             Bonferroni.adj.Pval = top_LimmaALL$adj.P.Val,
                             UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), subset_paths, dbannot= "uniprot" ,collapse = T),
                             GO = get_pathways_annotations(rownames(top_LimmaALL), subset_paths, dbannot= "GO" ,collapse = T), 
                             stringsAsFactors = F)
      
      write.table(anot_all, file = here("results", "KO_drugs", folder, paste0(drug, "_KO_deregCir_validation_afterLimmaAnnot_diagnosis.tsv")), sep = "\t", quote = F,
                  row.names = F, col.names = T)
      
    }else{
      
      message("No significant circuits detected after bonferroni p-value adjustment")
      
    }
    
    
    # drug_FC <- apply(path_vals_KO, 2, mean)
    # FC <- log2(drug_FC[1:1441] /drug_FC[1442:length(drug_FC)])ç
    
    # WE WILL CALCULATE THE EUCLIDEAN DISTANCE
    
    disease_pvals <- path_vals_phy_KO[,1:1441]
    controls_pvals <- path_vals_phy_KO[,1442:dim(path_vals_phy_KO)[2]]
    
    
    FC <- (disease_pvals - controls_pvals)**2 %>% apply(., 2, sum) %>% sqrt(.)
    
    names(FC) <- gsub(paste0("_",drug), "", names(FC))
    
    FC_df <- data.frame(patients = names(FC),
                        Euclidean_Distance = FC,
                        Ranked_Patients = rank(FC),
                        Diagnosis = metadata_init$diagnosis[match(names(FC), metadata_init$samples)] ,
                        Cluster =  as.character(metadata_init$cluster[match(names(FC), metadata_init$samples)]) ,
                        stringsAsFactors = F) %>% .[order(.$Ranked_Patients, decreasing = F),]
    
    message("Plotting the Euclidean distances ")  
    ## PLOT IT
    
    ## Calculate the median of the Euclidean distance to divide the reponsive / non responsive groups
    mediana <- median(FC_df$Euclidean_Distance)
    
    FCs <- ggplot(data = FC_df,
                  mapping = aes(x = Ranked_Patients,
                                y = Euclidean_Distance,
                                color = Diagnosis)) +
      ggtitle(paste0('Euclidean distance of the circuits activity absolute values \non the ranked patients after simulation of ',drug ,' drug effect' )) +
      geom_point(size = 3) +
      geom_hline(yintercept = mediana, linetype = "dashed") +
      theme_minimal()+
      theme(axis.title.x = element_text(size = 28),
            axis.title.y = element_text(size = 28),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.text = element_text(size = 28),
            legend.title =  element_text(size = 32),
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
    
    bars_df_res <- data.frame(table(FC_df$Diagnosis[FC_df$Euclidean_Distance >= mediana]))
    
    bars_df_NOres <- data.frame(table(FC_df$Diagnosis[FC_df$Euclidean_Distance < mediana]))
    
    bars_responders <- ggplot(data = bars_df_res, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
      ggtitle("High response patients")+
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            axis.text.y = element_text(face="bold",size= 26),
            axis.text.x = element_text(face="bold",size= 20),
            legend.position = "none",
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5))+
      scale_y_continuous(limits = c(0, max(bars_df_NOres$Freq, bars_df_res$Freq))) +
      coord_flip()
    
    
    bars_nonresponders <- ggplot(data = bars_df_NOres, mapping = aes(x = Var1, y = Freq, fill = Var1)) +
      ggtitle("Low response patients")+
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.title = element_blank(),
            axis.text.y = element_text(face="bold",size= 26),
            axis.text.x = element_text(face="bold",size= 20),
            legend.position = "none",
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5)) +
      scale_y_continuous(limits = c(0, max(bars_df_NOres$Freq, bars_df_res$Freq))) +
      coord_flip()
    
    png(filename = here("results", "KO_drugs", folder, paste0("EuDistplot_validation",drug,"_KO_diagnosis.png")), width = 5000, height = 3000, res = 200)
    print((FCs | bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1))) 
    dev.off()

    return(FC_df)
    
    message("Saving all the plots and tables.. DONE!!")
    
  }else{
    
    message(paste0('The target of ', drug, ' could not be found in Hipathia'))
    
  }
  
}



############################################################################################
#### Function to do histogram plot of the X top changing circuits (after KO_effect) after measured by a t-test, comparing 
#### drug_KO effect vs Controls vs Diseases  14.04.2021
###############################################


plot_hist_drugEffvsC_topX <- function(pathvals_controls, pathvals_disease, pathvals_drugEff, folder , metadata_table = metadata ,top_LimmaALL_table = top_LimmaALL, cluster2 = T,
                                       totalDeregcir_plot = 100 ){
  
  require(ggplot2)
  
  ## We extract drugs effect.
  drug_name <- unique(str_split(colnames(pathvals_drugEff), "_")[[1]][2])
  
  ## Clean up the names from the drug for later filtering samples.
  colnames(pathvals_drugEff) <- gsub(paste0("_", drug_name), "", colnames(pathvals_drugEff))
  
  
  ## Select the number of top deregulated circuits from top_LimmaALL, which are the top deregulated circuits ordered by Bonferroni corrected p-values 
  ## from the comparison between all disease samples versus Controls.
  if(is.numeric(totalDeregcir_plot)){
    
    message(paste0("Calculating histograms comparisons plots for ", drug_name, " in top ", totalDeregcir_plot, " deregulated circuits"))
    
  }else{
    message("Incorrect value for the attribute totalDeregcir_plot, the value must be a positive integer")
  }
  
  ## Read table of top deregulated circuits 
  
  top_dereg <- rownames(top_LimmaALL_table)[1:totalDeregcir_plot]
  # top_dereg <- rownames(top_LimmaALL_table)
  
  ## Do a t-test to check which differences are relevant
  
  vals_dis <- pathvals_disease[top_dereg,]
  vals_drugeff <- pathvals_drugEff[top_dereg,]
  
  diff <- data.frame(cir = rownames(vals_dis), 
                     sig = sapply(1:nrow(vals_dis), function(x) t.test(vals_dis[x,], vals_drugeff[x,])$p.value), stringsAsFactors = F) 
  
  diff <- diff[diff$sig < 0.05, ] %>%  .[order(.$sig, decreasing = F),] 
  diff$cir_names <- get_path_names(pathways, as.character(diff$cir))
  
  write.table(diff, file = here("results", "KO_drugs", folder, paste0(drug_name,"Ttest_significantCir.tsv")), 
              row.names = F, col.names = T, quote = F, sep = "\t")                              
  
  ## Ordering most changinf circuits after drug effect based in the differences in the median (NOT THE BEST, BETTER DOING A t-test)
  # diff_dis_drug <- data.frame(cir = rownames(pathvals_disease),diff = abs(apply(pathvals_disease, 1, median) - apply(pathvals_drugEff, 1, median)), stringsAsFactors = F)
  # diff_dis_drug <- diff_dis_drug[diff_dis_drug$cir %in% top_dereg, ]
  # diff_dis_drug <- diff_dis_drug[order(diff_dis_drug$diff, decreasing = F),]
  # diff_dis_drug <- diff_dis_drug[diff_dis_drug$diff > 0.0001,]
  # diff_dis_drug$cir_names <- get_path_names(pathways, as.character(diff_dis_drug$cir))
  
  
  ## Read de controls and disease data and reshape to stack it!!
  idxC <- which(rownames(pathvals_controls) %in% diff$cir)
  dat_c <- stack(as.data.frame(t(pathvals_controls[idxC,])))
  dat_c <- dat_c[order(dat_c$values, decreasing = F ),]
  dat_c$group <- "Control"
  
  
  if(cluster2 == F){
    
    idx <- which(rownames(pathvals_drugEff) %in%  diff$cir)
    dat_drug <- pathvals_drugEff[idx, colnames(pathvals_drugEff) %in% metadata_table$samples[metadata_table$cluster != 2]]
    dat_drug <- stack(as.data.frame(t(dat_drug)))
    dat_drug$group <- "Drug-treated"
    
    idxD <- which(rownames(pathvals_disease) %in%  diff$cir)
    dat_d <- stack(as.data.frame(t(pathvals_disease[idxD, colnames(pathvals_disease) %in% metadata_table$samples[metadata_table$cluster != 2]])))
    dat_d$group <- "Disease"
    
    message("NOT INCLUDING CLUSTER 2") 
    
    title_plot <- "no_Cluster2"
    
    C2 <- "without Cluster 2"
    
  } else {
    
    idx <- which(rownames(pathvals_drugEff) %in% diff$cir )
    dat_drug <- pathvals_drugEff[idx, ]
    dat_drug <- stack(as.data.frame(t(dat_drug)))
    dat_drug$group <- "Drug-treated"
    
    idxD <- which(rownames(pathvals_disease) %in% diff$cir )
    dat_d <- stack(as.data.frame(t(pathvals_disease[idxD,])))
    dat_d$group <- "Disease"
    
    message("INCLUDING CLUSTER 2") 
    
    title_plot <- "with_Cluster2"
    
    C2 <- ""
    
  }
  
  
  dat <- rbind(dat_c, dat_d, dat_drug)
  colnames(dat) <- c("activation_values", "circuits", "Group") 
  dat_plot <- dat
  dat_plot$circuits <- get_path_names(pathways, as.character(dat_plot$circuits))
  dat_plot$circuits <- factor(dat_plot$circuits, levels = diff$cir_names, ordered = T)
  
  
  
  p <- ggplot(dat_plot, aes(x = circuits, y = activation_values, fill = Group)) + 
    ggtitle(paste0('Distribution of the circuits activity absolute values \n from the most affected circuits after simulation of ',drug_name ,' effect ', C2)) +
    geom_boxplot(outlier.fill = "black", outlier.size = 0.025, varwidth=T, lwd = 0.05 ) +
    # theme_minimal() +
    theme( axis.title.y = element_blank(),
           panel.grid = element_line(),
           axis.title.x = element_blank(),
           axis.text.y = element_text(size= 18),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0, size = 4),
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 18),
           plot.title = element_text(size = 22, face = "bold", hjust = 0.5))
  
  
  # png(filename = here("results", "KO_drugs", folder, paste0(drug_name, title_plot,"vsC_boxplots.png")), width = 12000, height = 7000, res = 600)
  ggsave(filename = here("results", "KO_drugs", folder, paste0(drug_name, title_plot,"vsC_boxplots.png")), plot = p, dpi = 400, width = 18, height = 8)
  # dev.new()
  # p    
  # dev.off()
  
  return(dat_plot)
  
}




############################################################################################
#### Function to get a resume table and do boxplots after drugs simulations on Disease samples of the activity values of the circuits that changed towards the Controls activity values 
#### from the top deregulated circuits after measured by a t-test, comparing drug_KO effect vs Controls vs Diseases  05.08.2021
###############################################

plot_boxplots_drugEffvsCvsD_topDeregulatedCir_withoppositeDirection <- function(pathvals_controls, pathvals_disease, pathvals_drugEff, folder , 
                                                                                metadata_table = metadata ,top_LimmaALL_table = top_LimmaALL, cluster2 = T){
  
  require(ggplot2)
  
  ## Extract the drug effect and the name of the drug
  drug_name <- unique(str_split(colnames(pathvals_drugEff), "_")[[1]][2])
  
  ## Clean up the names from the drug for later filtering samples.
  colnames(pathvals_drugEff) <- gsub(paste0("_", drug_name), "", colnames(pathvals_drugEff))
  
  
  ## Read table of top deregulated circuits 
  
  top_dereg <- rownames(top_LimmaALL_table)#[1:totalDeregcir_plot]
  
  
  ## Do a t-test to check which differences are relevant
  
  vals_dis <- pathvals_disease[top_dereg,]
  vals_drugeff <- pathvals_drugEff[top_dereg,]
  
  diff <- data.frame(cir = rownames(vals_dis), 
                     Boferroni.pval.Ttest = sapply(1:nrow(vals_dis), function(x) t.test(vals_dis[x,], vals_drugeff[x,])$p.value), 
                     Bonferroni.pval.DvsC = top_LimmaALL_table$adj.P.Val[match(rownames(vals_dis), rownames(top_LimmaALL_table))] ,
                     stringsAsFactors = F)
  
  diff$Boferroni.pval.Ttest <- stats::p.adjust(diff$Boferroni.pval.Ttest, method = "bonferroni" )
  
  diff <- diff[diff$Boferroni.pval.Ttest < 0.05, ] %>%  .[order(.$Boferroni.pval.Ttest, decreasing = F),] 
  diff$cir_names <- get_path_names(pathways, as.character(diff$cir))
  diff$median_C <- apply(pathvals_controls[diff$cir,], 1, median)
  diff$median_D <- apply(pathvals_disease[diff$cir,], 1, median)
  diff$median_T <- apply(pathvals_drugEff[diff$cir,], 1, median)
  
  
  for (i in 1: dim(diff)[1]){
    
    if(log2(diff$median_D[i]) - log2(diff$median_C[i])>=0){
      
      diff$FC.DvsC <- log2(diff$median_D[i]) - log2(diff$median_C[i])
      diff$FC.TvsD <- log2(diff$median_T[i]) - log2(diff$median_D[i])
      diff$disease_regulation[i] <- "UP-REGULATION"
      
      
    }else{
      
      diff$FC.DvsC <- log2(diff$median_D[i]) - log2(diff$median_C[i])
      diff$FC.TvsD <- log2(diff$median_T[i]) - log2(diff$median_D[i])
      diff$disease_regulation[i] <- "DOWN-REGULATION"
      
    }
    
    
  }
  
  
  for (i in 1: dim(diff)[1]){
    
    if(log2(diff$median_D[i]) - log2(diff$median_T[i])>=0){
      
      diff$treatment_direction[i] <- "DOWN-REGULATION"
      
    }else{
      
      diff$treatment_direction[i] <- "UP-REGULATION"
      
    }
    
    
  }
  
  
  ### FILTER OUT ONLY THOSE THAT GO IN THE CORRECT DIRECTION
  
  diff <- diff[which(diff$disease_regulation != diff$treatment_direction), ]
  
  
  write.table(diff, file = here("results", "KO_drugs", folder, paste0(drug_name,"Ttest_treatment_direction.tsv")), 
              row.names = F, col.names = T, quote = F, sep = "\t")   
  
  ## Read de controls and disease data and reshape to stack it!!
  idxC <- which(rownames(pathvals_controls) %in% diff$cir)
  dat_c <- stack(as.data.frame(t(pathvals_controls[idxC,])))
  dat_c <- dat_c[order(dat_c$values, decreasing = F ),]
  dat_c$group <- "Control"
  
  
  if(cluster2 == F){
    
    idx <- which(rownames(pathvals_drugEff) %in%  diff$cir)
    dat_drug <- pathvals_drugEff[idx, colnames(pathvals_drugEff) %in% metadata_table$samples[metadata_table$cluster != 2]]
    dat_drug <- stack(as.data.frame(t(dat_drug)))
    dat_drug$group <- "Drug-treated"
    
    idxD <- which(rownames(pathvals_disease) %in%  diff$cir)
    dat_d <- stack(as.data.frame(t(pathvals_disease[idxD, colnames(pathvals_disease) %in% metadata_table$samples[metadata_table$cluster != 2]])))
    dat_d$group <- "Disease"
    
    message("NOT INCLUDING CLUSTER 2") 
    
    title_plot <- "no_Cluster2"
    
    C2 <- "without Cluster 2"
    
  } else {
    
    idx <- which(rownames(pathvals_drugEff) %in% diff$cir )
    dat_drug <- pathvals_drugEff[idx, ]
    dat_drug <- stack(as.data.frame(t(dat_drug)))
    dat_drug$group <- "Drug-treated"
    
    idxD <- which(rownames(pathvals_disease) %in% diff$cir )
    dat_d <- stack(as.data.frame(t(pathvals_disease[idxD,])))
    dat_d$group <- "Disease"
    
    message("INCLUDING CLUSTER 2") 
    
    title_plot <- "with_Cluster2"
    
    C2 <- ""
    
  }
  
  
  dat <- rbind(dat_d, dat_c, dat_drug)
  colnames(dat) <- c("activation_values", "circuits", "Group") 
  dat_plot <- dat
  dat_plot <- dat_plot[order(dat_plot$activation_values, decreasing = F), ]
  dat_plot$circuits <- get_path_names(pathways, as.character(dat_plot$circuits))
  # dat_plot$circuits <- factor(dat_plot$circuits, levels =  diff$cir_names[match(top_dereg, diff$cir)], ordered = T) ## Order by higer adj.pvalue
  dat_plot$circuits <- factor(dat_plot$circuits, levels =  unique(dat_plot$circuits), ordered = T)
  dat_plot$Group <- factor(dat_plot$Group, levels = unique(dat_plot$Group))
  
  
  p <- ggplot(dat_plot, aes(x = circuits, y = activation_values, fill = Group)) + 
    ggtitle(paste0('Distribution of the circuits activity absolute values \n from the most affected circuits after simulation of ',drug_name ,' effect ', C2)) +
    geom_boxplot(outlier.fill = "black", outlier.size = 0.025, varwidth=T, lwd = 0.05 ) +
    # theme_minimal() +
    theme( axis.title.y = element_blank(),
           panel.grid = element_line(),
           axis.title.x = element_blank(),
           axis.text.y = element_text(size= 18),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1, size = 8, face = "bold", margin = margin(0)),
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 18),
           plot.title = element_text(size = 22, face = "bold", hjust = 0.5))
  
  
  # png(filename = here("results", "KO_drugs", folder, paste0(drug_name, title_plot,"vsC_boxplots.png")), width = 12000, height = 7000, res = 600)
  ggsave(filename = here("results", "KO_drugs", folder, paste0(drug_name, title_plot,"vsC_boxplots.png")), plot = p, dpi = 400, width = 18, height = 8)
  # dev.new()
  # p    
  # dev.off()
  
  return(dat_plot)
  
}


