      ##########################################################################################
#### Calculate the effect of the Drug by using already calculated drug simulations on disease  ######
  ####  with Hipathia in patients data and now calculate eucl.distance with Controls      #####
      ######################################################################################


library(pacman)
  
devtools::install_github("thomasp85/patchwork")

pacman::p_load("here", "hipathia", "utils", "stringr", "genefilter", "tidyr","data.table","limma", "ggplot2", "patchwork", "magrittr")


### Read Hipathia results from the validation dataset (2 dataset send by Guillermo) ####

results <- readRDS(file = here("rds","validation_cohort_Hiresults.rds"))

metadata <- fread(file = here("data", "feedback_guillermo_Marzo2021", "metadata_reformatted.tsv"))

top_LimmaALL<- readRDS(here("rds", "top_LimmaALL_validationCohort_allvsC.rds"))

## Loading Pathways (all and only physiological)

path_list <- read.table(file = here("data", "physiological_paths.tsv"), sep = "\t") #physiological_pathways
pathways <- load_pathways("hsa")
pathways_phy <- load_pathways("hsa", pathways_list = path_list$V2 )

subpathways.list <- lapply(pathways_phy$pathigraphs, function(x){names(x$effector.subgraphs)})

subpathways_phy <- data.frame(unlist(subpathways.list), stringsAsFactors = F)

colnames(subpathways_phy) <- "hipathia"

path_vals <- get_paths_data(results, matrix = TRUE)

path_vals_phy <- path_vals[which(rownames(path_vals) %in% subpathways_phy$hipathia),]
path_vals_phy <- normalize_paths(path_vals_phy, pathways)

path_vals_phy_C <- path_vals_phy[, which(colnames(path_vals_phy) %in% metadata$samples[metadata$diagnosis == "CTRL"])]

path_vals_phy_D <- path_vals_phy[, which(colnames(path_vals_phy) %in% metadata$samples[metadata$diagnosis != "CTRL"])]

### Load Drugs .rds results and Hipathia validation cohort results ####

  drugs_files <-  list.files(path = here("rds/KO_drugs/on_2dataset_bydiagnosis"), pattern = "results.rds", full.names = TRUE)
  
  drug_Hi <- lapply(drugs_files, readRDS)
  
  names(drug_Hi) <- sapply(drugs_files, function(x){str_split(str_split(x, "_")[[1]][5], "/")[[1]][2]})
  
  drug_pathvals <- lapply(drug_Hi, function(x){ path_vals <- get_paths_data(x, matrix = TRUE)
                                                path_vals_phy <- path_vals[which(rownames(path_vals) %in% subpathways_phy$hipathia), grep("_", colnames(path_vals))]                                                      
                                                path_vals_phy <- normalize_paths(path_vals_phy, pathways)
                                                 return(path_vals_phy)}) 
 
    
### RESHAPE AND PLOT COMPARISONS ####
  

  plot_hist_drugEffvsC_top40 <- function(pathvals_controls, pathvals_disease, pathvals_drugEff, folder , metadata_table = metadata ,top_LimmaALL_table = top_LimmaALL, cluster2 = T,
                                         totalDeregcir_plot = 100 ){
    
    require(ggplot2)
    
    ## Extract the drug effect and the name of the drug
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

    
    ## Do a t-test to check which differences are relevant
     
     vals_dis <- pathvals_disease[top_dereg,]
     vals_drugeff <- pathvals_drugEff[top_dereg,]
     
     diff <- data.frame(cir = rownames(vals_dis), 
                        sig = sapply(1:nrow(vals_dis), function(x) t.test(vals_dis[x,], vals_drugeff[x,])$p.value), stringsAsFactors = F) 
     
     diff <- diff[diff$sig < 0.05, ] %>%  .[order(.$sig, decreasing = F),] 
     diff$cir_names <- get_path_names(pathways, as.character(diff$cir))
      
     write.table(diff, file = here("results", "KO_drugs", folder, paste0(drug_name,"Ttest_significantCir.tsv")), 
                 row.names = F, col.names = T, quote = F, sep = "\t")                              
                                                                
    ## Ordering most changing circuits after drug effect based in the differences in the median (NOT THE BEST, BETTER DOING A t-test)
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

    
    return(dat_plot)
                 
  }

  a <- lapply(drug_pathvals, function(x){plot_hist_drugEffvsC_top40(pathvals_controls = path_vals_phy_C, pathvals_disease = path_vals_phy_D, pathvals_drugEff = x, folder = "on_2dataset_bydiagnosis", cluster2 = T,
                                  totalDeregcir_plot = 339)})  

    
  b <- lapply(drug_pathvals, function(x){plot_hist_drugEffvsC_top40(pathvals_controls = path_vals_phy_C, pathvals_disease = path_vals_phy_D, pathvals_drugEff = x, folder = "on_2dataset_bydiagnosis", cluster2 = F,
                                                                         totalDeregcir_plot = 100)})    
  
  
  ### 02.08.2021 Extract the circuits that change in the correct direction  after the Drug simulation ####

  
  plot_boxplots_drugEffvsCvsD_topDeregulatedCir_withoppositeDirection <- function(pathvals_controls, pathvals_disease, pathvals_drugEff, folder , folder2 , pathways = pathways,
                                                                                  metadata_table = metadata ,top_LimmaALL_table = top_LimmaALL, cluster2 = T){
    
    require(ggplot2)
    require(hipathia)
    
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
    diff$mean_C <- apply(pathvals_controls[diff$cir,], 1, mean)
    diff$mean_D <- apply(pathvals_disease[diff$cir,], 1, mean)
    diff$mean_T <- apply(pathvals_drugEff[diff$cir,], 1, mean)
    
    
    for (i in 1: dim(diff)[1]){
      
      if(log2(diff$mean_D[i]) - log2(diff$mean_C[i])>=0){
        
        diff$FC.DvsC <- log2(diff$mean_D[i]) - log2(diff$mean_C[i])
        diff$FC.TvsD <- log2(diff$mean_T[i]) - log2(diff$mean_D[i])
        diff$disease_regulation[i] <- "UP-REGULATION"
        
        
      }else{
        
        diff$FC.DvsC <- log2(diff$mean_D[i]) - log2(diff$mean_C[i])
        diff$FC.TvsD <- log2(diff$mean_T[i]) - log2(diff$mean_D[i])
        diff$disease_regulation[i] <- "DOWN-REGULATION"
        
      }
      
      
    }
    
    
    for (i in 1: dim(diff)[1]){
      
      if(log2(diff$mean_D[i]) - log2(diff$mean_T[i])>=0){
        
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
    dat_plot$circuits <- factor(dat_plot$circuits, levels =  unique(dat_plot$circuits), ordered = T) ## Order by activity values
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
    
    #### 05.08.2021 Do the FC table for each patient before and after drug treatments, 
    #### and DvsC FC and plot the FC_disease - FC_treatment in a manhattan plot style STORED IN FOLDER2 !!!! ####
    
    ## Create the FC table with the 339 top dereg circuits
    
    FC_df <- data.frame(circuit = rownames(top_LimmaALL_table),
                        cir_name = get_path_names(pathways,rownames(top_LimmaALL_table)),
                        FC_DvsC = top_LimmaALL_table$logFC,
                        FC_TvsD =  apply(pathvals_drugEff[rownames(top_LimmaALL_table),], 1,
                                         function(x) log2(mean(x))) - apply(pathvals_disease[rownames(top_LimmaALL_table),], 1,
                                                                              function(x) log2(mean(x))))
    
    FC_df <- FC_df[-which(FC_df$FC_TvsD == 0),]
    
    FC_df$sumFCs <- sapply(1:dim(FC_df)[1], function(x) { FC_df$FC_DvsC[x] + FC_df$FC_TvsD[x] })
    
    write.table(FC_df, file = here("results", "KO_drugs", folder2, paste0(drug_name,"FC_circuit_Oriented.tsv")), 
                row.names = F, col.names = T, quote = F, sep = "\t") 
    
    
    ## Create the FC table with the 339 top deregulated circuits per patient. 
    ## For this purpose we will ONLY USE THE SINGULAR MEAN VALUE OF THE log2FC of the pathvals from controls,  for the FC_DvsC.
    
    ### UPDATED ON 10/09/2021 with MARIAs FEEDBACK: First calculate the log2FC on each comparison
    
    ## Create the simulated control patient sample with the top-dereg circuits from top_limmaALL_table
    
    control_sample <- sapply(pathvals_controls[rownames(top_LimmaALL_table),],1, FUN = mean) ## Median value of each circuit across control samples
    
    FC_patients<- data.frame(patient = colnames(pathvals_disease),
                               FC_DvsC = apply(pathvals_disease[rownames(top_LimmaALL_table),], 2, function(x) mean(x - control_sample)) ,
                               FC_TvsD = apply(pathvals_drugEff[rownames(top_LimmaALL_table),]-pathvals_disease[rownames(top_LimmaALL_table),],2, mean), 
                               FC_TvsC = apply(pathvals_drugEff[rownames(top_LimmaALL_table),], 2, function(x) mean(x - control_sample)),
                               stringsAsFactors = T)
    
    # FC_patients <- data.frame(patient = colnames(pathvals_disease),
    #                           FC_DvsC = apply(pathvals_disease[rownames(top_LimmaALL_table),], 2,
    #                                           function(x) log2(x)) - apply(pathvals_controls[rownames(top_LimmaALL_table),],2,
    #                                                                              function(x) log2(mean(x))) ,
    #                           FC_TvsD =  apply(pathvals_drugEff[rownames(top_LimmaALL_table),], 2,
    #                                            function(x) log2(mean(x))) - apply(pathvals_disease[rownames(top_LimmaALL_table),], 2,
    #                                                                               function(x) log2(mean(x))))
     
    
    FC_patients$sumFC_DvsCplusTvsD <- sapply(1:dim(FC_patients)[1], function(x) { FC_patients$FC_DvsC[x] + FC_patients$FC_TvsD[x] })
    
    
    write.table(FC_patients, file = here("results", "KO_drugs", folder2, paste0(drug_name,"FC_patient_Oriented.tsv")), 
                row.names = F, col.names = T, quote = F, sep = "\t") 
    
    
    FC_plot_patient <- data.frame(patients = FC_patients$patient,
                        FC_difference = FC_patients$sumFCs,
                        ranked_patients = as.character(rank(FC_patients$sumFCs)),
                        Diagnosis = metadata_table$diagnosis[match(FC_patients$patient, metadata_table$samples)] ,
                        Cluster =  as.character(metadata_table$cluster[match(FC_patients$patient, metadata_table$samples)]) ,
                        stringsAsFactors = F) %>% .[order(.$ranked_patients, decreasing = F),]
    
    message("Plotting the LogFC_diff ")  
    ## PLOT IT
    FCs <- ggplot(data = FC_plot_patient,
                  mapping = aes(x =  Diagnosis,
                                y = logFC_difference,
                                color =  Cluster)) +
      ggtitle(paste0('Log2 FC resulting from adding the logFC for each patient before (D = disease) and after simulation of ',drug_name ,' drug effect (T = treated) to log2 FC D vs C' )) +
      geom_point(size = 2) +
      theme_bw()+
      theme(axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.title =  element_text(size = 24),
            plot.title = element_text(size = 20, face = "bold"))
    
    ggsave(filename = here("results", "KO_drugs", folder2, paste0(drug_name, title_plot,"FCplot_patients.png")), plot = FCs, dpi = 400, width = 18, height = 8)
    
    return(dat_plot)
    
  }
  
  
  
  test <- lapply(drug_pathvals, function(x){plot_boxplots_drugEffvsCvsD_topDeregulatedCir_withoppositeDirection(pathvals_controls = path_vals_phy_C, pathvals_disease = path_vals_phy_D, pathvals_drugEff = x, pathways = pathways,
                                                                                                                folder = "boxplots_TtestRelevant_correctDirectionCir",folder2 = "FC_manhattanPlot_tables",  cluster2 = T)})     
  saveRDS(test, file = here("rds", "test.rds"))
  
  test2 <- lapply(drug_pathvals, function(x){plot_boxplots_drugEffvsCvsD_topDeregulatedCir_withoppositeDirection(pathvals_controls = path_vals_phy_C, pathvals_disease = path_vals_phy_D, pathvals_drugEff = x, pathways = pathways,
                                                                                                                 folder = "boxplots_TtestRelevant_correctDirectionCir",folder2 = "FC_manhattanPlot_tables",  cluster2 = F)})     
  
  saveRDS(test2, file = here("rds", "test2.rds"))
  