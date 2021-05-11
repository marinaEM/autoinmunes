      ##########################################################################################
#### Calculate the effect of the Drug by using already calculated drug simulations on disease  ######
  ####  with Hipathia in patients data and now calculate eucl.distance with Controls      #####
      ######################################################################################


library(pacman)
  
devtools::install_github("thomasp85/patchwork")

pacman::p_load("here", "hipathia", "utils", "stringr", "genefilter", "tidyr","data.table","limma", "ggplot2", "patchwork")


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

  a <- lapply(drug_pathvals, function(x){plot_hist_drugEffvsC_top40(pathvals_controls = path_vals_phy_C, pathvals_disease = path_vals_phy_D, pathvals_drugEff = x, folder = "on_2dataset_bydiagnosis", cluster2 = T,
                                  totalDeregcir_plot = 339)})  

    
  b <- lapply(drug_pathvals, function(x){plot_hist_drugEffvsC_top40(pathvals_controls = path_vals_phy_C, pathvals_disease = path_vals_phy_D, pathvals_drugEff = x, folder = "on_2dataset_bydiagnosis", cluster2 = F,
                                                                         totalDeregcir_plot = 100)})    
    