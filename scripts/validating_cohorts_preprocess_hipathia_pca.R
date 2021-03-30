      ######################################################################
############ Load normalization and Hipathia preprocess for new data ###########
   ##########   after Guillermo's mail, test and validation cohorts ######
      ##################################################################

      library(pacman)
      pacman::p_load("here", "dplyr","reutils", "hipathia", "utils", "stringr", "edgeR","data.table",
                     "AnnotationDbi", "org.Hs.eg.db", "tidyr", "preprocessCore", "ggplot2", "xlsx")
      
      
      expreset_raw <- fread(file = "/home/m3m/INFO_PROYECTO/autoinmunes/data/feedback_guillermo_Marzo2021/CS.Counts.csv",header = T) %>% as.data.frame(.)
      
      metadata <- fread(file = "/home/m3m/INFO_PROYECTO/autoinmunes/data/feedback_guillermo_Marzo2021/CS.Metadata.csv")
      metadata$Center <- gsub(" ", "_", metadata$Center)
      
      rownames(expreset_raw) <- expreset_raw$ID
      expreset_raw <- expreset_raw[,-1]
      print("read...done")
      
      ### 1. Explore and Normalize RNA-seq data wiht edgeR ####
      # Explore how our data is organized
      hist(as.numeric(expreset_raw[2,]),breaks=3)
      var(as.numeric(expreset_raw[2,]))
      
      getVari <- apply(expreset_raw, 1, var)
      hist(getVari,100)
      
      #  NOOOO, la varianza dew nuestros datos es bajita. Clean all low variance rows . Var < 1 -> Como la varianza media de nuestros datos es muy alta aquella que sean menores que uno no nos va a aportar mucha infromación.
      # Sólo quitaremos las varianzas pequeñas y las observaciones que tengas todos sus valores = 0 a cuando estemos mirando la expresion diferencial, si no, nos estaríamos quitando datos importantes.
      
      # idx <- which(getVari<1)
      # length(idx)
      # # expreset <- expreset_raw[-idx, ] -> sólo cuando miremos expresión diferencial
      # print("idex to clean...done")
      
      print("clean...done")
      
      # Normalization by TMM with "edgeR" package
      dge <- DGEList(counts=expreset_raw)
      print("dge...done")
      tmm <-  calcNormFactors(dge, method="TMM")
      print("tmm...done")
      logcpm <- cpm(tmm, prior.count= 0.5, log=TRUE)
      print("dge...done")
      # eliminate from rownames the ".number", beacuse Hipathia do not process them well
      # rownames(logcpm)<-gsub("\\..*", "", rownames(logcpm))
      # print("normalization 1...done")
      
      saveRDS(logcpm, file = here("rds", "logcpm_validating_cohorts.rds"))
      
      ### 2. Run HiPathia and LIMMA for all samples ####
      
      ## PREPROCESSMENT: Data scaling and normalization
      
      trans_data <- translate_data(logcpm, "hsa")
      exp_data <- normalize_data(trans_data)
      
      ## Loading Pathways (all and only physiological)
      
      path_list <- read.table(file = here("data", "physiological_paths.tsv"), sep = "\t") #physiological_pathways
      pathways <- load_pathways("hsa")
      pathways_phy <- load_pathways("hsa", pathways_list = path_list$V2 )
      
      subpathways.list <- lapply(pathways_phy$pathigraphs, function(x){names(x$effector.subgraphs)})
      
      subpathways_phy <- data.frame(unlist(subpathways.list), stringsAsFactors = F)
      
      colnames(subpathways_phy) <- "hipathia"
      
      ## Using Hipathia to compute the signal (subpathways) 
      
      results <- hipathia(exp_data, pathways, decompose = FALSE, verbose=FALSE)
      saveRDS(results, file = here("rds","validation_cohort_Hiresults.rds"))
      
      path_vals <- get_paths_data(results, matrix = TRUE)
      
      path_vals_phy <- path_vals[which(rownames(path_vals) %in% subpathways_phy$hipathia),]
      
      ## Load updated annotations and annotation object to subtitute
      hp <- hipathia:::hub()
      annot_uniprot_hsa <- hp[[names(hp)[hp$title == "annot_uniprot_hsa.rda"]]]
      
      update_anot <- read.delim(file = "/home/m3m/INFO_PROYECTO/Hipathia/results/hipathia_FULLannotations_26032021.tsv", header = T)
      update_anot <- update_anot[,c(2,3)]
      colnames(update_anot) <- c("gene", "function")
      
      hp[[names(hp)[hp$title == "annot_uniprot_hsa.rda"]]] <- update_anot #### NO SIRVE HABLAR CON KINZA PARA ARREGLARLO!!!!
      
      
      ## LIMMA FOR ALL SAMPLES 
      
      target_matrix <- data.frame(samples = colnames(path_vals), cluster = metadata$Cluster[match(colnames(path_vals), metadata$OMICID)],
                                  diagnosis = metadata$Diagnosis[match(colnames(path_vals), metadata$OMICID)],
                                  age = metadata$Age[match(colnames(path_vals), metadata$OMICID)],
                                  gender = metadata$Gender[match(colnames(path_vals), metadata$OMICID)],
                                  race = metadata$Race[match(colnames(path_vals), metadata$OMICID)],
                                  center = metadata$Center[match(colnames(path_vals), metadata$OMICID)],
                                  cohort = metadata$Cohort[match(colnames(path_vals), metadata$OMICID)],
                                  pool = metadata$POOL[match(colnames(path_vals), metadata$OMICID)],
                                  RIN = metadata$RIN[match(colnames(path_vals), metadata$OMICID)],
                                  stringsAsFactors = F)
      
      target_matrix$type <- target_matrix$diagnosis %>% gsub(pattern = "SSc|SjS|SLE|PAPs|MCTD|UCTD|RA", replacement ="DISEASE", x = .)
      
      write.table(target_matrix, file = here("data", "feedback_guillermo_Marzo2021", "metadata_reformatted.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
      
      sample_group_limma <- target_matrix[, c(1,4,5,7,9,10,11)] # Note that for Limma the Control has to be first in the levels
      
      
      design <- model.matrix( ~ 0  + type + age + gender + center + pool + RIN ,  data = sample_group_limma)
      colnames(design) <- gsub("typeCTRL", "C", colnames(design)) 
      colnames(design) <- gsub("typeDISEASE", "disease", colnames(design)) 
      rownames(design) <- sample_group_limma$samples
      cont.matrix <- makeContrasts( diseasevsC=disease-C, levels=design)
      
      fit <- lmFit(path_vals_phy , design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit3 <- eBayes(fit2)
      
      tableLimma <- topTable(fit3, number = rownames(path_vals), adjust.method="bonferroni", sort.by="p")
      top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
      top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
      
      anot_all <-data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
                            logFC = top_LimmaALL$logFC,
                            P.val = top_LimmaALL$P.Value,
                            FDR.Pval = top_LimmaALL$adj.P.Val,
                            UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
                            GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
                            stringsAsFactors = F)
      
      write.xlsx(anot_all, file = here( "results", "validation_allvsC_afterlimma_annot.xlsx"))
      
      pathways_deregAll <- str_split(anot_all$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.) 
      
      length(pathways_deregAll) # 61
      
      write.xlsx(pathways_deregAll, file = here("results","validation_pathways_allvsC.xlsx"),row.names = F)
      
      ### 3. PCA over pathvals ####
      
      ## Reorder de data to have the same sample order as in the colnames of our data
      
      metadata_pca <- metadata[order(match(metadata$OMICID,colnames(path_vals))), ] %>% mutate_all(as.character) ## reorder
      
      ## Quality of the data : PCA over pathvals
      
      library(ggplot2)
      
      start_time <- Sys.time()
      
      pathvals_pca <- prcomp(t(path_vals[,which(!colnames(path_vals) %in% metadata_pca$OMICID[metadata_pca$Cluster %in% c("0")])]),scale. = F)
      
      percentVar <- round(100*pathvals_pca$sdev^2/sum(pathvals_pca$sdev^2),1)
      sd_ratio <- sqrt(percentVar[2] / percentVar[1])
      
      dataGG <- data.frame(PC1 = pathvals_pca$x[,1], PC2 = pathvals_pca$x[,2],
                           Cluster = metadata_pca$Cluster[!metadata_pca$Cluster%in% c("0")])
      # Center = metadata_pca$PublicCenter)
      dev.new()
      ggplot(dataGG, aes(PC1, PC2)) +
         geom_point(aes(colour = Cluster), size = 3) +
         ggtitle("PCA of pathways activity values") +
         xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
         ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
         theme_bw()+
         theme(plot.title = element_text(hjust = 0, face = "bold", size = 34),
               legend.title = element_text (size = 30, face = "bold"),
               legend.text = element_text(size = 22),
               axis.title.y = element_text(hjust = 0.5 ,size = 22),
               axis.title.x = element_text(hjust = 0.5 ,size = 22))
      # axix.text.y = element_blank(),
      # axix.text.x = element_blank()) 
      # coord_fixed(ratio = sd_ratio) +
      # scale_shape_manual(values = c(4,15))+
      # scale_color_manual(values = c("darkorange2", "dodgerblue4"))+
      
       ggsave(file = here("results","PCA_pathvals_validationAll_withoutC0.png"), plot = last_plot(), dpi = 300)
      
      # ggsave(file = here("results","PCA_pathvals_validationAll_withoutC2yC0.png"), plot = last_plot(), dpi = 300)
      
      end_time <- Sys.time()
      
      end_time - start_time
      
      ### 4. PCA over expreset ####
      
      ## Quality of the data 
      # 
      # library(ggplot2)
      # 
      # start_time <- Sys.time()
      # 
      # logcpm_pca <- prcomp(t(logcpm),scale. = F)
      # 
      # percentVar <- round(100*logcpm_pca$sdev^2/sum(logcpm_pca$sdev^2),1)
      # sd_ratio <- sqrt(percentVar[2] / percentVar[1])
      # 
      # dataGG <- data.frame(PC1 = logcpm_pca$x[,1], PC2 = logcpm_pca$x[,2],
      #                      Cluster = metadata_pca$DISEASE)
      # dev.new()
      # ggplot(dataGG, aes(PC1, PC2)) +
      #    geom_point(aes(colour = Cluster), size = 3) +
      #    ggtitle("PCA of expressionSet over clusters") +
      #    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
      #    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
      #    theme_bw()+
      #    theme(plot.title = element_text(hjust = 0, face = "bold", size = 34),
      #          legend.title = element_text (size = 30, face = "bold"),
      #          legend.text = element_text(size = 22),
      #          axis.title.y = element_text(hjust = 0.5 ,size = 22),
      #          axis.title.x = element_text(hjust = 0.5 ,size = 22))
      # # axix.text.y = element_blank(),
      # # axix.text.x = element_blank()) 
      # #coord_fixed(ratio = sd_ratio) +
      # # scale_shape_manual(values = c(4,15))+
      # # scale_color_manual(values = c("darkorange2", "dodgerblue4"))+
      # 
      # 
      # # ggsave(file = here("results","PCA_expreset_bloodTransc.png"), plot = last_plot(), dpi = 300)
      # 
      # end_time <- Sys.time()
      # 
      # end_time - start_time
      
      