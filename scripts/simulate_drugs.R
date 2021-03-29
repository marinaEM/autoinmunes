#############################################################################
#### Simulate the effect of Anakinra drug on each sample and compare ######
#############################################################################

library(pacman)
devtools::install_github("thomasp85/patchwork")

pacman::p_load("here", "hipathia", "utils", "stringr", "genefilter", "tidyr","data.table", "xlsx","limma", "ggplot2")

### 1. Load data and perform gene action of ANAKINRA --> IL1R1 entrezID:3554 ####

## Load normalized expression data
logcpm <- readRDS(file = here("rds", "logcpm_bloodTranscriptome.rds"))
metadata_init <- fread(file = "/home/m3m/INFO_PROYECTO/autoinmunes/data/Blood.Transcriptome/CS.Transcriptome.Metadata.csv")

## Subset only disease samples
logcpm_dis <- logcpm[, which(colnames(logcpm) %in% metadata_init$PublicID[metadata_init$DISEASE != "CONTROLS"])]

metadata_init <- metadata_init[metadata_init$DISEASE != "CONTROLS",]

metadata_init$type <- "controls"

## Loading Pathways (all and only physiological)

path_list <- read.table(file = here("data", "physiological_paths.tsv"), sep = "\t") #physiological_pathways
pathways <- load_pathways("hsa")
pathways_phy <- load_pathways("hsa", pathways_list = path_list$V2 )

subpathways.list <- lapply(pathways_phy$pathigraphs, function(x){names(x$effector.subgraphs)})

subpathways_phy <- data.frame(unlist(subpathways.list), stringsAsFactors = F)

colnames(subpathways_phy) <- "hipathia"

"3554" %in% pathways$all.genes ## Checking that ANAKIRNA's target is in Hipathia

metadata_KO <- metadata_init
metadata_KO$PublicID <- paste0(metadata_KO$PublicID, "_3354KO")
metadata_KO$type <- "treated"

#### 2.  KO of gene "3554" ANAKINRA --> IL1R1 entrezID:3554 ######

## Run HiPathia and LIMMA for all samples  besides controls ###

trans_data <- translate_data(logcpm_dis, "hsa")
"3554" %in% rownames(trans_data) ## Checking that ANAKIRNA's target is in the dataset

trans_data_KO <- trans_data
trans_data_KO["3554", ] <- 0.0001
colnames(trans_data_KO) <- paste0(colnames(trans_data), "_3354KO")

trans_data_KO_C <- cbind(trans_data_KO, trans_data)

exp_data_KO_C <- normalize_data(trans_data_KO_C)


    ## Using Hipathia to compute the signal (subpathways) 
    
    results_KO <- hipathia(exp_data_KO_C, pathways, decompose = FALSE, verbose=FALSE)
    saveRDS(results_KO, file = here("rds","KO_anakinra_C_bloodTranscriptome_Hiresults.rds"))
    
    path_vals_KO <- get_paths_data(results_KO, matrix = TRUE)
    
    path_vals_phy_KO <- path_vals_KO[which(rownames(path_vals_KO) %in% subpathways_phy$hipathia),]

    hist(path_vals_phy_KO)
    
    anot_cir <- data.frame(Circuit = get_path_names(pathways, rownames(path_vals_phy_KO)),
                           UniprotKB = get_pathways_annotations(rownames(path_vals_phy_KO), pathways, dbannot= "uniprot" ,collapse = T),
                           GO = get_pathways_annotations(rownames(path_vals_phy_KO), pathways, dbannot= "GO" ,collapse = T), 
                           stringsAsFactors = F)
    
    
    ## LIMMA FOR ALL SAMPLES 
    
    metadata <- rbind(metadata_KO, metadata_init)
    
    target_matrix_KO <- data.frame(samples = colnames(path_vals_phy_KO), 
                                cluster = metadata$DISEASE[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                diagnosis = metadata$Diagnosis[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                age = metadata$Age[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                gender = metadata$Gender[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                antimalarials = metadata$Antimalarials[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                immunosuppresants = metadata$Immunosuppresants[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                biologicals = metadata$Biologicals[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                steroids = metadata$Steroids[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                systemic.antibiotics = metadata$Systemic.Antibiotics[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                public.center =metadata$PublicCenter[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                pool = metadata$Pool[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                type = metadata$type[match(colnames(path_vals_phy_KO), metadata$PublicID)],
                                stringsAsFactors = F)

    
    sample_group_limma <- target_matrix_KO[, c(1,3,4,5,6,7,8,9,10,13)] # Note that for Limma the Control has to be first in the levels
    
    
    design <- model.matrix( ~ 0  + type + age + gender + antimalarials + immunosuppresants + biologicals + steroids + systemic.antibiotics ,  data = sample_group_limma)
    colnames(design) <- c( "C", "treated", "age", "gender", "antimalarials", "immunosuppresants","biologicals", "steroids","systemic.antibiotics")
    rownames(design) <- sample_group_limma$samples
    cont.matrix <- makeContrasts( treatedvsC=treated-C, levels=design)
    
    fit <- lmFit(path_vals_phy_KO , design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)
    
    tableLimma <- topTable(fit3, number = rownames(path_vals_phy_KO), adjust.method="fdr", sort.by="p")
    top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    
    anot_all <-data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
                          logFC = top_LimmaALL$logFC,
                          P.val = top_LimmaALL$P.Value,
                          FDR.Pval = top_LimmaALL$adj.P.Val,
                          UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
                          GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
                          stringsAsFactors = F)
    
        # write.xlsx(anot_all, file = here( "results", "anakinra_KO_bloodtranscriptome_allvsC_afterlimma_annot.xlsx")) 
    
    #### 3. Make comparison plots of the significant deregulate circuits ###
    
    dereg_circuits_KO <- rownames(top_LimmaALL)
    
    dereg_cir_pathvals <- t(path_vals_phy_KO[dereg_circuits_KO,]) %>% as.data.frame(.)
    
    dereg_cir_pathvals$type <- c(rep("treated", 1181), rep("disease", 1181))

    # boxplot(P-hsa04064-115 ~ type, data = dereg_cir_pathvals)    
    
    
    #### 4. Calculate  and plot the logFC for each sample after KO ###
    
    library(patchwork)
    
    anakinra <- apply(path_vals_phy_KO, 2, mean)
    
    FC_anakinra <- log2(anakinra[1:918] /anakinra[919:length(anakinra)])
    
    names(FC_anakinra) <- gsub("_3354KO", "", names(FC_anakinra))
    
    FC_anakinra <- data.frame(patients = names(FC_anakinra),
                              logFC = abs(FC_anakinra),
                              ranked_patients = rank(FC_anakinra),
                              Diagnosis = metadata_init$Diagnosis[match(names(FC_anakinra), metadata_init$PublicID)] ,
                              Cluster =  metadata_init$DISEASE[match(names(FC_anakinra), metadata_init$PublicID)] ,
                              stringsAsFactors = F) %>% .[order(.$ranked_patients, decreasing = F),]
    
    

    FCs <- ggplot(data = FC_anakinra,
           mapping = aes(x = ranked_patients,
                         y = logFC,
                         color = Diagnosis)) +
            ggtitle('Log FC of the circuits activity absolute values on the ranked patients after simulation of Anakinra drug effect') +
            geom_point(size = 3) +
            geom_hline(yintercept = 0.0003, linetype = "dashed") +
            theme_minimal()+
            theme(axis.title.x = element_text(size = 18),
                  axis.title.y = element_text(size = 18),
                  legend.text = element_text(size = 18),
                  legend.title =  element_text(size = 24),
                  plot.title = element_text(size = 18, face = "bold"))
    
    bars_df_res <- data.frame(table(FC_anakinra$Diagnosis[FC_anakinra$logFC >= 0.0003]))
    
    bars_df_NOres <- data.frame(table(FC_anakinra$Diagnosis[FC_anakinra$logFC < 0.0003]))
        
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
     
    dev.new()
    (FCs |bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1))
    
    ggsave(last_plot(), filename = here("results", "FC_anakinra_diagnosis.png"), dpi = 400 )

    
    ### Do same graph for the clusters insted of the diagnosis
    
    
    FCs <- ggplot(data = FC_anakinra,
                  mapping = aes(x = ranked_patients,
                                y = logFC,
                                color = Cluster)) +
        ggtitle('Log FC of the circuits activity absolute values on the ranked patients after simulation of Anakinra drug effect') +
        geom_point(size = 3) +
        geom_hline(yintercept = 0.0003, linetype = "dashed") +
        theme_minimal()+
        theme(axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text = element_text(size = 18),
              legend.title =  element_text(size = 24),
              plot.title = element_text(size = 18, face = "bold"))
    
    bars_df_res <- data.frame(table(FC_anakinra$Cluster[FC_anakinra$logFC >= 0.0003]))
    
    bars_df_NOres <- data.frame(table(FC_anakinra$Cluster[FC_anakinra$logFC < 0.0003]))
    
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
    
    dev.new()
    (FCs | bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1)) #+ plot_layout(guides = 'collect')  
    
    ggsave(last_plot(), filename = here("results", "FC_anakinra_clusters.png"), dpi = 400 )
    
    

#### 3. SIMULATION OF BELIMUMAB DRUG ( KO OF TARGET TNFSF13B ENTREZID = "10673") ######
    
    ## Run HiPathia and LIMMA for all samples  besides controls ###
    
    metadata_KO2 <- metadata_init
    metadata_KO2$PublicID <- paste0(metadata_KO2$PublicID, "_10673KO")
    metadata_KO2$type <- "treated"
    
    
    trans_data <- translate_data(logcpm_dis, "hsa")
    "10673" %in% rownames(trans_data) ## Checking that BELIMUMAB's target is in the dataset
    
    trans_data_KO2 <- trans_data
    trans_data_KO2["10673", ] <- 0.0001
    colnames(trans_data_KO2) <- paste0(colnames(trans_data), "_10673KO")
    
    trans_data_KO2_C <- cbind(trans_data_KO2, trans_data)
    
    exp_data_KO2_C <- normalize_data(trans_data_KO2_C)
    
    
    ## Using Hipathia to compute the signal (subpathways) 
    
    results_KO2 <- hipathia(exp_data_KO2_C, pathways, decompose = FALSE, verbose=FALSE)
    saveRDS(results_KO2, file = here("rds","KO_belumimab_C_bloodTranscriptome_Hiresults.rds"))
    
    path_vals_KO2 <- get_paths_data(results_KO2, matrix = TRUE)
    
    path_vals_phy_KO2 <- path_vals_KO2[which(rownames(path_vals_KO2) %in% subpathways_phy$hipathia),]
    
    hist(path_vals_phy_KO2)
    
     ## LIMMA FOR ALL SAMPLES 
    
    metadata <- rbind(metadata_KO2, metadata_init)
    
    target_matrix_KO2 <- data.frame(samples = colnames(path_vals_phy_KO2), 
                                   cluster = metadata$DISEASE[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   diagnosis = metadata$Diagnosis[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   age = metadata$Age[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   gender = metadata$Gender[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   antimalarials = metadata$Antimalarials[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   immunosuppresants = metadata$Immunosuppresants[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   biologicals = metadata$Biologicals[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   steroids = metadata$Steroids[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   systemic.antibiotics = metadata$Systemic.Antibiotics[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   public.center =metadata$PublicCenter[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   pool = metadata$Pool[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   type = metadata$type[match(colnames(path_vals_phy_KO2), metadata$PublicID)],
                                   stringsAsFactors = F)
    
    
    sample_group_limma <- target_matrix_KO2[, c(1,3,4,5,6,7,8,9,10,13)] # Note that for Limma the Control has to be first in the levels
    
    
    design <- model.matrix( ~ 0  + type + age + gender + antimalarials + immunosuppresants + biologicals + steroids + systemic.antibiotics ,  data = sample_group_limma)
    colnames(design) <- c( "C", "treated", "age", "gender", "antimalarials", "immunosuppresants","biologicals", "steroids","systemic.antibiotics")
    rownames(design) <- sample_group_limma$samples
    cont.matrix <- makeContrasts( treatedvsC=treated-C, levels=design)
    
    fit <- lmFit(path_vals_phy_KO2 , design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)
    
    tableLimma <- topTable(fit3, number = rownames(path_vals_phy_KO2), adjust.method="fdr", sort.by="p")
    top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    
    anot_all <-data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
                          logFC = top_LimmaALL$logFC,
                          P.val = top_LimmaALL$P.Value,
                          FDR.Pval = top_LimmaALL$adj.P.Val,
                          UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
                          GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
                          stringsAsFactors = F)
    
     write.xlsx(anot_all, file = here( "results", "belimumab_KO_bloodtranscriptome_allvsC_afterlimma_annot.xlsx")) 
    
    #### 3. Make comparison plots of the significant deregulated circuits ###
    
    # dereg_circuits_KO <- rownames(top_LimmaALL)
    # 
    # dereg_cir_pathvals <- t(path_vals_phy_KO[dereg_circuits_KO,]) %>% as.data.frame(.)
    # 
    # dereg_cir_pathvals$type <- c(rep("treated", 1181), rep("disease", 1181))
    # 
    # # boxplot(P-hsa04064-115 ~ type, data = dereg_cir_pathvals)    

    
    #### 4. Calculate  and plot the logFC for each sample after KO ###
    
    library(patchwork)
    
    belimumab <- apply(path_vals_phy_KO2, 2, mean)
    
    FC_belimumab <- log2(belimumab[1:918] /belimumab[919:length(anakinra)])
    
    names(FC_belimumab) <- gsub("_10673KO", "", names(FC_belimumab))
    
    FC_belimumab <- data.frame(patients = names(FC_belimumab),
                              logFC = abs(FC_belimumab),
                              ranked_patients = rank(FC_belimumab),
                              Diagnosis = metadata_init$Diagnosis[match(names(FC_belimumab), metadata_init$PublicID)] ,
                              Cluster =  metadata_init$DISEASE[match(names(FC_belimumab), metadata_init$PublicID)] ,
                              stringsAsFactors = F) %>% .[order(.$ranked_patients, decreasing = F),]
    
    
    
    FCs <- ggplot(data = FC_belimumab,
                  mapping = aes(x = ranked_patients,
                                y = logFC,
                                color = Diagnosis)) +
        ggtitle('Log FC of the circuits activity absolute values on the ranked patients after simulation of Belimumab drug effect') +
        geom_point(size = 3) +
        geom_hline(yintercept = 0.003, linetype = "dashed") +
        theme_minimal()+
        theme(axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text = element_text(size = 18),
              legend.title =  element_text(size = 24),
              plot.title = element_text(size = 18, face = "bold"))
    
    bars_df_res <- data.frame(table(FC_belimumab$Diagnosis[FC_belimumab$logFC >= 0.003]))
    
    bars_df_NOres <- data.frame(table(FC_belimumab$Diagnosis[FC_belimumab$logFC < 0.003]))
    
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
    
    dev.new()
    (FCs |bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1))
    
    ggsave(last_plot(), filename = here("results", "FC_belimumab_diagnosis.png"), dpi = 400 )
    
    
    ### Do same graph for the clusters insted of the diagnosis
    
    
    FCs <- ggplot(data = FC_belimumab,
                  mapping = aes(x = ranked_patients,
                                y = logFC,
                                color = Cluster)) +
        ggtitle('Log FC of the circuits activity absolute values on the ranked patients after simulation of Belimumab drug effect') +
        geom_point(size = 3) +
        geom_hline(yintercept = 0.003, linetype = "dashed") +
        theme_minimal()+
        theme(axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text = element_text(size = 18),
              legend.title =  element_text(size = 24),
              plot.title = element_text(size = 18, face = "bold"))
    
    bars_df_res <- data.frame(table(FC_belimumab$Cluster[FC_belimumab$logFC >= 0.003]))
    
    bars_df_NOres <- data.frame(table(FC_belimumab$Cluster[FC_belimumab$logFC < 0.003]))
    
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
    
    dev.new()
    (FCs |bars_responders/plot_spacer()+ bars_nonresponders) + plot_layout(widths = c(3, 1))
    
    ggsave(last_plot(), filename = here("results", "FC_belimumab_cluster.png"), dpi = 400 )
        
   
    
    
    
    