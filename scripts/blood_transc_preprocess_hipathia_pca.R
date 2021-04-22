#############################################################################
#### Load normalization and Hipathia preprocess for Marta's group data ######
#############################################################################

library(pacman)
pacman::p_load("here", "rentrez","reutils", "hipathia", "biomaRt", "utils", "stringr", "genefilter", "edgeR",
               "AnnotationDbi", "org.Hs.eg.db", "hgfocus.db", "tidyr", "preprocessCore", "data.table", "ggplot2")

library("openxlsx")

expreset_raw <- fread(file = "/home/m3m/INFO_PROYECTO/autoinmunes/data/Blood.Transcriptome_Julio2020/CS.Transcriptome.Counts.csv") %>% as.data.frame(.)

metadata <- fread(file = "/home/m3m/INFO_PROYECTO/autoinmunes/data/Blood.Transcriptome_Julio2020/CS.Transcriptome.Metadata.csv")

rownames(expreset_raw) <- expreset_raw$V1
expreset_raw <- expreset_raw[,-1]
print("read...done")

# Clean from 0 value expression genes and low variance

# Because low variance implies all 0  we don't need to clean from all 0 expression genes but it would be like:
# cl <-which(apply(expreset_raw, 1, function(x) all(x==0)))
# length(cl)
# expreset <- expreset_raw[-cl, ]
# param <- 1e-4
# expreset <- expreset[getVari > param & !is.na(getVar), ] 

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
  
  # saveRDS(logcpm, file = here("rds", "logcpm_bloodTranscriptome.rds"))

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
    # saveRDS(results, file = here("rds","marta_bloodTranscriptome_Hiresults.rds"))
    
    path_vals <- get_paths_data(results, matrix = TRUE)
    
    path_vals_phy <- path_vals[which(rownames(path_vals) %in% subpathways_phy$hipathia),]

    
    
    ## LIMMA FOR ALL SAMPLES 
    
    target_matrix <- data.frame(samples = colnames(path_vals), cluster = metadata$DISEASE[match(colnames(path_vals), metadata$PublicID)],
                                diagnosis = metadata$Diagnosis[match(colnames(path_vals), metadata$PublicID)],
                                age = metadata$Age[match(colnames(path_vals), metadata$PublicID)],
                                gender = metadata$Gender[match(colnames(path_vals), metadata$PublicID)],
                                antimalarials = metadata$Antimalarials[match(colnames(path_vals), metadata$PublicID)],
                                immunosuppresants = metadata$Immunosuppresants[match(colnames(path_vals), metadata$PublicID)],
                                biologicals = metadata$Biologicals[match(colnames(path_vals), metadata$PublicID)],
                                steroids = metadata$Steroids[match(colnames(path_vals), metadata$PublicID)],
                                systemic.antibiotics = metadata$Systemic.Antibiotics[match(colnames(path_vals), metadata$PublicID)],
                                public.center =metadata$PublicCenter[match(colnames(path_vals), metadata$PublicID)],
                                pool = metadata$Pool[match(colnames(path_vals), metadata$PublicID)],
                                RIN = metadata$RIN[match(colnames(path_vals), metadata$PublicID)],
                                stringsAsFactors = F)
    
    target_matrix$type <- target_matrix$cluster %>% gsub(pattern = "CLUSTER1|CLUSTER2|CLUSTER3|CLUSTER4", replacement ="DISEASE", x = .)
    
    sample_group_limma <- target_matrix[, c(1,3,4,5,6,7,8,9,10,13)] # Note that for Limma the Control has to be first in the levels

    
    design <- model.matrix( ~ 0  + type + age + gender + antimalarials + immunosuppresants + biologicals + steroids + systemic.antibiotics + RIN + pool ,  data = target_matrix)
    # colnames(design) <- c( "C", "disease", "age", "gender", "antimalarials", "immunosuppresants","biologicals", "steroids","systemic.antibiotics")
    colnames(design) <- gsub("typeCONTROLS", "C", colnames(design)) 
    colnames(design) <- gsub("typedisease", "disease", colnames(design)) 
    rownames(design) <- sample_group_limma$samples
    cont.matrix <- makeContrasts( diseasevsC=disease-C, levels=design)
    
    fit <- lmFit(path_vals_phy , design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)
    
    #tableLimma <- topTable(fit3, number = rownames(path_vals), adjust.method="fdr", sort.by="p") 
    tableLimma <- topTable(fit3, number = rownames(path_vals), adjust.method="bonferroni", sort.by="p")  ## Lo hcemos con bonferroni a petición de Guillermo
    top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    
    anot_all <-data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
                          logFC = top_LimmaALL$logFC,
                          P.val = top_LimmaALL$P.Value,
                          FDR.Pval = top_LimmaALL$adj.P.Val,
                          UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
                          GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
                          stringsAsFactors = F)
   
    # write.xlsx(anot_all, file = here( "results", "bloodtranscriptome_allvsC_afterlimma_annot.xlsx"))
    
     write.xlsx(anot_all, file = here( "results", "bloodtranscriptome_allvsC_afterlimma_annot_bonferroni.xlsx"))
    
    pathways_deregAll <- str_split(anot_all$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.) ## 61
    
    # write.xlsx(pathways_deregAll, file = here("results","pathways_allvsC.xlsx"),row.names = F)
    
    write.xlsx(pathways_deregAll, file = here("results","pathways_allvsC_bonferroni.xlsx"),row.names = F)
    
### 3. PCA over pathvals ####
    
    ## Reorder de data to have the same sample order as in the colnames of our data
    
    metadata_pca <- metadata[order(match(metadata$PublicID,colnames(path_vals))), ] ## reorder
    
    ## Quality of the data : PCA over pathvals
    
    library(ggplot2)
    
    start_time <- Sys.time()
    
    pathvals_pca <- prcomp(t(path_vals[,which(!colnames(path_vals) %in% metadata_pca$PublicID[metadata_pca$DISEASE == "CLUSTER2"])]),scale. = F)
    
    percentVar <- round(100*pathvals_pca$sdev^2/sum(pathvals_pca$sdev^2),1)
    sd_ratio <- sqrt(percentVar[2] / percentVar[1])
    
    dataGG <- data.frame(PC1 = pathvals_pca$x[,1], PC2 = pathvals_pca$x[,2],
                         Cluster = metadata_pca$DISEASE[metadata_pca$DISEASE!= "CLUSTER2"])
                         # Center = metadata_pca$PublicCenter)
    dev.new()
    ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(colour = Cluster), size = 3) +
      ggtitle("PCA of pathvals") +
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
    #coord_fixed(ratio = sd_ratio) +
    # scale_shape_manual(values = c(4,15))+
    # scale_color_manual(values = c("darkorange2", "dodgerblue4"))+
    
    ggsave(file = here("results","PCA_pathvals_bloodTransc_witoutC2.png"), plot = last_plot(), dpi = 300)
    
    # ggsave(file = here("results","PCA_pathvals_bloodTransc.png"), plot = last_plot(), dpi = 300)
    
    end_time <- Sys.time()
    
    end_time - start_time
    
### 4. PCA over expreset ####
    
    ## Quality of the data 
    
    library(ggplot2)
    
    start_time <- Sys.time()
    
    logcpm_pca <- prcomp(t(logcpm),scale. = F)
  
    percentVar <- round(100*logcpm_pca$sdev^2/sum(logcpm_pca$sdev^2),1)
    sd_ratio <- sqrt(percentVar[2] / percentVar[1])
    
    dataGG <- data.frame(PC1 = logcpm_pca$x[,1], PC2 = logcpm_pca$x[,2],
                         Cluster = metadata_pca$DISEASE)
    dev.new()
    ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(colour = Cluster), size = 3) +
      ggtitle("PCA of expressionSet over clusters") +
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
    #coord_fixed(ratio = sd_ratio) +
    # scale_shape_manual(values = c(4,15))+
    # scale_color_manual(values = c("darkorange2", "dodgerblue4"))+
    
    
    # ggsave(file = here("results","PCA_expreset_bloodTransc.png"), plot = last_plot(), dpi = 300)
    
    end_time <- Sys.time()
    
    end_time - start_time
    

### 5. Run Hipathia on each disease vs control ####
    
    ## LIMMA FOR EACH DISEASE
    
    ### First we divide the matrix for each disease and controls and then paste the controls to each disease for the comparisons
     
    target_matrixes_diseases <- list()
    
    for (i in names(table(target_matrix$diagnosis))) {
      
      target_matrixes_diseases[[i]] <- target_matrix[grep(pattern = i, x = target_matrix$diagnosis), ] 
  
    }
    
    controls <- target_matrixes_diseases[["CTRL"]]
    
    target_matrixes_diseases <- target_matrixes_diseases[-1]
    
    target_matrixes_diseases <- lapply(target_matrixes_diseases, function (x) {rbind(x, controls)}) ## Pasting the controls to the diaseases
    

  ## Testing data
    
    # pathvals <- path_vals_phy
    # matrix_metadata <- target_matrixes_diseases[[2]]
    
 ## Making a function to apply to each disease vs control comparison 
    
    #     get_diff_expressed_pathways_diseases <- function(matrix_metadata, pathvals){
    #     
    #     disease <- names(table(matrix_metadata$diagnosis))[names(table(matrix_metadata$diagnosis)) != "CTRL"]
    #     message(paste0("performing limma and extracting significant circuits for... ", disease))
    #   
    #     ## Select the columns corresponding to each disease
    #     samples <- matrix_metadata$samples
    #     idx <- which(colnames(pathvals) %in% samples)
    #     disease_path_vals <- pathvals[, idx]
    #     
    #     ##  Disease vs control limma
    #     sample_group_limma <- matrix_metadata[, c(1,3,4,5,6,7,8,9,10,13)]
    #     design <- model.matrix( ~ 0  + type + age + gender, data = sample_group_limma ) #,+ antimalarials + immunosuppresants + biologicals + steroids + systemic.antibiotics ,  data = sample_group_limma)
    #     colnames(design) <- c( "C", "disease", "age", "gender") #, "antimalarials", "immunosuppresants","biologicals", "steroids","systemic.antibiotics")
    #     rownames(design) <- sample_group_limma$samples
    #     cont.matrix <- makeContrasts( diseasevsC=disease-C, levels=design)
    #     
    #     fit <- lmFit(disease_path_vals , design)
    #     fit2 <- contrasts.fit(fit, cont.matrix)
    #     fit3 <- eBayes(fit2)
    #     tableLimma <- topTable(fit3, number = rownames(path_vals), adjust.method="fdr", sort.by="p")
    #     top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
    #     top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
    #     
    #     anot_all <- data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
    #                            logFC = top_LimmaALL$logFC,
    #                            P.val = top_LimmaALL$P.Value,
    #                            FDR.Pval = top_LimmaALL$adj.P.Val,
    #                             UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
    #                             # GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
    #                            stringsAsFactors = F)
    #     
    #     # write.xlsx(anot_all, file = here( "results", paste0(disease, "vsC_afterlimma_annot.xlsx")),row.names = F)  
    #     
    #     pathways_dereg <- str_split(anot_all$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.)
    #     
    #     message(length(pathways_dereg))
    #     message(length(anot_all$Circuit))
    #     
    #     write.xlsx(pathways_dereg, file = here( "results", paste0(disease, "pathways_dereg.xlsx")),row.names = F, col.names = F) 
    #     
    #   return(anot_all)
    # }
    
 
    diseases_limma_pathvals <- lapply(target_matrixes_diseases, function(x){get_diff_expressed_pathways_diseases(x, path_vals_phy)}) 
    # saveRDS(diseases_limma_pathvals, file = here("rds", "diseases_limma_pathvals.rds"))
    
    pathways_deP_SADs <- lapply(diseases_limma_pathvals, function(x){ str_split(x$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.)})

### 6. Analyze the differences and common deregulated circuits among comparisons ####
    
    ### Get only the circuits that are uniquely deregulated in each disease vs control comparison and common ones
    
     table_circuits <- data.frame(code = rownames(path_vals), name = get_path_names(pathways, rownames(path_vals)))
    
     common_cir_between_dis <- Reduce(intersect, list(diseases_limma_pathvals[[1]]$Circuit, diseases_limma_pathvals[[2]]$Circuit, diseases_limma_pathvals[[3]]$Circuit,
                                    diseases_limma_pathvals[[4]]$Circuit,diseases_limma_pathvals[[5]]$Circuit,diseases_limma_pathvals[[6]]$Circuit,
                                    diseases_limma_pathvals[[7]]$Circuit))
     
    length(common_cir_between_dis) ## 107 common circuits
      
    disease_limma_common <- lapply(diseases_limma_pathvals, function(x){x[x$Circuit %in% common_cir_between_dis, ]})
    
    disease_limma_unique <- lapply(diseases_limma_pathvals, function(x){x[!x$Circuit %in% common_cir_between_dis, ]})
    
   
    ## Get the pathways that are unique and common per disease
    
    pathways_common <- str_split(common_cir_between_dis, ":") %>% sapply(., function(x){x[1]}) %>% unique(.)
    length(pathways_common) ## 29 common pathways
    
    pathways_perDisease <- lapply(disease_limma_unique, function(x){
      
      str_split(x$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.)
      
    })
    
    pathways_common2 <- pathways_common[pathways_common %in% pathways_deregAll]
    
    
### 7. Run Hipathia on each cluster vs control ####
    
    ## LIMMA FOR EACH CLUSTER
    
    ### First we divide the matrix for each disease and controls and then paste the controls to each disease for the comparisons
    
    target_matrixes_clusters <- list()
    
    for (i in names(table(target_matrix$cluster))) {
      
      target_matrixes_clusters[[i]] <- target_matrix[grep(pattern = i, x = target_matrix$cluster), ] 
      
    }
    
    controls <- target_matrixes_clusters[["CONTROLS"]]
    
    target_matrixes_clusters <- target_matrixes_clusters[-5] ### Delete controls data to paste them later in each matrix
    
    target_matrixes_clusters <- lapply(target_matrixes_clusters, function (x) {rbind(x, controls)}) ## Pasting the controls to the diaseases
    
    
    ## Apply  the function to calculate dereg circuits for each cluster vs control comparison 
 
    diseases_limma_pathvals <- lapply(target_matrixes_clusters, function(x){get_diff_expressed_pathways_cluster(x, path_vals_phy)}) 
    saveRDS(cluster_limma_pathvals, file = here("rds", "clusters_limma_pathvals.rds"))
    
    pathways_deP_SADs <- lapply(diseases_limma_pathvals, function(x){ str_split(x$Circuit, ":") %>% sapply(., function(x){x[1]}) %>% unique(.)})
    
    