#################################################################
#### Histograms of the distribution of relevant circuits ########
################################################################

library(pacman)

pacman::p_load("here", "hipathia", "utils", "stringr", "genefilter", "tidyr","data.table","limma", "ggplot2", "patchwork")


### Read Hipathia results from the validation dataset (2 dataset send by Guillermo) ####

  results <- readRDS(file = "/home/m3m/INFO_PROYECTO/autoinmunes/rds/validation_cohort_Hiresults.rds")

  metadata <- fread(file = here("data", "feedback_guillermo_Marzo2021", "metadata_reformatted.tsv"))
  
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
  
  ## LIMMA FOR ALL SAMPLES 
  
  target_matrix <- metadata
  
  sample_group_limma <- target_matrix[, c(1,4,5,7,9,10,11)] # Note that for Limma the Control has to be first in the levels
  
  
  design <- model.matrix( ~ 0  + type + age + gender + center + pool + RIN ,  data = sample_group_limma)
  colnames(design) <- gsub("typeCTRL", "C", colnames(design)) 
  colnames(design) <- gsub("typeDISEASE", "disease", colnames(design)) 
  rownames(design) <- sample_group_limma$samples
  cont.matrix <- makeContrasts( diseasevsC=disease-C, levels=design)
  
  fit <- lmFit(path_vals_phy , design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit3 <- eBayes(fit2)
  
  tableLimma <- topTable(fit3, number = rownames(path_vals_phy), adjust.method="bonferroni", sort.by="p")
  top_LimmaALL <- tableLimma[tableLimma$adj.P.Val < 0.05,] 
  top_LimmaALL <- top_LimmaALL[order(top_LimmaALL$adj.P.Val, decreasing = F),]
  
  saveRDS(top_LimmaALL, file = here("rds", "top_LimmaALL_validationCohort_allvsC.rds"))

  anot_all <- data.frame(Circuit = get_path_names(pathways, rownames(top_LimmaALL)),
                        logFC = top_LimmaALL$logFC,
                        P.val = top_LimmaALL$P.Value,
                        FDR.Pval = top_LimmaALL$adj.P.Val,
                        UniprotKB = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "uniprot" ,collapse = T),
                        GO = get_pathways_annotations(rownames(top_LimmaALL), pathways, dbannot= "GO" ,collapse = T), 
                        stringsAsFactors = F)
  
  ##### SELECT DEREGULATED CIRCUITS AND CONTROLS FOR HISTOGRAMS DISTRIBUTIONS ##########
  
  idx1 <- which(rownames(path_vals_phy) %in% rownames(top_LimmaALL))
  idx2 <- which(colnames(path_vals_phy) %in% metadata$samples[metadata$diagnosis == "CTRL"])
  
  path_vals_hist <- path_vals_phy[idx1, idx2 ]
  all(rownames(path_vals_hist) %in% rownames(top_LimmaALL))
  
  # rownames(path_vals_hist) <- get_path_names(pathways, rownames(path_vals_hist))
  
  top20_min20 <- rownames(top_LimmaALL)[c(1:20, 314:333)] ## Para ver solo los 20 mas dereg y menos dereg
  path_vals_hist <- path_vals_hist[top20_min20,]
  
  # mean_order <- data.frame(cir = colnames(path_vals_hist),mean = apply(path_vals_hist, 2, mean))
  # mean_order <- mean_order[order(mean_order$mean, decreasing = F), ]
  
  # path_vals_hist <- path_vals_hist[, mean_order$cir]
  numSamples <- nrow(path_vals_hist)
  colors <- rainbow(numSamples)

  dev.new()
  boxplot(t(path_vals_hist),
          main = "Deregulated circuits distributions in Control samples",
          col =  colors,
          cex.axis = 0.7,
          las = 2 )
  

  
  sub1 <- colnames(path_vals_hist)[colnames(path_vals_hist) %in% metadata$samples[metadata$diagnosis == "CTRL"]] 
  sub2 <- colnames(path_vals_hist)[colnames(path_vals_hist) %in% metadata$samples[metadata$diagnosis != "CTRL"]] 
  

  dev.print(svg, file = here("results","boxplotDeregulatedCir_controls_topmin20_norm.svg"))
  dev.print(svg, file = here("results","boxplotDeregulatedCir_controls.svg"))
  dev.print(svg, file = here("results","boxplotDeregulatedCir_controls_ordered.svg"))
  
  
  ##### SELECT DEREGULATED CIRCUITS AND DISEASE FOR HISTOGRAMS DISTRIBUTIONS ##########
  
  idx3 <- which(rownames(path_vals_phy) %in% rownames(top_LimmaALL))
  idx4 <- which(colnames(path_vals_phy) %in% metadata$samples[metadata$diagnosis != "CTRL"])
  
  path_vals_hist_dis <- path_vals_phy[idx3, idx4 ]
  all(rownames(path_vals_hist_dis) %in% rownames(top_LimmaALL))
  
  # rownames(path_vals_hist_dis) <- get_path_names(pathways, rownames(path_vals_hist_dis))
  top20_min20 <- rownames(top_LimmaALL)[c(1:20, 314:333)] ## Para ver solo los 20 mas dereg y menos dereg
  path_vals_hist_dis <- path_vals_hist_dis[top20_min20,]
  # path_vals_hist_dis <- data.frame(t(path_vals_hist_dis), stringsAsFactors = F)
  
  
  # mean_order <- data.frame(cir = colnames(path_vals_hist_dis),mean = apply(path_vals_hist_dis, 2, mean))
  # mean_order <- mean_order[order(mean_order$mean, decreasing = F), ]
  
  ## Ximo quiere que use el orden de menor a mayor de los controles
  # path_vals_hist_dis <- path_vals_hist_dis[, mean_order$cir]
  
  numSamples <- nrow(path_vals_hist_dis)
  colors <- rainbow(numSamples)
  
  dev.new()
  boxplot(t(path_vals_hist_dis),
          main = "Deregulated circuits distributions in Disease samples",
          col =  colors,
          cex.axis = 0.7,
          las = 2 )
  
  dev.print(svg, file = here("results","boxplotDeregulatedCir_diseases_topmin20_norm.svg"))
  # dev.print(svg, file = here("results","boxplotDeregulatedCir_disease.svg"))
  dev.print(svg, file = here("results","boxplotDeregulatedCir_disease_ordered.svg"))
  
  
  # ##### MAKE A FUNCTION: SELECT DEREGULATED CIRCUITS AND DISEASE FOR HISTOGRAMS DISTRIBUTIONS ##########
  # 
  # select_plot_circuits_samples <- function(path_vals_phy, diagnosis ,  top_LimmaALL = top_LimmaALL, metadata = metadata, metaginfo= pathways){
  # 
  #   idx1 <- which(rownames(path_vals_phy) %in% rownames(top_LimmaALL))
  #   idx2 <- which(colnames(path_vals_phy) %in% metadata$samples[metadata$diagnosis != diagnosis])
  # 
  #   path_vals_hist_dis <- path_vals_phy[idx1, idx2 ]
  #   all(rownames(path_vals_hist_dis) %in% rownames(top_LimmaALL))
  # 
  #   rownames(path_vals_hist_dis) <- get_path_names(metaginfo, rownames(path_vals_hist_dis))
  #   path_vals_hist_dis <- data.frame(t(path_vals_hist_dis), stringsAsFactors = F)
  # 
  #   mean_order <- data.frame(cir = colnames(path_vals_hist_dis),mean = apply(path_vals_hist_dis, 2, mean))
  #   mean_order <- mean_order[order(mean_order$mean, decreasing = F), ]
  # 
  #   path_vals_hist_dis <- path_vals_hist_dis[, mean_order$cir]
  # 
  # numSamples <- ncol(exprs(rawData))
  # colors <- rainbow(numSamples)  
  #   p <- boxplot(path_vals_hist_dis,
  #                main = paste0("Deregulated circuits distributions in ", diagnosis, " samples"),
  #                col =  colors,
  #                cex.axis = 0.7,
  #                las = 2 )
  # 
  #   png(filename = here("results", paste0("boxplotDeregulatedCir_",diagnosis,".png")), width = 5000, height = 3000, res = 200)
  #   plot(p)
  #   dev.off()
  # 
  # }
  # 
  # select_plot_circuits_samples(path_vals_phy = path_vals_phy , diagnosis = "MCTD")
  # 
  # sapply( names(table(metadata$diagnosis)), function(x) {select_plot_circuits_samples(path_vals_phy = path_vals_phy, diagnosis = x)})

  ###### DO ALL BOXPLOTS OF TOP DEREG CIRCUITS ONE BESIDES ANOTHER ####
  
  dat_c <- stack(as.data.frame(t(path_vals_hist)))
  dat_c$group <- "Control"
  
  # disease_noC2 <- path_vals_hist_dis[, colnames(path_vals_hist_dis) %in% metadata$samples[metadata$cluster != 2]]
  dat_d <- stack(as.data.frame(t(path_vals_hist_dis)))
  # dat_d <- stack(as.data.frame(t(disease_noC2)))
  dat_d$group <- "Disease"
  
  dat <- rbind(dat_c,dat_d)
  colnames(dat) <- c("activation_values", "circuits", "Group")
  
  subset <- which(dat$circuits %in% rownames(top_LimmaALL)[c(1:40)])
  dat_plot <- dat[subset,]
  dat_plot <- dat_plot[order(dat_plot$activation_values, decreasing = F),]
  dat_plot$circuits <- factor(dat_plot$circuits, levels=unique(dat_plot$circuits))
  # dat_plot$circuits <- get_path_names(pathways, as.character(dat_plot$circuits))
  
  dev.new()
  ggplot(dat_plot, aes(x = circuits, y = activation_values, fill = Group)) + 
    geom_boxplot()+
    theme( axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size= 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0, size = 12),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18))
      # facet_wrap( ~ group, scales = 'free_x') ## ESTO SERIA PARA VARIOS PLOTS 
  
  ggsave(filename = here("results","boxplotDeregulatedCir_top40max_noC2.png"), plot = last_plot())
  ggsave(filename = here("results","boxplotDeregulatedCir_top40max.png"), plot = last_plot())

  