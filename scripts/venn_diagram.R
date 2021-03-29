#################################################
#### Venn diagrams of the data for report  ######
################################################

library(pacman)
pacman::p_load("here", "hipathia", "biomaRt", "utils", "stringr", "RColorBrewer",
               "AnnotationDbi", "tidyr","data.table", "xlsx", "gplots")

# Load the data
 diseases_limma_pathvals <- readRDS(file = file = here("rds", "diseases_limma_pathvals.rds"))
 circuits_perdisease <- lapply(diseases_limma_pathvals, function(x){x$Circuit})

# Chart

 myCol <- brewer.pal(5, "Set2") 
 
venn.diagram(
  x = circuits_perdisease[c(1,2,5,6,7)],
  category.names = names(circuits_perdisease)[c(1,2,5,6,7)],
  filename = 'venn_diagramm_12-567.png',

    # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1300 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.085, 0.085, 0.085, 0.085, 0.085),
  cat.fontfamily = "sans",
  # rotation = 1 only for 3 bent diagrams
)


# gplots::venn(data = circuits_perdisease[1:5])
