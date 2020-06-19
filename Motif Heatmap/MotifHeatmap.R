library(Seurat)
library(Signac)
library(viridis)

# Funtion to simplify isolating the EC clusters
GetGroups <-function(
  object,
  group.by,
  idents
) {
  if (is.null(x = group.by)) {
    obj.groups <- Idents(object = object)
  } else {
    obj.md <- object[[group.by]]
    obj.groups <- obj.md[, 1]
    names(obj.groups) <- rownames(x = obj.md)
  }
  if (!is.null(idents)) {
    obj.groups <- obj.groups[obj.groups %in% idents]
  }
  return(obj.groups)
}

# Load the object and only keep the EC clusters
object <- readRDS("tempp.rds")
object$cell_id <- object@active.ident
group.by <- "cell_id"
obj.groups <- GetGroups(
  object = object,
  group.by = group.by,
  idents = NULL
)
obj.groups.levels <- levels(obj.groups)
obj.groups.levels.EC <- obj.groups.levels[grepl("EC", obj.groups.levels)]
object <- subset(object, idents = obj.groups.levels.EC)

# Objects which will be used later
cell.id = obj.groups.levels.EC
matrix_list <- 1
top_genes = c()

# For each EC cluster...
for(idx in c(1:length(cell.id))) {
  message(paste0('- Getting fold values for ', cell.id[idx]))

  # Find the motif peaks relative to the other EC clusters
  da_peaks <- FindMarkers(
    object = object,
    ident.1 = cell.id[idx],
    only.pos = TRUE,
    test.use = 'LR',
    latent.vars = 'nCount_peaks'
  )

  # Test the differentially accessible peaks for overrepresented motifs
  enriched.motifs <- FindMotifs(
    object = object,
    features = head(rownames(da_peaks), 1000)
  )
  
  # Remove P-Value < 0.05
  enriched.motifs <- enriched.motifs[enriched.motifs[, 7] < 0.05, ]
  
  # Sort by P-Value and then motif FC
  enriched.motifs <- enriched.motifs[order(enriched.motifs[, 7], -enriched.motifs[, 6]), ]
  
  # Edit dataframe
  row.names(enriched.motifs) <- enriched.motifs$motif.name
  
  # Keep top 20 motifs
  top_genes <- c(top_genes, row.names(head(enriched.motifs, n = 20)))
  enriched.motifs <- enriched.motifs[, c(6, 8)]
  names(enriched.motifs)[1] <- cell.id[idx]
  
  # Merge all dataframes together
  if(matrix_list == 1){
    matrix_list <- enriched.motifs
  } else {
    matrix_list <- merge(matrix_list, enriched.motifs)
  }
}

# Display the head of the compiled list
head(matrix_list)

# Remove repeating motifs
top_genes_unique <- unique(top_genes)

# Only keep significant motifs
matrix_list_unique <- matrix_list[matrix_list$motif.name %in% top_genes_unique,]

# Sort by EC1, then EC2, ect.
matrix_list_unique <- matrix_list_unique[order(
  -matrix_list_unique[, 2], 
  -matrix_list_unique[, 3], 
  -matrix_list_unique[, 4], 
  -matrix_list_unique[, 5], 
  -matrix_list_unique[, 6], 
  -matrix_list_unique[, 7], 
  -matrix_list_unique[, 8], 
  -matrix_list_unique[, 9]
), ]

# Make motif names the row names to make the matrix numeric
row.names(matrix_list_unique) <- matrix_list_unique$motif.name
matrix_list_unique <- matrix_list_unique[, -c(1)]

# Make the dataframe a matrix and reverse it
matrix <- data.matrix(matrix_list_unique)
matrix <- apply(matrix, 2, rev)

# Produce heatmap
heatmap <- heatmap(matrix, Rowv=NA, Colv=NA, col = viridis::plasma(256), scale="column", margins=c(0,3))
