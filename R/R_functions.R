
#' Mitochondrial-ribosomal RNA ratio
#'
#' This function calculates the mitochondrial-ribosomal RNA ratio. It assumes that the
#' contains the mitichondrial gene IDs begin with the prefix "MT-" and that the 
#' ribosomal gene IDs  begin with the prefix "RP-". If this is not the case, then
#' the gene IDs can be specified.
#'
#' @param object = seurat object
#' @param mitrib = ribosomal and mitichondrial gene IDs
#' @param mito = mitichondrial gene IDs
#' @return An array of mitochondrial-ribosomal RNA ratios per cell.
#' @export
GetMitRibRatio <- function(object, mitrib=NULL, mito = NULL){
  if(length(mitrib)==0){mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = object[["RNA"]]), value = TRUE)}
  if(length(mito)==0){mito <- grep(pattern = "^MT-", x = rownames(x = object[["RNA"]]), value = TRUE)}
  mitribcounts<- object@assays$RNA@counts[which(rownames(object@assays$RNA@counts) %in% mitrib), ]
  mitribcounts <- Seurat::CreateSeuratObject(mitribcounts)
  n.calc <- Seurat:::CalcN(object = mitribcounts)
  names(x = n.calc) <- paste(names(x = n.calc), "RNA", sep = '_mitrib_')
  mitoribo_ratio <- colSums(x = as.data.frame(Seurat::GetAssayData(object = mitribcounts, assay = "RNA", slot = "counts")[mito, , drop = FALSE]))/mitribcounts[[paste0("nCount_RNA")]]
  return(mitoribo_ratio)
}

#' Cell type score
#'
#' This function calculates the module score for each cluster/cell type, using 
#' the Seurat AddModuleScore function. 
#'
#' @param object = seurat object
#' @param cell_type_label = the label for the cell types/cluster IDs to be used for calculating the module score
#' @return A matrix of the module scores per cell types/cluster per cell.
#' @export
Get_cell_type_score<-function(object, cell_type_label = "cell_type_overall"){
  Idents(object = object) <- cell_type_label
  cell_type = Seurat::Idents(object = object) 
  cell_types = sort(unique(object@meta.data$cell_type_overall))
  subset_ids = NULL
  for(c in c(1:length(cell_types))){
    w = which(cell_type == cell_types[[c]])
    w = rownames(data@meta.data)[w]
    if(length(w)>400){
      w = sample(w, 400)
    }
    subset_ids = c(subset_ids, w)
  }
  data.small <- subset(x = object, cells = subset_ids)
  data.markers <- Seurat::FindAllMarkers(data.small, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25,test.use = "roc")
  markers_use = NULL
  for(c in c(1:length(cell_types))){
    w = data.markers[which(data.markers[,"cluster"]==cell_types[c]),]
    w1 = w[rev(order(w[,"myAUC"])),"gene"][c(1:10)]
    data <- Seurat::AddModuleScore(object = object,features = list(w1),ctrl = 5,name = concat(c("ModuleScore_",cell_types[c])))
  }
  marker_array = object@meta.data[,grep("ModuleScore", colnames(object@meta.data))]
  return(marker_array)
}

#' Get normalised variables
#'
#' This function calculates the normalised variables for inclusion into the doublet detection 
#' model. This performs log-centred normalisation per sample
#'
#' @param object = seurat object
#' @param variate_names = the list of labels (from the meta.data of the seurat object)
#' @return A matrix of the normalised variables per cell.
#' @export
Get_normalised_variables<-function(object, variate_names){
  variables = object@meta.data[, variate_names]
  variables_norm = variables
  library(compositions)
  idents = sort(unique(object@meta.data[,"orig.ident"]))
  for(i in c(1:length(variables[1,]))){
    x = variables[,i]
    names(x) = rownames(object@meta.data)
    for(s in c(1:length(idents))){
      w = rownames(object@meta.data)[which(object@meta.data[,"orig.ident"]== idents[s])]
      if(length(w)>0){
        clr1 = clr(x[w])
        clr1  = scale(as.numeric(clr1), center = TRUE, scale = TRUE)
        variables_norm[w,i] = as.numeric(clr1)
        print(concat(c(i, " ", idents[s] )))
      }}
  }
  return(variables_norm)
}

#' Predict doublets
#'
#' This function builds classifier based on "variate_names" plus "marker_array' (module scores 
#' for each cell type), and runs the classifier over the whole dataset. This assumes that the 
#' doublets_training are all cell IDs wihtin the seurat object. 
#'
#' @param object = seurat object
#' @param marker_array = a matrix of the module scores per cell types/cluster per cell.
#' @param doublets_training = an array of the training doublet/multiplet set (e.g. from VDJ-seq, CITE-seq or other).
#' @param variables_norm = a matrix of the normalised variables per cell.
#' @param training_trim = a parameter for filtering outliers from the training doublet/multiplet set (default = 1).
#' @return doublet object (list of 1. the training doublet prediction table, 2. the doublet probabilities excluding the 
#' @return training doublets, 3. the doublet probabilities including the training doublets, 4. the list of predicted doublets). 
#' @export
Predict_doublets<-function(object ,marker_array, doublets_training, variables_norm, training_trim = 1){
  ids = rownames(object@meta.data)
  variables_norm1 = cbind(variables_norm[ids,], marker_array[ids,])
  double_pos = doublets_training
  single_pos = setdiff(ids , doublets_training_use)
  
  doublet_probabilities = rep(-1, length(cell_ids))
  names(doublet_probabilities) = cell_ids
  doublet_probabilities_without_sample = doublet_probabilities
  ### exclude outliers from doublet detection
  exclude_doublets = NULL
  exclude_singlets = NULL
  for(h in c(1:length(variables_norm[1,]))){
    var = variables_norm[double_pos,h]
    exclude_doublets  = c(exclude_doublets , double_pos [c(which(var<mean(var)-training_trim*sd(var)), which(var>mean(var)+ training_trim*sd(var)))])
    var = variables_norm[single_pos,h]
    exclude_singlets  = c(exclude_singlets , single_pos [c(which(var<mean(var)-training_trim*sd(var)), which(var>mean(var)+ training_trim*sd(var)))])
  }
  exclude_doublets = unique(exclude_doublets)
  exclude_singlets = unique(exclude_singlets)
  double_pos1 = setdiff(double_pos, exclude_doublets)
  single_pos1 = setdiff(single_pos, exclude_singlets)
  
  cells_sub = c(double_pos1, single_pos1)
  x = data.frame(variables_norm1[cells_sub,]) ### subset of data
  x_all = data.frame(variables_norm1) ### all data
  x1=x
  class = c(rep("doublet", length(double_pos1)), rep("singlet", length(single_pos1)))
  x1$class = factor(class)
  colnames(x_all) = colnames(x1)[which(colnames(x1)!="class")]
  x = x1[,which(colnames(x1)!="class")]
  
  model <- glm(class ~ ., x1,family=binomial(link='logit'))
  p <- predict(model, newdata=x, type="response")
  glm.pred=rep("doublet", length(p))
  glm.pred[p>0.5]="singlet"
  length(which(glm.pred==class))*100/length(p)
  
  prediction_table0 = table(glm.pred, class) #### columns = train class, rows = pred class
  
  p <- predict(model, newdata=x_all, type="response")
  glm.pred=rep("doublet", length(p))
  glm.pred[p>0.5]="singlet"
  pred_doublet = cell_ids[which(glm.pred=="doublet")]
  pred_singlet = cell_ids[which(glm.pred=="singlet")]
  perc_correct =length(which(double_pos  %in% pred_doublet))*100/length(double_pos)
  print (perc_correct)
  
  doublet_probabilities_without_sample = p
  doublet_probabilities_with_sample = p
  doublet_probabilities_with_sample[doublets_training] = 0
  print(table(glm.pred)*100/length(glm.pred))
  predicted_doublets = names(which(doublet_probabilities_with_sample<0.5))

  doublet_predictions = c(list(prediction_table0), list(doublet_probabilities_without_sample), list(doublet_probabilities_with_sample),list(predicted_doublets))
  names(doublet_predictions) = c("prediction_table",  "doublet_probabilities_excluding_training","doublet_probabilities_including_training","predicted_doublets")
  return(doublet_predictions)
}

