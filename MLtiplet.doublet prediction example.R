# MLtuplet R script v1.0 (last updated 19/02/2021)
# Author: Rachael Bashford-Rogers (Wellcome Centre for Human Genetics, University of Oxford, rbr1@well.ox.ac.uk)
# This script will predict scRNA-seq doublets/multiplets based on the features of known doublets/multiplet features (from VDJ-seq and/or CITE-seq information)

##### before loading R, ensure that all dependancies are correct
# module purge
# module load HDF5/1.10.5-gompi-2019a
# module load umap-learn/0.3.10-foss-2019a-Python-3.7.2
# module load Seurat/3.1.2-foss-2019a-R-3.6.0
# module load Harmony/1.0.0-foss-2019a-R-3.6.0
# R

library("Seurat")
library('harmony') 
library(ggplot2)
library(pryr)
library(future) 


concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha)) }

GetMitRibRatio <- function(object) {
  mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = object[["RNA"]]), value = TRUE)
  mito <- grep(pattern = "^MT-", x = rownames(x = object[["RNA"]]), value = TRUE)
  mitribcounts<- object@assays$RNA@counts[which(rownames(object@assays$RNA@counts) %in% mitrib), ]
  mitribcounts <- Seurat::CreateSeuratObject(mitribcounts)
  n.calc <- Seurat:::CalcN(object = mitribcounts)
  names(x = n.calc) <- paste(names(x = n.calc), "RNA", sep = '_mitrib_')
  mitoribo_ratio <- colSums(x = as.data.frame(Seurat::GetAssayData(object = mitribcounts, assay = "RNA", slot = "counts")[mito, , drop = FALSE]))/
    mitribcounts[[paste0("nCount_RNA")]]
  return(mitoribo_ratio)
}

Get_cell_type_score<-function(data, cell_type_label = "cell_type_overall"){
	Idents(object = data) <- cell_type_label
	cell_type = Idents(object = data) 
	cell_types = sort(unique(data@meta.data$cell_type_overall))
	subset_ids = NULL
	for(c in c(1:length(cell_types))){
		w = which(cell_type == cell_types[[c]])
		w = rownames(data@meta.data)[w]
		if(length(w)>400){
			w = sample(w, 400)
		}
		subset_ids = c(subset_ids, w)
	}
	data.small <- subset(x = data, cells = subset_ids)
	data.markers <- FindAllMarkers(data.small, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25,test.use = "roc")
	markers_use = NULL
	for(c in c(1:length(cell_types))){
		w = data.markers[which(data.markers[,"cluster"]==cell_types[c]),]
		w1 = w[rev(order(w[,"myAUC"])),"gene"][c(1:10)]
		data <- AddModuleScore(object = data,features = list(w1),ctrl = 5,name = concat(c("ModuleScore_",cell_types[c])))
		print(c)
	}
	head(data@meta.data)
	marker_array = data@meta.data[,grep("ModuleScore", colnames(data@meta.data))]
	return(marker_array)
}

Get_normalised_variables<-function(data, variate_names){
	variables = data@meta.data[, variate_names]
	variables_norm = variables
	library(compositions)
	idents = sort(unique(data@meta.data[,"orig.ident"]))
	for(i in c(1:length(variables[1,]))){
		x = variables[,i]
		names(x) = rownames(data@meta.data)
		for(s in c(1:length(idents))){
			w = rownames(data@meta.data)[which(data@meta.data[,"orig.ident"]== idents[s])]
			if(length(w)>0){
				clr1 = clr(x[w])
				clr1  = scale(as.numeric(clr1), center = TRUE, scale = TRUE)
				variables_norm[w,i] = as.numeric(clr1)
				print(concat(c(i, " ", idents[s] )))
			}}
	}
	return(variables_norm)
}


Predict_doublets<-function(marker_array, data, double_pos, single_pos, variate_names ){
	variables_norm = Get_normalised_variables(data, variate_names)
	ids = rownames(data@meta.data)
	variables_norm = cbind(variables_norm[ids,], marker_array[ids,])

	######### predict doublets
	cell_ids = rownames(variables_norm)
	doublet_probabilities = rep(-1, length(cell_ids))
	names(doublet_probabilities) = cell_ids
	doublet_probabilities_without_sample = doublet_probabilities
		
	summary(variables_norm[double_pos, variate_names[1]])
	summary(variables_norm[single_pos,variate_names[1]])
	### exclude outliers from doublet detection
	
	exclude_doublets = NULL
	exclude_singlets = NULL
	for(h in c(1:length(variables_norm[1,]))){
		var = variables_norm[double_pos,h]
		exclude_doublets  = c(exclude_doublets , double_pos [c(which(var<mean(var)-2*sd(var)), which(var>mean(var)+2*sd(var)))])
		var = variables_norm[single_pos,h]
		exclude_singlets  = c(exclude_singlets , single_pos [c(which(var<mean(var)-2*sd(var)), which(var>mean(var)+2*sd(var)))])
	}
	exclude_doublets = unique(exclude_doublets)
	exclude_singlets = unique(exclude_singlets)
	double_pos1 = setdiff(double_pos, exclude_doublets)
	single_pos1 = setdiff(single_pos, exclude_singlets)
		
	summary(variables_norm[double_pos1,variate_names[1]])
	summary(variables_norm[single_pos1,variate_names[1]])
	
	cells_sub = c(double_pos1, single_pos1)
		
	## model 0
	library(MASS)
	variables_norm1 = variables_norm
	x = data.frame(variables_norm1[cells_sub,]) ### subset of data
	x_all = data.frame(variables_norm1) ### all data
	x1=x
	class = c(rep("doublet", length(double_pos1)), rep("singlet", length(single_pos1)))
	x1$class = factor(class)
	colnames(x_all) = colnames(x1)[which(colnames(x1)!="class")]
	x = x1[,which(colnames(x1)!="class")]
	
	model <- glm(class ~ ., x1,family=binomial(link='logit'))
	# model <- glm(class ~ ., x1,family= quasibinomial)
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
	
	summary(variables_norm[which(glm.pred =="doublet"),variate_names[1]])
	summary(variables_norm[which(glm.pred =="singlet"),variate_names[1]])	
	print (perc_correct)
	
	doublet_probabilities_without_sample = p
	doublet_probabilities_with_sample = p
	doublet_probabilities_with_sample[double_pos1] = 1
	print(table(glm.pred)*100/length(glm.pred))
	
	doublet_predictions = c(list(prediction_table0), list(doublet_probabilities_without_sample), list(doublet_probabilities_with_sample))
	names(doublet_predictions) = c("prediction_table0",  "doublet_probabilities_without_sample","doublet_probabilities_with_sample")
	return(doublet_predictions)
}


############ example on human PBMC samples
batch = "CITE_CC1"
out_dir = "CITE_CC/"
out_dir_raw = "CITE_CC/"
PLOTS = "PLOTS/"


concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}

############# load scRNA-seq data
pbmc  = readRDS(file=concat(c(out_dir,PLOTS,"/Seurat_harmonised_pre_predicate_", batch,".pbmc")))### rds object of seurat RNA-seq object
pbmc = UpdateSeuratObject(pbmc)


############# get features for learning doublet/multiplet features
cell_ids = rownames(pbmc@meta.data)
scale_data = pbmc@assays$ RNA@ scale.data
raw_data = pbmc@assays$ RNA@ counts
genes = rownames(scale_data)
orig.ident = pbmc@meta.data$orig.ident
names(orig.ident) = cell_ids
cell_type = pbmc@meta.data$cell_type
names(cell_type) = cell_ids
raw_data = raw_data[,colnames(scale_data)]
cluster = cell_type
clusters = sort(unique(cluster))
orig.idents = sort(unique(orig.ident))

data = pbmc
############### ad mito-ribo ratio to meta data
mitoribo_ratio = GetMitRibRatio(data)
data <- AddMetaData(object = data, metadata = as.data.frame(mitoribo_ratio) ,col.name = "mito.ribo_ratio")

######################## doublets from BCR/TCR information
VDJ_doublets_summary  = readRDS(file=concat(c(out_dir,PLOTS,"/Seurat_VDJ_doublets_", batch,".rds")))
names(VDJ_doublets_summary[["VDJ_intra_BCR_TCR_doublets"]])
use = c("2x IGHs", "2x TRBs","2x IGH and 2x IGK/L","2x TRAs and 2x TRBs")
VDJ_doublets = sort(unique(c(unlist(VDJ_doublets_summary[["VDJ_BCR_TCR_doublets"]]), unlist(VDJ_doublets_summary[["VDJ_intra_BCR_TCR_doublets"]][use]))))

######################## doublets from CITE-seq information
list_doublets_CITE = readRDS(file=concat(c(out_dir,PLOTS,"/Seurat_CITE_doublets_", batch,".rds")))
CITE_doublets = sort(unique(unlist(list_doublets_CITE)))
######################## doublets from DoubletFinder
list_doublets_DoubletFinder = readRDS(file=concat(c(out_dir,PLOTS,"/Seurat_harmonised_pre_predicate_", batch,".DoubletFinder1")))
DoubletFinder_doublets = names(which(list_doublets_DoubletFinder=="Doublet"))

combined = unique(c(VDJ_doublets, CITE_doublets, DoubletFinder_doublets))
combined1 = unique(c(VDJ_doublets, CITE_doublets))
training_doublets_list = c(list(VDJ_doublets), list(CITE_doublets), list(DoubletFinder_doublets), list(combined1), list(combined))
names(training_doublets_list) = c("VDJ-training","CITE-training","DoubletFinder-training","VDJ & CITE-training", "VDJ, CITE & DoubletFinder-training")
saveRDS(file=concat(c(out_dir,PLOTS,"/Seurat_training_", batch,".rds")), training_doublets_list)

############################ get module scores for each cell type
marker_array = Get_cell_type_score(data, cell_type_label = "cell_type_overall")

##################### PREDICT_DOUBLETS #########################
####### example run a titration of training sets
fraction_training = c(1:5)*2/10
cell_ids = rownames(data@meta.data)
for(typ in c(1:length(training_doublets_list))){
	doublets_training = training_doublets_list[[typ]]
	for(f in c(1:length(fraction_training))){
		doublets_training_use = doublets_training
		if(fraction_training[f]<1){
			n = ceiling(length(doublets_training)* fraction_training[f])
			doublets_training_use = sample(doublets_training,n)
		}
		single_pos = setdiff(cell_ids,doublets_training_use)
		double_pos = doublets_training_use
		#### build classifier based on "variate_names" plus "marker_array' (module scores for each cell type), and run the classifier over the whole dataset
		doublet_predictions = Predict_doublets(marker_array, data, double_pos, single_pos, variate_names = c("nUMI","mito.ribo_ratio"))
		saveRDS(file=concat(c(out_dir,PLOTS,"/Seurat_training_", batch,"_", typ,"_", f,".rds")), doublet_predictions)
	}	
}

##### get characteristics of each prediction 
cell_types = sort(unique(cell_type))
count_predicted_doublets = NULL
names_pred_doublets = NULL
for(typ in c(1:length(training_doublets_list))){
	doublets_training = training_doublets_list[[typ]]
	for(f in c(1:length(fraction_training))){
		doublet_predictions = readRDS(file=concat(c(out_dir,PLOTS,"/Seurat_training_", batch,"_", typ,"_", f,".rds")))
		prob = doublet_predictions $ doublet_probabilities_without_sample
		predict_doublets = names(prob[which(prob<0.5)])
		t1 = table(cell_type[predict_doublets])
		t2 = table(orig.ident[predict_doublets])
		t = rep(0,length(c("total", cell_types, orig.idents)))
		names(t) = c("total",cell_types, orig.idents)
		t[names(t1)] = t1
		t[names(t2)] = t2
		t["total"] = ceiling(length(doublets_training)* fraction_training[f])
		names_pred_doublets = c(names_pred_doublets, concat(c( names(training_doublets_list) [typ], ".", fraction_training[f],"x" )))
		if(length(count_predicted_doublets)==0){count_predicted_doublets = t
		}else{count_predicted_doublets = rbind(count_predicted_doublets, t)}
	}	
}
rownames(count_predicted_doublets) = names_pred_doublets
filtering_file = concat(c(out_dir_raw,"/Doublet_predictions_titration_", batch,".txt"))
write.table(count_predicted_doublets, file = filtering_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
print(concat(c("scp -p mfj169@rescomp2.well.ox.ac.uk:", filtering_file," ./ " )))


##################### PREDICT_DOUBLETS #########################
####### example run a only one sample for the training sets
cell_ids = rownames(data@meta.data)
for(typ in c(1:length(training_doublets_list))){
	doublets_training = training_doublets_list[[typ]]
	for(f in c(1:length(orig.idents))){
		w = names(which(orig.ident== orig.idents[f]))
		doublets_training_use = intersect(doublets_training, w)
		single_pos = setdiff(cell_ids,doublets_training_use)
		double_pos = doublets_training_use
		#### build classifier based on "variate_names" plus "marker_array' (module scores for each cell type), and run the classifier over the whole dataset
		doublet_predictions = Predict_doublets(marker_array, data, double_pos, single_pos, variate_names = c("nUMI","mito.ribo_ratio"))
		saveRDS(file=concat(c(out_dir,PLOTS,"/Seurat_training_", batch,"_", typ,"_", orig.idents[f],".rds")), doublet_predictions)
	}	
}


##### get characteristics of each prediction 
cell_types = sort(unique(cell_type))
count_predicted_doublets = NULL
names_pred_doublets = NULL
for(typ in c(1:length(training_doublets_list))){
	doublets_training = training_doublets_list[[typ]]
	for(f in c(1:length(orig.idents))){
		doublet_predictions = readRDS(file=concat(c(out_dir,PLOTS,"/Seurat_training_", batch,"_", typ,"_", orig.idents[f],".rds")))
		prob = doublet_predictions $ doublet_probabilities_without_sample
		predict_doublets = names(prob[which(prob<0.5)])
		t1 = table(cell_type[predict_doublets])
		t2 = table(orig.ident[predict_doublets])
		t = rep(0,length(c("total", cell_types, orig.idents)))
		names(t) = c("total",cell_types, orig.idents)
		t[names(t1)] = t1
		t[names(t2)] = t2
		t["total"] = ceiling(length(doublets_training)* fraction_training[f])
		names_pred_doublets = c(names_pred_doublets, concat(c( names(training_doublets_list) [typ], ".", orig.idents[f],"x" )))
		if(length(count_predicted_doublets)==0){count_predicted_doublets = t
		}else{count_predicted_doublets = rbind(count_predicted_doublets, t)}
	}	
}
rownames(count_predicted_doublets) = names_pred_doublets
filtering_file = concat(c(out_dir_raw,"/Doublet_predictions_per_sample_training_", batch,".txt"))
write.table(count_predicted_doublets, file = filtering_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
