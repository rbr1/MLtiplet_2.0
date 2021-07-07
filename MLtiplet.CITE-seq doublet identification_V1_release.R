# RECITublet R script v1.0 (last updated 02/09/2020)
# Author: Rachael Bashford-Rogers (Wellcome Centre for Human Genetics, University of Oxford, rbr1@well.ox.ac.uk)
# This script will identify scRNA-seq droplets that have features of doublets/multiplets based on CITE-seq information

# run as an interactive job via: qlogin -pe shmem 4 -P immune-rep.prjc -q short.qc
# module loadings: 
# qlogin -pe shmem 4 -P immune-rep.prjc -q short.qc
module purge
module load HDF5/1.10.5-gompi-2019a
module load umap-learn/0.3.10-foss-2019a-Python-3.7.2
module load Seurat/3.1.2-foss-2019a-R-3.6.0
module load Harmony/1.0.0-foss-2019a-R-3.6.0
R

#### identifying VDJ doublets/multiplets
#### load requirements
library("Seurat")
library('harmony') 
library(ggplot2)


concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}

########################################################
#### load in object with gene expression and ADT data
#### example is run of 3 healthy PBMC samples (from 10X website)
batch = "CC2" #### batch name for all outputs
pbmc = readRDS(file=concat(c("Seurat_object_CITE_CC2.seurat"))) #### downloadable from: https://www.dropbox.com/s/tuwamdw4zz5mswb/Seurat_object_CITE_CC2.seurat?dl=0

### GEX QC (standard)
QC_GEX<-function(object, prefix_figures){
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = object@assays$ RNA@ counts), value = TRUE)
	length(mito.genes)
	percent.mito <- Matrix::colSums(object@assays$ RNA@ counts[mito.genes, ]) / Matrix::colSums(object@assays$ RNA@ counts)
	object <- AddMetaData(object = object,metadata = percent.mito,col.name = "percent.mito")
	
	matrix = object@assays$ RNA@ counts
	mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = matrix), value = TRUE, ignore.case = TRUE)
	mito <- grep(pattern = "^MT-", x = rownames(x = matrix), value = TRUE, ignore.case = TRUE)
	mitribcounts<- matrix[which(rownames(matrix) %in% mitrib), ]
	mitoribo_ratio <- Matrix::colSums(mitribcounts[mito, , drop = FALSE])/Matrix::colSums(mitribcounts)
	object <- AddMetaData(object = object, metadata = as.data.frame(mitoribo_ratio) ,col.name = "mito.ribo_ratio")
	
	fileout1=concat(c(prefix_figures,"1.pdf"))
	w=5
	pdf(file=fileout1, height=w*3, width=w*5)
	par(mfrow= c(1,1), mar = c(8,4,4,4))
	VlnPlot(object,features = c("nFeature_RNA", "nCount_RNA","nCount_ADT", "percent.mito"),ncol = 1, group.by  = "orig.ident" )
	dev.off()

	object <- subset(object, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mito < 0.20)
	table(object@ meta.data$orig.ident)
	
	object <- NormalizeData(object)
	object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
	object <- ScaleData(object = object)
	return(object)		
}

prefix_figures = concat(c("Seurat_filtering_", batch,"_"))
pbmc = QC_GEX(pbmc, prefix_figures)

######### batch correct samples and UMAP
Batch_corrections_and_dimensionality_reduction<-function(object, prefix_figures){
	object <- RunPCA(object, verbose = FALSE)
	object <- RunHarmony(object, c("orig.ident"),plot_convergence = TRUE, nclust = 60, max.iter.cluster = 100, max.iter.harmony = 10)

	library(cowplot)
	fileout1=concat(c(prefix_figures,"2.pdf"))
	w=6
	pdf(file=fileout1, height=w*1, width=w*3)
	par(mfrow= c(1,1), mar = c(4,4,4,4))	
	p1 <- DimPlot(object = object, reduction = "harmony", pt.size = .1, group.by = "orig.ident", do.return = TRUE)
	p2 <- VlnPlot(object = object, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
	plot_grid(p1,p2)
	dev.off()

	object = RunUMAP(object ,reduction = "harmony", dims = 1:20)
	object = FindNeighbors(object ,reduction = "umap", dims = 1:2)
	object = FindClusters(object,reduction = "umap" ,resolution = 0.1)
	
	fileout1=concat(c(prefix_figures,"3.pdf"))
	w=6
	pdf(file=fileout1, height=w*1, width=w*2.6)
	par(mfrow= c(1,1), mar = c(4,4,4,4))
	p1 <- DimPlot(object = object, reduction = "umap", pt.size = .1, group.by = "orig.ident", do.return = TRUE)
	p2 <- DimPlot(object = object, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, do.return = TRUE)
	plot_grid(p1,p2)
	dev.off()
	return(object)
}

pbmc = Batch_corrections_and_dimensionality_reduction(pbmc, prefix_figures)

##### load in labels for cell types
cell_type = readRDS(file=concat(c("Seurat_harmonised_pre_predicate_CITE_CC1.all_cell_types"))) #### downloadable from: https://www.dropbox.com/s/uswu0jdnbgrycof/Seurat_harmonised_pre_predicate_CITE_CC1.all_cell_types?dl=0
cell_id = rownames(pbmc@meta.data)
cell_id_ct = names(cell_type)
samples = sort(unique(pbmc@meta.data$orig.ident))
for(i in c(1:length(samples))){
	cell_id_ct = gsub(concat(c("||", samples[i])), concat(c("_",i)), cell_id_ct, fixed = T)
}
length(intersect(cell_id, cell_id_ct))
length(setdiff(cell_id ,cell_id_ct))
cell_type_object = rep("-", length(cell_id))
names(cell_type_object) = cell_id
cell_type_object[cell_id_ct] = cell_type
pbmc@meta.data$cell_type= cell_type_object[cell_id]

fileout1=concat(c(prefix_figures,"4.pdf"))
w=6
pdf(file=fileout1, height=w*1, width=w*2.1)
par(mfrow= c(1,1), mar = c(4,4,4,4))
p1 <- DimPlot(object = pbmc, reduction = "umap", group.by = "cell_type", pt.size = .1, do.return = TRUE)
plot_grid(p1)
dev.off()

saveRDS(file=concat(c("Seurat_object_", batch,".seurat_post_filtering")), pbmc)

########################################################
######### normalise ADT
pbmc = readRDS(file=concat(c("Seurat_object_", batch,".seurat_post_filtering")))

Normalise_ADT<-function(object){
	object <- NormalizeData(object, assay = "ADT", normalization.method = "CLR")
	object <- ScaleData(object, assay = "ADT")
	return(object)}
	
pbmc=Normalise_ADT(pbmc)
	
############# CITE-seq - gene matching files
CITE_seq = rownames(pbmc@assays$ ADT@counts)
ref_file = "CITE_Seq_gene_name.txt" #### downloadable from: https://www.dropbox.com/s/rqdazl5lhrtqwgf/CITE_Seq_gene_name.txt?dl=0
### file containing the mapping of CITE-seq names with gene/protein names
p <- as.matrix(read.csv(ref_file, head=T, sep="\t"))
gene_CITE =p[match(CITE_seq, p[,1]),2]
names(CITE_seq) = gene_CITE
CITE_seq_protein = p[,3]
names(CITE_seq_protein) = p[,1]
gene_CITEs = sort(unique(gene_CITE))
orig.ident = pbmc@meta.data$orig.ident
orig.idents = sort(unique(orig.ident))
CITE_data_all =pbmc@assays$ ADT@scale.data
GEX_data_all =pbmc@assays$ RNA@counts

############# functions: threshold CITE-seq data to determine which cells are expressing each marker by taking the best CITE-seq divider between GEX postivie and GEX negative.
Plot_CITE_seq_versus_gene_expression<-function(gene_CITE, gene_CITEs, CITE_data_all,orig.ident, orig.idents, CITE_seq,  GEX_data_all, prefix_figures){
	for(i in c(1:length(CITE_seq))){
		fileout1=concat(c(prefix_figures, CITE_seq[i],".pdf"))
		library(RColorBrewer)
		cols =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
		w=3.25
		pdf(file=fileout1, height=w*5, width=w*6)
		par(mfrow= c(5,6), mar = c(5,5,3,3))
		for(s in c(1:length(orig.idents))){
			w = which(orig.ident== orig.idents[s])
			x = CITE_data_all [CITE_seq[i], w]
			y = GEX_data_all[gene_CITE[i], w]
			w1 = intersect(which(is.na(x)==F),which(is.na(y)==F))
			x = x[w1]
			y = y[w1]				
			if(length(x)>100){
				print (length(x))
				plot (y, x, col = cols[1], pch = 21, bg =  cols[1],main = concat(c(orig.idents[s],"\n", CITE_seq[i])), xlab = gene_CITE[i], ylab =  CITE_seq[i])
			}}
		dev.off()
	}
}
Get_mean_gene_expression<-function(gene_CITE, cluster, gene_CITEs, GEX_data_all){
	clusters = sort(unique(cluster))
	m_mean_expression = t(matrix(data = 0,nrow = length(clusters),ncol = length(gene_CITEs), dimnames = c(list(clusters),list(gene_CITEs))))
	m_quantile_expression = t(matrix(data = 0,nrow = length(clusters),ncol = length(gene_CITEs), dimnames = c(list(clusters),list(gene_CITEs))))
	for(c in c(1:length(clusters))){
		w = which(cluster== clusters[c])
		m_mean_expression[,clusters[c]] = apply(GEX_data_all[gene_CITEs,w], 1, mean)
		m_quantile_expression[,clusters[c]] = apply(GEX_data_all[gene_CITEs,w], 1, function(x){quantile(x,prob = 0.9)})
		print(c)	}
	CITE_gene_expression = c(list(m_mean_expression), list(m_quantile_expression))
	names(CITE_gene_expression) = c("mean", "90th quantile")
	return(CITE_gene_expression)
}
CITE_seq_level_per_sample<-function(orig.ident, orig.idents, CITE_seq,  cluster, CITE_data_all){
	clusters = sort(unique(cluster))
	m_CITE_expression_by_sample = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	m_CITE_expression_by_cluster = t(matrix(data = 0,nrow = length(clusters),ncol = length(CITE_seq), dimnames = c(list(clusters),list(CITE_seq))))
	non_na = NULL
	for(i in c(1:length(CITE_seq))){
		non_na = c(non_na, list(which(is.na(CITE_data_all[CITE_seq[i],])==F)))}
	names(non_na) = CITE_seq
	for(c in c(1:length(orig.idents))){
		w = which(orig.ident== orig.idents[c])
		print (c)
		for(i in c(1:length(CITE_seq))){
			w1 = intersect(non_na[[i]],w )
			m_CITE_expression_by_sample[CITE_seq[i], orig.idents[c]] = mean(CITE_data_all[CITE_seq[i],w1])
		}}
	for(c in c(1:length(clusters))){
		print(c)
		w = which(cluster== clusters[c])
		for(i in c(1:length(CITE_seq))){
			w1 = intersect(non_na[[i]],w )
			m_CITE_expression_by_cluster[CITE_seq[i],clusters[c]] = mean(CITE_data_all[CITE_seq[i],w1])
		}}
	CITE_seq_expression = c(list(m_CITE_expression_by_sample), list(m_CITE_expression_by_cluster))
	names(CITE_seq_expression) = c("m_CITE_expression_by_sample", "m_CITE_expression_by_cluster")
	return(CITE_seq_expression)
}
Threshold_CITE_seq_positivity<-function(orig.idents, orig.ident, CITE_seq, m_CITE_expression_by_sample, CITE_data_all, CITE_gene_expression, cell_ids){
	library(MASS)
	thresholds = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	prediction_accuracy = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	pos_accuracy = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	neg_accuracy = t(matrix(data = 0,nrow = length(orig.idents),ncol = length(CITE_seq), dimnames = c(list(orig.idents),list(CITE_seq))))
	cell_ids_CITE_GEX = intersect(cell_ids, colnames(CITE_data_all))
	for(i in c(1:length(CITE_seq))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample[CITE_seq[i],]!=0))
		quantile = CITE_gene_expression $ '90th quantile'[names(CITE_seq)[i],]
		quantile = quantile[which(is.na(quantile)==F)]
		mean = CITE_seq_expression$m_CITE_expression_by_cluster[CITE_seq[i],]
		mean = mean[which(is.na(mean)==F)]
		pos = intersect(names(which(quantile>= quantile(quantile, 0.5))), names(which(mean >= quantile(mean, 0.5))))
		neg = setdiff(names(which(quantile <=  quantile(quantile, 0.1))), pos)
		w_use_na = names(which(is.na(CITE_data_all[CITE_seq[i],])==F))
		# neg = setdiff(clusters, pos)
		for(s in c(1:length(sample_use))){
			w = cell_ids [which(orig.ident== sample_use[s])]
			w1 = intersect(w,  cell_ids [which(cluster %in% neg)])
			w2 = intersect(w,  cell_ids [which(cluster %in% pos)])
			w1 = intersect(intersect(w1, cell_ids_CITE_GEX), w_use_na)
			w2 = intersect(intersect(w2, cell_ids_CITE_GEX), w_use_na)
			if(min(c(length(w1), length(w2)))>5){
				x_pos = CITE_data_all[CITE_seq[i],w2]
				x_neg = CITE_data_all[CITE_seq[i],w1]
				x = c(x_neg, x_pos)
				if(length(x)>10){
					if(length(unique(x))>2){
						x1 = c(x_neg, x_pos)
						class = c(rep("neg", length(x_neg)), rep("pos", length(x_pos)))
						ids = c(w1,w2)
						class1 = data.frame(ids )
						class1$class = factor(class)
						class1$x = x
						model <- lda(class ~x, data = class1)
						predictions <- model %>% predict(class1)
						if(min(table(predictions $class))>5){
							x2 = data.frame(seq(from = max(x[which(predictions $class=="neg")]), to = min(x[which(predictions $class=="pos")]), length = 10000))
							colnames(x2) = 'x'
							pred = predict(model, newdata = x2)$class
							thr = mean(c(min(x2[which(pred=="pos"),]), max(x2[which(pred=="neg"),]))) 
							acc = length(which(predictions$class==class))*100/length(class)
							pacc = length(intersect(which(class=="pos"),which(x1>= thr)))*100/length(which(class=="pos"))
							nacc = length(intersect(which(class=="neg"),which(x1< thr)))*100/length(which(class=="neg"))
							thresholds[CITE_seq[i],sample_use[s]] =  thr
							prediction_accuracy [CITE_seq[i],sample_use[s]] = acc
							pos_accuracy [CITE_seq[i],sample_use[s]] =  pacc
							neg_accuracy [CITE_seq[i],sample_use[s]] =  nacc
							print(i)
			}}}}
		}
	}
	CITE_seq_positivity_thresholds = c(list(thresholds), list(prediction_accuracy), list(pos_accuracy), list(neg_accuracy))
	names(CITE_seq_positivity_thresholds) = c(("thresholds"),("prediction_accuracy"),("pos_accuracy"),("neg_accuracy"))
	return(CITE_seq_positivity_thresholds)
}
Choose_CITE_seq_markers<-function(CITE_seq, CITE_seq_positivity_thresholds, Tpos_accuracy, Tneg_accuracy, Tprediction_accuracy){
	headers = c("thresholds", "prediction_accuracy", "pos_accuracy", "neg_accuracy")
	thrs = matrix(data = 0,nrow = length(CITE_seq),ncol = length(headers), dimnames = c(list(CITE_seq),list(headers)))
	for(i in c(1:length(CITE_seq))){
		w1 = intersect(which(CITE_seq_positivity_thresholds$thresholds[CITE_seq[i],]!=0), which(is.finite(CITE_seq_positivity_thresholds$thresholds[i,])))
		w2 = intersect(which(CITE_seq_positivity_thresholds$prediction_accuracy[i,]!=0), which(is.finite(CITE_seq_positivity_thresholds$prediction_accuracy[i,])))
		w3 = intersect(which(CITE_seq_positivity_thresholds$pos_accuracy[i,]!=0), which(is.finite(CITE_seq_positivity_thresholds$pos_accuracy[i,])))
		w4 = intersect(which(CITE_seq_positivity_thresholds$neg_accuracy[i,]!=0), which(is.finite(CITE_seq_positivity_thresholds$neg_accuracy[i,])))
		thrs[i,] = c(mean(CITE_seq_positivity_thresholds$thresholds[i,w1]), mean(CITE_seq_positivity_thresholds$prediction_accuracy[i,w2]), mean(CITE_seq_positivity_thresholds$pos_accuracy[i,w3]), mean(CITE_seq_positivity_thresholds$neg_accuracy[i,w4]))
	}
	CITE_seq_use = CITE_seq[intersect (intersect(which(thrs[,'pos_accuracy']> Tpos_accuracy), which(thrs[,'neg_accuracy']> Tneg_accuracy)), intersect(which(is.finite(thrs[,'thresholds'])==T), which(thrs[,'prediction_accuracy']> Tprediction_accuracy)))]
	info_CITE_seq_use = thrs [CITE_seq_use,]
	CITE_seq_choose = c(list(thrs), list(CITE_seq_use), list(info_CITE_seq_use))
	names(CITE_seq_choose) = c("thrs","CITE_seq_use","info_CITE_seq_use")
	return(CITE_seq_choose)
}
Plot_CITE_seq_thresholds<-function(CITE_seq_choose, prefix_figures, CITE_seq_positivity_thresholds, CITE_data_all, GEX_data_all){
	library(RColorBrewer)
	cols =  rep(add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5), 100)
	w=4
	fileout1 = concat(c(prefix_figures,"1.pdf"))
	pdf(file=fileout1, height=w*4, width=w*4)
	par(mfrow= c(5,5), mar = c(5,5,3,3))
	sample_match = match(orig.ident, orig.idents)
	for(i in c(1:length(CITE_seq_choose$CITE_seq_use))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample[CITE_seq_choose$CITE_seq_use[i],]!=0))
		if(length(sample_use)>0){
			for(s in c(1:length(sample_use))){
				w = which(orig.ident== sample_use[s])
				x = GEX_data_all[names(CITE_seq_choose$CITE_seq_use)[i],w]
				y = CITE_data_all[(CITE_seq_choose$CITE_seq_use[i]),w]
				w1 = intersect(which(is.na(x)==F),which(is.na(y)==F))
				x = x[w1]
				y = y[w1]
				if(length(unique(y))>10){
					plot (x, y, col = cols[sample_match[w]], pch = 21, bg =  cols[sample_match[w]],main = concat(c(sample_use[s],"\n", CITE_seq_choose$CITE_seq_use[i])), 
					xlab = names(CITE_seq_choose$CITE_seq_use)[i], ylab = (CITE_seq_choose$CITE_seq_use[i]))
					segments(-100, CITE_seq_positivity_thresholds $thresholds[CITE_seq_choose$CITE_seq_use[i], sample_use[s]], 100000, CITE_seq_positivity_thresholds $thresholds[CITE_seq_choose$CITE_seq_use[i], sample_use[s]], lwd = 2, col = "red")
				}}
			}
		print (i)
		}
	dev.off()
}
Score_CITE_positivity<-function(CITE_seq_choose, CITE_seq_positivity_thresholds, CITE_data_all, cell_ids){
	CITE_score = (CITE_data_all*NA)
	for(i in c(1:length(CITE_seq_choose$CITE_seq_use))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample [CITE_seq_choose$CITE_seq_use[i],]!=0))
		if(length(sample_use)>0){
			for(s in c(1:length(sample_use))){
				w = cell_ids [which(orig.ident== sample_use[s])]
				threshold = CITE_seq_positivity_thresholds $thresholds[CITE_seq_choose$CITE_seq_use[i], sample_use[s]]
				w1 = intersect(w, names (which(CITE_data_all[CITE_seq_choose$CITE_seq_use[i],]> threshold)))
				if(length(w1)>0){
					x = CITE_data_all[CITE_seq_choose$CITE_seq_use[i],w1]
					#score = rank(x)/length(x)
					score = (x-min(x))
					score = (score/max(score))+1
					CITE_score[CITE_seq_choose$CITE_seq_use[i],w1] = score}}}
		print (i)
		}
	return(CITE_score)}
Plot_CITE_scores<-function(CITE_score, CITE_data_all, CITE_seq_choose, orig.ident, orig.idents, prefix_figures,CITE_seq_expression){
	library(RColorBrewer)
	cols =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	w=4 
	fileout1 = concat(c(prefix_figures,"2.txt"))
	pdf(file=fileout1, height=w*4, width=w*4)
	par(mfrow= c(5,5), mar = c(5,5,3,3))
	sample_match = match(orig.ident, orig.idents)
	for(i in c(1:length(CITE_seq_choose$CITE_seq_use))){
		sample_use = names(which(CITE_seq_expression$m_CITE_expression_by_sample[CITE_seq_choose$CITE_seq_use[i],]!=0))
		if(length(sample_use)>0){
			for(s in c(1:length(sample_use))){
				w = which(orig.ident== sample_use[s])
				scores = CITE_score[CITE_seq_choose$CITE_seq_use[i],w]
				scaled_expression = CITE_data_all[CITE_seq_choose$CITE_seq_use[i],w]
				scores[which(is.na(scores))] = 0
				plot (scaled_expression, scores, col = cols[1], pch = 21, bg =  cols[1],main = concat(c(sample_use[s],"\n", CITE_seq_choose$CITE_seq_use[i])))
				}}
		print (i)
		}
	dev.off()
}

############# run thresholding functions
prefix_figures = concat(c("CITE_expression", batch,"_"))
Plot_CITE_seq_versus_gene_expression(gene_CITE, gene_CITEs, CITE_data_all,orig.ident, orig.idents, CITE_seq,  GEX_data_all, prefix_figures)

cell_type = pbmc@meta.data$cell_type
cluster = pbmc@meta.data$seurat_clusters
CITE_seq_expression  = CITE_seq_level_per_sample(orig.ident, orig.idents, CITE_seq, cluster= cluster, CITE_data_all)
CITE_gene_expression = Get_mean_gene_expression(gene_CITE, cluster, gene_CITEs, GEX_data_all)

############# predict thresholds of positive CITE seq cells
cell_ids = rownames(pbmc@meta.data)
CITE_seq_positivity_thresholds = Threshold_CITE_seq_positivity(orig.idents, orig.ident, CITE_seq, CITE_seq_expression, CITE_data_all, CITE_gene_expression,cell_ids)

############# choose CITE seq markers that have good predictions
CITE_seq_choose = Choose_CITE_seq_markers(CITE_seq, CITE_seq_positivity_thresholds,Tpos_accuracy = 20, Tneg_accuracy = 30, Tprediction_accuracy = 60)

############# plot CITE-seq thresholds
Plot_CITE_seq_thresholds(CITE_seq_choose, prefix_figures, CITE_seq_positivity_thresholds, CITE_data_all, GEX_data_all)

############# generate CITE-seq score per positive cell and save output

CITE_score = Score_CITE_positivity(CITE_seq_choose, CITE_seq_positivity_thresholds, CITE_data_all, cell_ids)
saveRDS(file=concat(c("CITE_seq_scale.CITE_score")), CITE_score)

############# plot scores versus expression of CITE seq
prefix_figures = concat(c("CITE_expression", batch,"_"))
Plot_CITE_scores(CITE_score, CITE_data_all, CITE_seq_choose, orig.ident, orig.idents, prefix_figures,CITE_seq_expression)


############# identify doublets/multiplets based on mutually exclusive markers 
CITE_score = readRDS(file=concat(c("CITE_seq_scale.CITE_score")))
CITE_score = t(CITE_score) [cell_ids,]
CITE_genes = colnames(CITE_score)
names(CITE_genes) = names(CITE_seq)[match(CITE_genes , CITE_genes)]

ref_file = "CITE_seq_mutually_exclusive_expression.txt" ### from dropbox link: https://www.dropbox.com/s/8p7qq9o8qqrk83y/CITE_seq_mutually_exclusive_expression.txt?dl=0
p <- as.matrix(read.csv(ref_file, head=T, sep="\t"))
p=p[which(p[,3] %in% c("Mutually exclusive","Weak")),]
p=p[intersect(which(p[,1] !="CD56"),which(p[,2] !="CD56")), ]
p=p[which(p[,1] %in% CITE_seq_protein ==T),]
p=p[which(p[,2] %in% CITE_seq_protein ==T),]
############# check that all markers are in the CITE-seq list
p[,1][which(p[,1] %in% CITE_seq_protein ==T)]
p[,2][which(p[,2] %in% CITE_seq_protein ==T)]

mutually_exclusive_CITE_seq = NULL
mutually_exclusive_CITE_seq_names = NULL
for(i in c(1:length(p[,1]))){
	mutual1 = names(CITE_seq_protein) [which(CITE_seq_protein == p[i,1])]
	mutual2 = names(CITE_seq_protein) [which(CITE_seq_protein == p[i,2])]
	for(i1 in c(1:length(mutual1))){
		for(i2 in c(1:length(mutual2))){
			if(concat(c(mutual2[i2],":", mutual1[i1])) %in% mutually_exclusive_CITE_seq_names==F){
				mutually_exclusive_CITE_seq = c(mutually_exclusive_CITE_seq, list(c(mutual1[i1], mutual2[i2])))
				mutually_exclusive_CITE_seq_names=c(mutually_exclusive_CITE_seq_names, concat(c(mutual1[i1],":", mutual2[i2])))
}}}}

############# normalise the variables that will be used in the predictor
variables = pbmc@meta.data[,c("percent.mito", "nFeature_RNA","nCount_RNA","mito.ribo_ratio","nCount_ADT")]
rownames(variables)= cell_ids
variables_norm = variables
library(compositions)
for(i in c(1:length(variables[1,]))){
	x = variables[,i]
	names(x) = cell_ids
	for(s in c(1:length(orig.idents))){
		w = cell_ids[which(orig.ident== orig.idents[s])]
		if(length(w)>0){
			clr1 = x[w]/mean(x[w])
			variables_norm[w,i] = as.numeric(clr1)
			print(concat(c(i, " ", orig.idents[s] )))
		}}}

metaD = variables_norm

############# functions: identify double positive cells
Boxplot_nUMI_ngene_mt_cells1<-function(doublets_list, fileout1,metaD){
	headers = colnames(metaD)
	all_cells = rownames(metaD)
	library(RColorBrewer)
	cols = add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	w=2.5
	pdf(file=fileout1, height=w*1.6*3, width=w*4*1.8)
	par(mfrow= c(3,2), mar = c(18,4.5,3,1))
	for(h in c(1:length(headers))){
		groups = NULL
		for(i in c(1:length(doublets_list))){
			x = metaD[doublets_list[[i]], headers[h]]
			groups = c(groups, list(metaD[doublets_list[[i]], headers[h]]))}
		factors = names(doublets_list)
		boxplot(groups, ylab = "", col = cols,names= factors, las = 2, main= headers[h], outline =T, border = "grey",cex.lab = 0.9, lwd = 0.5, cex.main = 0.7)
	}
	dev.off()
}
Identify_CITE_seq_doublets <-function(fileout1, orig.idents, orig.ident, mutually_exclusive_CITE_seq_names, CITE_genes1, CITE_score, metaD){
	library(MASS)
	library(RColorBrewer)
	cols =  add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	variables = metaD
	variables_names = colnames(variables)
	doublets_list_all = list()
	for(s in c(1:length(orig.idents))){
		w = cell_ids [which(orig.ident== orig.idents[s])]
		double_pos = NULL
		single_pos = NULL
		double_pos_classification = NULL
		sample_positive_CITE = names(which(apply(CITE_score[w,],2,function(x){length(which(is.na(x)==F))})>0))
		for(cc in c(1:length(mutually_exclusive_CITE_seq))){
			cc1 = mutually_exclusive_CITE_seq[[cc]][1]
			cc2 = mutually_exclusive_CITE_seq[[cc]][2]
			cite_sub = CITE_genes[which(CITE_genes %in%(c(cc1,cc2)))]
			# plot_CITE = names(which(rowSums(raw_data[names(cite_sub),w])>0))
			plot_CITE = intersect(sample_positive_CITE,(cite_sub))
			if(length(plot_CITE)==2){
				x = CITE_score[w,plot_CITE[1]]
				y = CITE_score[w,plot_CITE[2]]
				w1 = w[intersect(which(is.finite(x)),which(is.finite(y)))]
				w1 = w1[intersect(which(x[w1]>0.01), which(y[w1]>0.01))]
				w2 = w[setdiff(which(is.finite(x)),which(is.finite(y)))]
				w3 = w[setdiff(which(is.finite(y)),which(is.finite(x)))]
				w23 = sort(unique(c(w2,w3)))
				if(min(c(length(w1),length(w23)))>2){
					double_pos = c(double_pos, w1)
					single_pos = c(single_pos, w23)
					double_pos_classification = c(double_pos_classification, rep(mutually_exclusive_CITE_seq_names[cc], length(w1)))
				}}
				}
		single_pos = setdiff(w, double_pos)
		doublets_list = NULL
		names = mutually_exclusive_CITE_seq
		var = variables[,"nCount_RNA"]
		names(var) = rownames(variables)
		summary = NULL
		for(cc in c(1:length(mutually_exclusive_CITE_seq))){
			x = var [double_pos[which(double_pos_classification== mutually_exclusive_CITE_seq_names[cc])]]
			doublets_list = c(doublets_list, list(names(x)))
			if(length(summary)==0){summary = summary(x)
			}else{summary = rbind(summary, summary(x))}
		}
		x = var [single_pos]
		summary = rbind(summary, summary(x))
		doublets_list = c(doublets_list, list(single_pos))		
		rownames(summary) = c(mutually_exclusive_CITE_seq_names,"singlets")
		names(doublets_list)= c(mutually_exclusive_CITE_seq_names,"singlets")
		w = intersect(which(summary[,'Median']<1.05*summary['singlets','Median']), which(rownames(summary)!="singlets"))
		w1 = intersect(which(summary[,'Median']<1.05*summary['singlets','Mean']), which(rownames(summary)!="singlets"))
		w = sort(unique(c(w,w1)))
		exclude_CITE = rownames(summary)[w]
		doublets_list[w] = list(NULL) 
		######################## plot
		names(doublets_list) = gsub("_TotalSeqC","",names(doublets_list) )
		fileout1 =concat(c("Seurat_CITE_nGene_CITEseq_doublets_", type,"_", batch,"_", orig.idents[s],".pdf"))
		Boxplot_nUMI_ngene_mt_cells1(doublets_list, fileout1,metaD)		
		w = which(double_pos_classification %in% exclude_CITE==F)
		doublet_info = c(list(cbind(double_pos[w], double_pos_classification[w])),list(single_pos), list(summary))
		names(doublet_info) = c("true_doublets", "singlets","summary_nUMIs")
		doublets_list_all = c(doublets_list_all, list(doublet_info))
	}
	names(doublets_list_all) = orig.idents
	return(doublets_list_all)
}
############# identify double positive cells
prefix_figures = concat(c("CITE_expression_", batch,"_"))
fileout1 = concat(c(prefix_figures,"3.txt"))
doublets_list_all =Identify_CITE_seq_doublets(fileout1, orig.idents, orig.ident, mutually_exclusive_CITE_seq_names,CITE_genes1, CITE_score, metaD)

types = NULL
for(i in c(1:length(doublets_list_all))){types = sort(unique(c(types, doublets_list_all[[i]][[1]][,2])))}
list_doublets_CITE = rep(list(c()), length(types))
names(list_doublets_CITE) = types
all_CITE_doublets = NULL
for(i in c(1:length(doublets_list_all))){
	if(length(doublets_list_all[[i]][[1]])>0){
		all_CITE_doublets = c(all_CITE_doublets, doublets_list_all[[i]][[1]][,1])}
	for(t in types){
		w = which(doublets_list_all[[i]][[1]][,2]==t)
		list_doublets_CITE[[t]] = c(list_doublets_CITE[[t]], doublets_list_all[[i]][[1]][w,1])}}
all_CITE_doublets= sort(unique(all_CITE_doublets))

saveRDS(file=concat(c("CITE_doublets_", batch,".rds")), all_CITE_doublets)

############# plot nGenes, nUMIs and mt% between these methods
Boxplot_nUMI_ngene_mt_doublets<-function(doublets_list, fileout1,metaD){
	headers = colnames(metaD)
	all_cells = rownames(metaD)
	library(RColorBrewer)
	cols = add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	w=2.5
	pdf(file=fileout1, height=w*1.5, width=w*4*0.9)
	par(mfrow= c(1,5), mar = c(18,4.5,3,1))
	for(h in c(1:length(headers))){
		groups = NULL
		names = NULL
		factors = c(names(doublets_list), "other")
		for(i in c(1:length(doublets_list))){
			x =metaD[doublets_list[[i]], headers[h]]
			if(length(x)>2){
				names = c(names, factors[i])
				groups = c(groups, list(x))}}
		x1 = unlist(groups)
		groups= c(groups, list(metaD[setdiff(all_cells, unlist(doublets_list)), headers[h]]))
		names = c(names, "other")
		x2 = metaD[setdiff(all_cells, unlist(doublets_list)), headers[h]]
		pval = wilcox.test(x1,y=x2,alternative = "two.sided")$p.value
		pval1 = signif(pval,digits =3)
		if(pval<1e-10){pval1 = "<1e-10"}
		boxplot(groups, ylab = "", col = cols,names= names, las = 2, main= concat(c(headers[h],"\np-value ",pval1)), outline =T, border = "grey",cex.lab = 0.9, lwd = 0.5, cex.main = 0.7, ylim = c(0,max(unlist(groups))))
	}
	dev.off()
}

fileout1 = concat(c("CITE_expression_", batch,"_doublets.pdf"))
metaD = variables_norm
names(list_doublets_CITE) = gsub("_TotalSeqC","",names(list_doublets_CITE))
Boxplot_nUMI_ngene_mt_doublets(list_doublets_CITE,  fileout1,metaD)

######################## Write output
saveRDS(file="Seurat_CITE_doublets.rds", list_doublets_CITE)


