# RECITublet R script v1.0 (last updated 07/07/2021)
# Author: Rachael Bashford-Rogers (Wellcome Centre for Human Genetics, University of Oxford, rbr1@well.ox.ac.uk)
# This script will identify scRNA-seq droplets that have features of doublets/multiplets based on VDJ-seq information

#### identifying VDJ doublets/multiplets
#### load requirements
library("Seurat")

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

cell_ids = rownames(pbmc@meta.data)

######################## identify doublets from BCR/TCR information
Get_BCR_TCR_doublets<-function(m_VDJ, threshold_umis, cell_type){
	cell_types = sort(unique(cell_type))
	B_cell_types = cell_types[grep("B cell", cell_types)]
	T_cell_types = cell_types[grep("T cell", cell_types)]
	non_B_T_cell_types = setdiff(cell_types, c(B_cell_types, T_cell_types))
	wb1 = suppressWarnings (rownames(m_VDJ$BCR) [intersect(which(m_VDJ$BCR[,"constant_region1"]!='-'), which(as.numeric(m_VDJ$BCR[,"n_umis1"])>= threshold_umis))])
	wb2 = suppressWarnings (rownames(m_VDJ$BCR) [intersect(which(m_VDJ$BCR[,"constant_region2"]!='-'), which(as.numeric(m_VDJ$BCR[,"n_umis2"])>= threshold_umis))])
	wt1 = suppressWarnings (rownames(m_VDJ$TCR) [intersect(which(m_VDJ$TCR[,"constant_region1"]!='-'), which(as.numeric(m_VDJ$TCR[,"n_umis1"])>= threshold_umis))])
	wt2 = suppressWarnings (rownames(m_VDJ$TCR) [intersect(which(m_VDJ$TCR[,"constant_region2"]!='-'), which(as.numeric(m_VDJ$TCR[,"n_umis2"])>= threshold_umis))])
	wb12 = intersect(wb1,wb2)
	wt12 = intersect(wt1,wt2)
	wb= unique(c(wb1,wb2))
	wt= unique(c(wt1,wt2))
	high_confidence_clashes = intersect(intersect(wb1,wb2), intersect(wt1,wt2))### both chains present for TCR/BCR
	mid_confidence_clashes = unique(c(wb12 [which(wb12 %in% c(wt1,wt2))],wt12 [which(wt12 %in% c(wb1,wb2))]))
	lower_confidence_clashes = intersect(wb, wt)
	
	non_B_cell_BCRs = intersect(rownames(m_VDJ$BCR)[which(cell_type %in% B_cell_types==F)], c(wb1, wb2))
	non_B_cell_BCRs = non_B_cell_BCRs[which(as.numeric (m_VDJ$BCR[non_B_cell_BCRs,"n_umis1"])> threshold_umis)]
	non_T_cell_TCRs = intersect(rownames(m_VDJ$TCR)[which(cell_type %in% T_cell_types==F)], c(wt1, wt2))
	non_T_cell_TCRs = non_T_cell_TCRs[which(as.numeric(m_VDJ$TCR[non_T_cell_TCRs,"n_umis1"])> threshold_umis)]

	VDJ_BCR_TCR_doublets = c(list(high_confidence_clashes), list(mid_confidence_clashes), list(lower_confidence_clashes),list(non_B_cell_BCRs), list(non_T_cell_TCRs))
	names(VDJ_BCR_TCR_doublets) = c("IGH+IGK/L+TRA+TRB","3 of [IGH,IGK/L,TRA,TRB]", "IGH or IGK/L + TRA or TRB","non-B cell clustered + BCR", "non-T cell clustered + TCR" )
	return(VDJ_BCR_TCR_doublets)
}
#
Count_BCR_TCR_doublets_by_cell_type<-function(cell_type, cell_types, VDJ_BCR_TCR_doublets){
	headers = c("all","4 chains","3 chains","2 chains","non-cell-type clustered")
	counts = matrix(data = 0,nrow = length(cell_types),ncol = length(headers), dimnames = c(list(cell_types),list(headers)))
	t = table(cell_type)
	counts[names(t), "all"] = t
	t = table(cell_type[ VDJ_BCR_TCR_doublets $ "IGH+IGK/L+TRA+TRB"])
	counts[names(t), "4 chains"] = t
	t = table(cell_type[ VDJ_BCR_TCR_doublets $ "3 of [IGH,IGK/L,TRA,TRB]"])
	counts[names(t), "3 chains"] = t
	t = table(cell_type[ VDJ_BCR_TCR_doublets $ "IGH or IGK/L + TRA or TRB"])
	counts[names(t), "2 chains"] = t
	t = table(cell_type[ unique(c(VDJ_BCR_TCR_doublets $ "non-B cell clustered",VDJ_BCR_TCR_doublets $ "non-T cell clustered"))])
	counts[names(t), "non-cell-type clustered"] = t
	return(counts)}
#
Get_intra_BCR_TCR_doublets <-function(m_VDJ, threshold_umis){
	w = which(m_VDJ$BCR[,"mixed_contig_n_umis1"]!='-')
	Bchain1_doublets = rownames(m_VDJ$BCR[w,]) [which(as.numeric(m_VDJ$BCR[w,"mixed_contig_n_umis1"])> threshold_umis)]
	w = which(m_VDJ$BCR[,"mixed_contig_n_umis2"]!='-')
	Bchain2_doublets = rownames(m_VDJ$BCR[w,]) [which(as.numeric(m_VDJ$BCR[w,"mixed_contig_n_umis2"])> threshold_umis)]
	
	w = which(m_VDJ$TCR[,"mixed_contig_n_umis1"]!='-')
	Tchain1_doublets = rownames(m_VDJ$TCR[w,]) [which(as.numeric(m_VDJ$TCR[w,"mixed_contig_n_umis1"])> threshold_umis)]
	w = which(m_VDJ$TCR[,"mixed_contig_n_umis2"]!='-')
	Tchain2_doublets = rownames(m_VDJ$TCR[w,]) [which(as.numeric(m_VDJ$TCR[w,"mixed_contig_n_umis2"])> threshold_umis)]

	Tchain12_doublets = intersect(Tchain1_doublets, Tchain2_doublets)
	Bchain12_doublets = intersect(Bchain1_doublets, Bchain2_doublets)
	
	VDJ_intra_BCR_TCR_doublets = c(list(Bchain1_doublets), list(Tchain2_doublets),list(Bchain12_doublets), list(Tchain12_doublets))
	names(VDJ_intra_BCR_TCR_doublets) = c("2x IGHs","2x TRBs","2x IGH and 2x IGK/L","2x TRAs and 2x TRBs")
	return(VDJ_intra_BCR_TCR_doublets)
}
#
Count_intra_BCR_TCR_doublets_by_cell_type <-function(cell_type, cell_types, VDJ_intra_BCR_TCR_doublets){
	headers = c("all",names(VDJ_intra_BCR_TCR_doublets))
	counts = matrix(data = 0,nrow = length(cell_types),ncol = length(headers), dimnames = c(list(cell_types),list(headers)))
	t = table(cell_type)
	counts[names(t), "all"] = t
	for(i in c(1:length(names(VDJ_intra_BCR_TCR_doublets)))){
		t = table(cell_type[ VDJ_intra_BCR_TCR_doublets [[i]]])
		counts [names(t),names(VDJ_intra_BCR_TCR_doublets)[i]] = t
	}
	return(counts)}
#
Boxplot_nUMI_ngene_mt_doublets<-function(doublets_list, fileout1,metaD){
	headers = colnames(metaD)
	all_cells = rownames(metaD)
	library(RColorBrewer)
	cols = add.alpha (brewer.pal(8, "Dark2"), alpha = 0.5)
	w=2.5
	pdf(file=fileout1, height=w*1.1, width=w*4*0.9)
	par(mfrow= c(1,5), mar = c(12,4.5,3,1))
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

########################
######### load VDJ information
m_VDJ = readRDS(file = "Seurat_CITE_CC2_VDJ.rds")### rds object of BCR/TCR object, which can be downloaded from: https://www.dropbox.com/s/01re50jmhiz0isq/Seurat_CITE_CC2_VDJ.rds?dl=0

m_VDJ_BCR  = m_VDJ$BCR
m_VDJ_TCR  = m_VDJ$TCR

######################## normalise the variables that will be used in the predictor
cell_ids = rownames(pbmc@meta.data)
orig.ident = pbmc@meta.data$orig.ident
orig.idents = sort(unique(orig.ident))
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

######################## Identify and count doublets/multiplets from the VDJ-seq data
cell_type = pbmc@meta.data$cell_type
names(cell_type) = cell_ids
cell_types = sort(unique(cell_types))
### note that B cell subsets should include the term "B cell" in annotation, and T cell subsets should include the term "T cell" in annotation
VDJ_BCR_TCR_doublets = Get_BCR_TCR_doublets(m_VDJ, threshold_umis = 10, cell_type)
VDJ_intra_BCR_TCR_doublets = Get_intra_BCR_TCR_doublets(m_VDJ, threshold_umis = 10)
counts_BCR_TCR_doublets = Count_BCR_TCR_doublets_by_cell_type(cell_type, cell_types, VDJ_BCR_TCR_doublets)
counts_intra_BCR_TCR_doublets = Count_intra_BCR_TCR_doublets_by_cell_type(cell_type, cell_types, VDJ_intra_BCR_TCR_doublets)

doublets = c(VDJ_BCR_TCR_doublets, VDJ_intra_BCR_TCR_doublets)
metaD = variables_norm

######################## plot nGenes, nUMIs and mt-rb-ratio between these methods
fileout1 = concat(c("Seurat_CITE_nGene_VDJ_doublets_", batch,".pdf"))
Boxplot_nUMI_ngene_mt_doublets(doublets, fileout1,metaD)

fileout1 = concat(c(out_dir, PLOTS,"Seurat_CITE_nGene_VDJ_doublets_summary_", batch,".pdf"))
VDJ_doublets1 = c(list(sort(unique(unlist(doublets[c("IGH+IGK/L+TRA+TRB","3 of [IGH,IGK/L,TRA,TRB]", "IGH or IGK/L + TRA or TRB")])))),
list(sort(unique(unlist(doublets[c("2x IGHs")])))),
list(sort(unique(unlist(doublets[c("2x TRBs")])))),
list(doublets[["non-B cell clustered + BCR"]]),list(doublets[["non-T cell clustered + TCR"]]))
names(VDJ_doublets1) = c("BCR + TCR","multiple BCRs","multiple TCRs","non-B cell clustered + BCR","non-T cell clustered + TCR")
Boxplot_nUMI_ngene_mt_doublets(VDJ_doublets1,  fileout1,metaD)

######################## Write output
VDJ_doublets_summary = c(list(VDJ_BCR_TCR_doublets),list(VDJ_intra_BCR_TCR_doublets))
names(VDJ_doublets_summary) = c('VDJ_BCR_TCR_doublets', 'VDJ_intra_BCR_TCR_doublets')
saveRDS(file=concat(c("Seurat_VDJ_doublets_", batch,".rds")), VDJ_doublets_summary)

