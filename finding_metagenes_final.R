########################################################
## Author: Chaitanya R. Acharya
##
## R code 1) to process single cell RNAseq data
##        2) to compute enrichment scores
##        3) to derive modules
##        4) to perform functional annotation of modules
########################################################

library(cellrangerRkit);library(data.table)
source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

### Preprocess sc-RNAseq data
###
genome <- "GRCh38"
brca = load_cellranger_matrix("~/Research/BRCA/Loi_BCTCells_NatMed_071018/",genome=genome)
dim(exprs(brca))

brca_dat = as.matrix(exprs(brca))
brca_dat = data.frame(brca_dat,check.names=F)

probes = fread("~/Research/pathways/Homo_sapiens.GRCh38.79_table.txt",header=F,data.table=F)
colnames(probes)=c("ENSEMBL","Genes","Coordinates","Strand")
probes$Chr = gsub('\\:.*', '', probes$Coordinates)
probes=probes[c(which(probes$Chr %in% c(1:22)),grep("^X$",probes$Chr),grep("^Y$",probes$Chr)),]
probes = probes[-grep("^MIR",probes$Genes),]
probes = probes[-grep("^RPS",probes$Genes),]
probes = probes[-grep("^RPL",probes$Genes),]

keep_probes = intersect(rownames(brca_dat),probes$ENSEMBL)
seu_exprs = brca_dat[match(keep_probes,rownames(brca_dat)),]
probes_exprs = probes[match(keep_probes,probes$ENSEMBL),]

seu_exprs$Genes = as.factor(as.character(probes_exprs$Genes))
dat = setDT(seu_exprs)[, lapply(.SD, max), by = Genes]
dat = as.data.frame(dat)
rownames(dat) = dat[,1]
dat = dat[,-1]
fwrite(dat,"BRCA_scRNAseq_all_patients_gene_level.dat",sep="\t",col.names=T,row.names=T,quote=F)

pdata = data.frame(barcode=colnames(dat),patientid = c(rep(1,5174),rep(2,1137)))
rownames(pdata) = pdata$barcode

################################################################################################
# VST transform data
################################################################################################
library(monocle)
gene_df = data.frame(gene_short_name=rownames(dat))
rownames(gene_df) = rownames(dat)
fd = new("AnnotatedDataFrame", data = gene_df)
mycds = newCellDataSet(as.matrix(dat),phenoData = new("AnnotatedDataFrame", data = pdata), featureData = fd, expressionFamily=negbinomial.size(),lowerDetectionLimit=0.5)
mycds = estimateSizeFactors(mycds)
mycds = estimateDispersions(mycds)
mycds = detectGenes(mycds, min_expr = 0.1)
expressed_genes = row.names(subset(fData(mycds),num_cells_expressed >= 10))
mycds = mycds[expressed_genes,]
x = pData(mycds)$num_genes_expressed
x_1 = (x - mean(x)) / sd(x)
df <- data.frame(x = x_1)
ggplot(df, aes(x)) + geom_histogram(bins = 50) + geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
disp_table <- dispersionTable(mycds)
pData(mycds)$UMI <- Matrix::colSums(exprs(mycds))
ggplot(pData(mycds), aes(num_genes_expressed, UMI)) + geom_point()
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mycds <- setOrderingFilter(mycds, unsup_clustering_genes$gene_id)
mon_exprs = t(t(exprs(mycds)) /  pData(mycds)[, 'Size_Factor'])
dat_vst = vstExprs(mycds,expr_matrix=mon_exprs)
fwrite(data.frame(dat_vst,check.names=F),"BRCA_scRNAseq_all_patients_VST.dat",col.names=T,row.names=T,quote=F,sep="\t")

################################################################################################
## Compute ssGSEA scores
################################################################################################
library(GSVA)
imm = read.gmt.file("~/Research/pathways/immune_system.gmt")
ssgsea_out = gsva(as.matrix(dat_vst),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
fwrite(data.frame(ssgsea_out,check.names=F),"BRCA_scRNAseq_immuneGS_ssGSEA.dat",col.names=T,row.names=T,quote=F,sep="\t")

ach = read.gmt.file("~/Research/pathways/acharya_genesets.gmt")
ssgsea_ach = gsva(as.matrix(dat_vst),gset.idx.list=ach$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
fwrite(data.frame(ssgsea_ach,check.names=F),"BRCA_scRNAseq_acharyaGS_ssGSEA.dat",col.names=T,row.names=T,quote=F,sep="\t")
ssgsea_ach = as.matrix(Matrix::t(scale(Matrix::t(ssgsea_ach))))

################################################################################################
## Clustering Enrichment score data
################################################################################################
library(Rtsne); library(densityClust); library(plyr); library(dplyr); library(ggthemes); library(ggpubr)
max_components = 2; num_dim = 5; num_clusters = 15; perp = 30

FM = ssgsea_out
FM = as.matrix(Matrix::t(scale(Matrix::t(FM))))

## PC reduction
set.seed(123)
irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
irlba_pca_res <- irlba_res$x
topDim_pca <- irlba_pca_res

## tSNE
set.seed(123)
tsne_res <- Rtsne(as.matrix(topDim_pca), dims = max_components,pca = F,perplexity=perp)
tsne_data <- tsne_res$Y[, 1:max_components]
colnames(tsne_data) = c("Component_1", "Component_2")
rownames(tsne_data) = colnames(FM)

## Density clustering
set.seed(123)
dataDist <- dist(tsne_data)
dataClust <- densityClust(dataDist, gaussian = T)
delta_rho_df <- data.frame(delta = dataClust$delta, rho = dataClust$rho)
rho_threshold <- 0
delta_threshold <- sort(delta_rho_df$delta, decreasing = T)[num_clusters] - .Machine$double.eps
dataClust <- densityClust::findClusters(dataClust, rho = rho_threshold, delta = delta_threshold)
tsne_data = data.frame(tsne_data)
tsne_data$Cluster = as.factor(dataClust$cluster)
tcenters = tsne_data[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))
tsne_data = cbind(tsne_data,t(FM))
ave_out = aggregate(. ~ Cluster,tsne_data,"mean")[,-c(1:3)]
corr_out = cor(ave_out)
corr_out1 = cor(t(ave_out))
centers = tsne_data[,c(1:3)] %>% dplyr::group_by(Cluster) %>% summarize_all(funs(median))

## Identifying extreme clusters based on average TCA enrichment scores
ggplot(data.frame(tsne_data),aes(x=Component_1,y=Component_2,colour=Cluster,label=Cluster)) + geom_point(size=1) + geom_point(data = tcenters, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) +
      geom_text(data=tcenters,mapping = aes(label = Cluster), colour="black",size = 6) + theme_pander(boxes=T)
ggline(tsne_data,x="Cluster",y=rownames(FM)[4],add = "mean_se")

min_cyt = which.min(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data,"median")[,2])
max_cyt = which.max(aggregate(TCELL_CYTOLYTIC_ACT ~ Cluster,tsne_data,"median")[,2])
samples_up = rownames(tsne_data[tsne_data$Cluster %in% max_cyt,])
samples_down = rownames(tsne_data[tsne_data$Cluster %in% min_cyt,])
samples = c(samples_down,samples_up)
pheno_cyt = c(rep(0,length(samples_down)),rep(1,length(samples_up)))
pheno_cyt[pheno_cyt>0]<- "ECA_UP"
pheno_cyt[pheno_cyt!="ECA_UP"]<- "ECA_DOWN"
pheno_cyt = as.factor(pheno_cyt)
pheno_cyt = relevel(pheno_cyt,ref="ECA_DOWN")
table(pheno_cyt)
ssgsea_out2 = ssgsea_out[,match(samples,colnames(ssgsea_out))]
ssgsea_out2 = ssgsea_out2[c(1,2,3,6,7,10:16,19,20),]

library(gplots)
col_breaks = c( seq(-1,0,length=100),               # for red
                seq(0.01,0.8,length=100),           # for yellow
                seq(0.81,1,length=100))             # for green
my_palette <- rev(colorRampPalette(c("red", "white", "blue"))(n = 299))
heatmap.2(as.matrix(ssgsea_out2),main="",col=my_palette,notecol="black",breaks=col_breaks,scale="row",key=F,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.5,cexCol=0.8,labCol="",srtCol=90,dendrogram=c("none"),Rowv=NULL,Colv=NULL,trace="none",margins=c(9,15))

################################################################################################
# Identify Modules
# Use VST transformed data as the normalized data
################################################################################################
library(Seurat)
seed = 123

### Subsetting data based on Cyt Act
dat_new = dat_vst[,match(samples,colnames(dat_vst))]
var_out = apply(dat_new,1,var)
dat_new = dat_new[which(var_out!=0),]
dat_raw = dat[,match(samples,colnames(dat))]
dat_raw = dat_raw[match(rownames(dat_new),rownames(dat)),]

set.seed(seed)
seu_new = CreateSeuratObject(raw.data = dat_raw, min.cells=1,is.expr = 0.5, project = "BRCA", meta.data = pData(mycds)[match(colnames(dat_new),rownames(pData(mycds))),])
seu_new@meta.data$pheno = pheno_cyt
seu_new@data = dat_new
seu_new@scale.data = as.matrix(Matrix::t(scale(Matrix::t(dat_new))))
seu_new = FindVariableGenes(seu_new, do.plot = F, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.015, x.high.cutoff = 3.2, y.cutoff = 0.5)
hist(colSums(seu_new@data),breaks = 100,main = "Total expression after normalisation",xlab = "Sum of expression")
seu_new = RunPCA(seu_new, pcs.print = 0,pc.genes = seu_new@var.genes)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2,group.by="pheno")
PCElbowPlot(seu_new)
seu_new = FindClusters(object = seu_new, reduction.type = "pca", dims.use = 1:10, resolution = 1, print.output = 0, save.SNN = TRUE)
seu_new = JackStraw(object = seu_new, num.replicate = 1000)
JackStrawPlot(object = seu_new, PCs = 1:20)
PCHeatmap(object = seu_new, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
PCHeatmap(object = seu_new, pc.use = 7:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
set.seed(seed)
seu_new = RunTSNE(seu_new, dims.use = 1:10,do.fast=T)
TSNEPlot(seu_new, do.label = TRUE, pt.size = 2.5)
TSNEPlot(seu_new, do.label = FALSE, pt.size = 2.5,group.by="pheno")
FeaturePlot(seu_new, features.plot = c("PRF1","GZMA","HAVCR2","C10orf54"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("CD8A","CD8B","CD4","FOXP3"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FM1 = FM[,match(rownames(seu_new@meta.data),colnames(FM))]
seu_new@meta.data = cbind(seu_new@meta.data,t(FM1))
FM2 = ssgsea_ach[,match(rownames(seu_new@meta.data),colnames(ssgsea_ach))]
seu_new@meta.data = cbind(seu_new@meta.data,t(FM2))
seu_new@meta.data = seu_new@meta.data[,-c(30:43)]
FeaturePlot(seu_new, features.plot = c("MACROPHAGE_ACTIVITY","Treg","STAT1","CP_INHIBITORY","CP_STIMULATORY"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(seu_new, features.plot = c("EMT","WNT BETA CATENIN SIGNALING","TGF BETA SIGNALING","P53 PATHWAY","PI3K-AKT-MTOR SIGNALING"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(seu_new, features.plot = c("HYPOXIA","INVASIVENESS GENE SIGNATURE","ANGIOGENESIS","ECM","CHROMOSOMAL INSTABILITY","FIBROBLAST GROWTH FACTOR RECEPTOR"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)

dat_cyt = data.frame(as.matrix(seu_new@data),check.names=F)
dat_cyt = dat_cyt[match(seu_new@var.genes,rownames(dat_cyt)),]
write.table(dat_cyt,"BRCA_scRNAseq_all_patients_ECAdata.dat",sep="\t",col.names=NA,quote=F)

library(qvalue)
pc1_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,1])$qval<=0.05))
pc2_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,2])$qval<=0.05))
pc3_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,3])$qval<=0.05))
pc4_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,4])$qval<=0.05))
pc5_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,5])$qval<=0.05))
pc6_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,6])$qval<=0.05))
pc7_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,7])$qval<=0.05))
pc8_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,8])$qval<=0.05))
pc9_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,9])$qval<=0.05))
pc10_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,10])$qval<=0.05))

pc_gene_list = list(pc1_genes,pc2_genes,pc3_genes,pc4_genes,pc5_genes,pc6_genes,pc7_genes,pc8_genes,pc9_genes,pc10_genes)
for(i in 1:length(pc_gene_list)){
    write.table(pc_gene_list[[i]],paste("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/new_results/MOD",i,"_genes.txt",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
}

################################################################################################
## Module description
################################################################################################

library(qgraph)
names(pc_gene_list) = paste("MOD",1:length(pc_gene_list),sep="")
overlap = crossprod(table(stack(pc_gene_list)))
qgraph(cor(overlap),layout="spring",legend.cex = 0.4)
title("Gene-set overlap",line=2.5)

library(clusterProfiler); library(org.Hs.eg.db)
OUTsig_fin_list = list()
OUTsig_path = list()
for(i in 1:length(pc_gene_list)){
  print(i)
  egmt = suppressMessages(enricher(pc_gene_list[[i]],pvalueCutoff = 1,qvalueCutoff =1, TERM2GENE=read.gmt("/Users/ca31/Research/pathways/c2.cp.v5.2.symbols.gmt")))
  OUTsig = egmt@result
  OUTsig = OUTsig[,-c(1,2,8)]
  OUTsig$bonf = p.adjust(OUTsig$pvalue,"bonferroni")
  OUTsig_fin_list[[i]] = OUTsig[OUTsig$bonf<=0.01,]
  write.table(OUTsig_fin_list[[i]],paste("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/new_results/MOD",i,"_enrichment_analysis.txt",sep=""),sep="\t",col.names=NA,quote=F)
  OUTsig_path[[i]] = rownames(OUTsig_fin_list[[i]])
}
names(OUTsig_fin_list) = names(OUTsig_path) = paste("MOD",1:length(OUTsig_path),sep="")
overlap1 = crossprod(table(stack(OUTsig_path)))

# library(GENIE3)
# set.seed(123)
# exprMatr = dat_cyt[match(pc1_genes,row.names(dat_cyt)),]
# weightMat = GENIE3(exprMatr)

################################################################################################
## Final Plots
################################################################################################

pdf("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/new_results/metagene_plots_new.pdf",width=11,height=7)
ggplot(data.frame(tsne_data),aes(x=Component_1,y=Component_2,colour=Cluster,label=Cluster)) + geom_point(size=1) + geom_point(data = tcenters, mapping = aes(x = Component_1, y = Component_2), size = 0, alpha = 0) +
      geom_text(data=tcenters,mapping = aes(label = Cluster), colour="black",size = 6) + theme_pander(boxes=T)
ggline(tsne_data,x="Cluster",y=rownames(FM)[4],add = "mean_se")
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2,group.by="pheno")
JackStrawPlot(object = seu_new, PCs = 1:20)
PCHeatmap(object = seu_new, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
PCHeatmap(object = seu_new, pc.use = 7:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
TSNEPlot(seu_new, do.label = TRUE, pt.size = 1.5)
TSNEPlot(seu_new, do.label = F, pt.size = 1.5,group.by="pheno")
FeaturePlot(seu_new, features.plot = c("PRF1","GZMA","HAVCR2","C10orf54"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("CD8A","CD4","FOXP3","ITGAE"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("PDCD1","CTLA4","TIGIT","CD69"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("MACROPHAGE_ACTIVITY","Treg","STAT1","CP_INHIBITORY","CP_STIMULATORY"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(seu_new, features.plot = c("EMT","WNT BETA CATENIN SIGNALING","TGF BETA SIGNALING","P53 PATHWAY","PI3K-AKT-MTOR SIGNALING"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(seu_new, features.plot = c("HYPOXIA","INVASIVENESS GENE SIGNATURE","ANGIOGENESIS","ECM","CHROMOSOMAL INSTABILITY","FIBROBLAST GROWTH FACTOR RECEPTOR"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
heatmap.2(as.matrix(ssgsea_out2),main="",col=my_palette,notecol="black",breaks=col_breaks,scale="row",key=F,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.5,cexCol=0.8,labCol="",srtCol=90,dendrogram=c("none"),Rowv=NULL,Colv=NULL,trace="none",margins=c(9,15))
qgraph(cor(overlap),layout="spring",legend.cex = 0.4)
dev.off()

save.image("finding_metagenes_final.rds")
