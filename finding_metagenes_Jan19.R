## Jan 19, 2019
## Chaitanya Acharya
## Finding metagenes

load("brca_preprocessing.rds")
devtools::install_github("JEFworks/MUDAN")

library(Seurat)
seu = CreateSeuratObject(raw.data = dat, min.cells = 3, min.genes = 200, is.expr = 0.5, project = "10X_BRCA",meta.data=pdata)
VlnPlot(object=seu,features.plot=c("nGene", "nUMI"),nCol=2)
GenePlot(object = seu, gene1 = "nUMI", gene2 = "nGene")
seu = FilterCells(object = seu, subset.names = c("nGene", "nUMI"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 7500))
seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

## GSVA
library(GSVA)
imm = read.gmt.file("~/Research/pathways/immune_system.gmt")
seu_ssgsea = gsva(as.matrix(seu@data),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
FM = seu_ssgsea
FM = as.matrix(Matrix::t(scale(Matrix::t(FM))))

######
#
library(MUDAN); library(irlba); library(Rtsne);library(plyr); library(dplyr); library(ggpubr); library(data.table)
max_components = 2; num_dim = 5; num_clusters = 15; perp = 30
set.seed(123)
irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
irlba_pca_res <- irlba_res$x
topDim_pca <- irlba_pca_res
rownames(topDim_pca)=colnames(FM)

# set.seed(123)
# pcs = getPcs(FM, nGenes=nrow(FM), nPcs=5, verbose=T)
# set.seed(123)
# emb = Rtsne::Rtsne(pcs, is_distance=FALSE, perplexity=30, num_threads=parallel::detectCores(),verbose=FALSE)$Y
# rownames(emb) <- rownames(pcs)
# emb = data.frame(emb)
# colnames(emb)=c("Component_1","Component_2")

set.seed(123)
tsne_res <- Rtsne(as.matrix(topDim_pca), dims = max_components,pca = F,perplexity=perp)
tsne_data <- tsne_res$Y[, 1:max_components]
colnames(tsne_data) = c("Component_1", "Component_2")
rownames(tsne_data) = colnames(FM)
tsne_data = data.frame(tsne_data)

set.seed(123)
com_graph = MUDAN::getComMembership(topDim_pca, k=30, method=igraph::cluster_louvain,verbose=T)
plotEmbedding(emb=tsne_data, groups=com_graph, main='Graph-based Community Detection', xlab="Component 1", ylab="Component 2", mark.clusters=TRUE, alpha=0.8, mark.cluster.cex=2,verbose=T,show.legend=F)

emb = data.frame(tsne_data,"Clusters"=com_graph,t(FM))
ggline(emb,x="Clusters",y=rownames(FM)[4],add = "mean_se")
emb_melt = data.table::melt(emb)
emb_melt = subset(emb_melt,variable %in% rownames(FM)[4])
cyt_ave = data.frame(aggregate(value ~ Clusters,emb_melt,"mean"),"Number"=as.numeric(table(com_graph)))
min_cyt = which.min(cyt_ave$value)
max_cyt = which.max(cyt_ave$value)
samples_up = rownames(emb[emb$Clusters %in% max_cyt,])
samples_down = rownames(emb[emb$Clusters %in% min_cyt,])
samples = c(samples_down,samples_up)
pheno_cyt = c(rep(0,length(samples_down)),rep(1,length(samples_up)))
pheno_cyt[pheno_cyt>0]<- "ECA_UP"
pheno_cyt[pheno_cyt!="ECA_UP"]<- "ECA_DOWN"
pheno_cyt = as.factor(pheno_cyt)
pheno_cyt = relevel(pheno_cyt,ref="ECA_DOWN")
table(pheno_cyt)
ssgsea_out2 = FM[,match(samples,colnames(FM))]
ssgsea_out2 = ssgsea_out2[c(1,2,3,6,7,10:16,19,20),]
library(gplots)
pdf("enrichment_score_matrix.pdf",width=11,height=6.71)
heatmap.2(as.matrix(ssgsea_out2),main="",Rowv=NULL,Colv=NULL,col=genespring.colors(2),scale="row",key=TRUE,keysize=1.2,symkey=FALSE,density.info="none",cexRow=1.0,labCol="",dendrogram=c("none"),trace="none")
dev.off()

###

seu_new = SubsetData(seu,cells.use=samples,do.clean=T,subset.raw=T)
seu_new = NormalizeData(seu_new, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(object=seu_new,features.plot=c("nGene", "nUMI"),nCol=2)
GenePlot(object = seu_new, gene1 = "nUMI", gene2 = "nGene")

seu_new = FindVariableGenes(seu_new, do.plot = T, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 0.5)
seu_new@meta.data$patientid = as.factor(seu_new@meta.data$patientid)
seu_new@meta.data$pheno = pheno_cyt
seu_new = ScaleData(object = seu_new)
seu_new = RunPCA(seu_new, pcs.print = 0, pc.genes = seu_new@var.genes)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2,group.by="pheno")
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2,group.by="patientid")
PCElbowPlot(seu_new)
PCHeatmap(object = seu_new, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
PCHeatmap(object = seu_new, pc.use = 7:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
seu_new = RunTSNE(seu_new, dims.use = 1:10,do.fast=T)
TSNEPlot(seu_new, do.label = TRUE, pt.size = 2.5)
TSNEPlot(seu_new, do.label = TRUE, pt.size = 2.5,group.by="pheno")

FeaturePlot(seu_new, features.plot = c("PRF1","GZMA","HAVCR2","C10orf54"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("CD8A","CD8B","CD4","FOXP3"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("PDCD1","CD274","TIGIT","ITGAE","CTLA4"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)

## Jackstraw analysis of the modules
seu_new = JackStraw(object = seu_new, num.replicate = 1000)
JackStrawPlot(object = seu_new, PCs = 1:20)

library(qvalue)
pc1_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,1])$qval<=0.05))
pc2_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,2])$qval<=0.05))
pc3_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,3])$qval<=0.05))
pc4_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,4])$qval<=0.05))
pc5_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,5])$qval<=0.05))
pc6_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,6])$qval<=0.05))
# pc7_genes = names(which(qvalue(seu_new@dr$pca@jackstraw@emperical.p.value[,7])$qval<=0.05))
pc_gene_list = list(pc1_genes,pc2_genes,pc3_genes,pc4_genes,pc5_genes,pc6_genes)

dat_cyt = data.frame(as.matrix(seu_new@data),check.names=F)
dat_cyt = dat_cyt[match(seu_new@var.genes,rownames(dat_cyt)),]

### PLOTS
pdf("metagene_plots.pdf",width=11,height=6.71)
plotEmbedding(emb=tsne_data, groups=com_graph, main='Graph-based Community Detection', xlab="Component 1", ylab="Component 2", mark.clusters=TRUE, alpha=0.8, mark.cluster.cex=2,verbose=T,show.legend=F)
ggline(emb,x="Clusters",y=rownames(FM)[4],add = "mean_se")
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2,group.by="pheno")
PCAPlot(object = seu_new, dim.1 = 1, dim.2 = 2,group.by="patientid")
JackStrawPlot(object = seu_new, PCs = 1:9)
PCHeatmap(object = seu_new, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE,cexRow=0.5)
# PCHeatmap(object = seu_new, pc.use = 7:8, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,use.full = FALSE)
TSNEPlot(seu_new, do.label = TRUE, pt.size = 1.5,group.by="pheno")
FeaturePlot(seu_new, features.plot = c("PRF1","GZMA","HAVCR2","C10orf54"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("CD8A","CD4","FOXP3","ITGAE"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(seu_new, features.plot = c("PDCD1","CTLA4","TIGIT","CD69"), nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
dev.off()

library(qgraph)
names(pc_gene_list) = paste("MOD",1:length(pc_gene_list),sep="")
overlap = crossprod(table(stack(pc_gene_list)))
qgraph(cor(table(stack(pc_gene_list))),legend.cex = 0.4)
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
  write.table(OUTsig_fin_list[[i]],paste("MOD",i,"_enrichment_analysis.txt",sep=""),sep="\t",col.names=NA,quote=F)
  OUTsig_path[[i]] = rownames(OUTsig_fin_list[[i]])
}
names(OUTsig_fin_list) = names(OUTsig_path) = paste("MOD",1:length(OUTsig_path),sep="")
overlap1 = crossprod(table(stack(OUTsig_path)))

library(GENIE3)
dat_mod1 = dat_cyt[match(pc1_genes,rownames(dat_cyt)),]
set.seed(123)
wt_mat1 = GENIE3(as.matrix(dat_mod1),nCores=4,verbose=T)

###########################################
save.image("finding_metagenes_Jan31.rds")
###########################################
