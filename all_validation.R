##############################################
## Author: Chaitanya R. Acharya
## Code to validate modules in multiple tumor datasets
#######################################################

load("rf_validations_downSampling.rds")

library(GSVA);
library(data.table);
library(caret);
library(ggplot2);
library(survminer);
library(survival);
library(survutils);
library(ggpubr);
library(ggthemes);

source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

imm = read.gmt.file("~/Research/pathways/immune_system.gmt")
pc_out_list = list(pc_fin_SKCM,pc_fin_SKCMmet,pc_fin_HNSC,pc_fin_COAD,pc_fin_KIRC,pc_fin_KIRP,pc_fin_CESC,pc_fin_BLCA,pc_fin_LIHC,pc_fin_GBM,pc_fin_LGG,pc_fin_UCEC)
tcga_dat_list = list(tcgaMELA,tcgaMELA_met,tcgaHNSC,tcgaCOAD,tcgaKIRC,tcgaKIRP,tcgaCESC,tcgaBLCA,tcgaLIHC,tcgaGBM,tcgaLGG,tcgaUCEC)
tcga_clin_list =  list(tcga_clinMELA,tcga_clinMELA_met,tcga_clinHNSC,tcga_clinCOAD,tcga_clinKIRC,tcga_clinKIRP,tcga_clinCESC,tcga_clinBLCA,tcga_clinLIHC,tcga_clinGBM,tcga_clinLGG,tcga_clinUCEC)
tcga_gsva_list = list(tcga_mela_gsva,tcga_mela_met_gsva,tcga_hnsc_gsva,tcga_COAD_gsva,tcga_KIRC_gsva,tcga_KIRP_gsva,tcga_CESC_gsva,tcga_BLCA_gsva,tcga_LIHC_gsva,tcga_GBM_gsva,tcga_LGG_gsva,tcga_UCEC_gsva)
title = c("SKCM","SKCM met","HNSC","COAD","KIRC","KIRP","CESC","BLCA","LIHC","GBM","LGG","UCEC")
tr_ctrl = trainControl(method="cv",number=10,classProbs=TRUE,verboseIter=F,sampling="down",allowParallel=F)
seed = 12345
pfs_out = pc_out = rep(0,length(pc_out_list))
predictors_out = list()
cox_out = list()
pdf("survival_plots_tumor_types.pdf")
for(i in 1:length(pc_out_list)){
  print(title[i])
  dat = pc_out_list[[i]]
  params = dat[dat[,1]<=0.05,]
  params = subset(params,PC %in% c(3,4,5))
  params = params[which.min(params[,1]),]
  ntrees = params$ntrees
  prop = params$train_prop/100
  PC = params$PC
  pc_genes = pc_gene_list[[PC]]
  pc_out[i] = PC
  classifier_genes = Reduce("intersect",list(pc_genes,rownames(tcgaLUAD)))

  dat_cyt_train = dat_cyt[match(classifier_genes,rownames(dat_cyt)),]
  test_data = tcga_dat_list[[i]][match(classifier_genes,rownames(tcga_dat_list[[i]])),]
  mtryGrid = expand.grid(mtry = if(length(classifier_genes) > 50){seq(50,min(length(classifier_genes),500),10)} else{seq(10,length(classifier_genes),10)})
  set.seed(seed)
  trainIndex = createDataPartition(pheno_cyt,prop,list=F,times=1)
  trainData = t(dat_cyt_train[,trainIndex])
  trainPheno = pheno_cyt[as.numeric(trainIndex)]
  set.seed(seed);
  train_out = train(trainData, trainPheno, method = "rf",importance=TRUE,metric="Accuracy",proximity=TRUE,trControl = tr_ctrl,tuneGrid = mtryGrid,ntree=ntrees)
  tcga_clin_list[[i]]$Probs = predict(train_out,apply(test_data,1,scale),type="prob")[,1]
  tcga_clin_list[[i]]$Preds = predict(train_out,apply(test_data,1,scale),type="raw")
  pfs_out[[i]] = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clin_list[[i]],p.adjust="none")$p.value)
  fitDSS_TCGA = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clin_list[[i]])
  print(ggsurvplot(fitDSS_TCGA,title=title[i],palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clin_list[[i]],ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title=""))
  cox = data.frame(tcga_clin_list[[i]],t(tcga_gsva_list[[i]]))
  terms = c("PFS_year","PFS_Ind","Preds","TCELL_CYTOLYTIC_ACT","Tfh_SIG","Th1_SIG","Treg","IMMUNE1","IMMUNE2","INTERFERON_SIG")
  cox = cox[,colnames(cox) %in% terms]
  coxx = data.frame(matrix(0,8,8))
  for(j in 3:ncol(cox)){
    coxx[(j-2),]=suppressMessages(get_cox_res(cox,"PFS_year","PFS_Ind",colnames(cox)[j]))
  }
  colnames(coxx) = c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
  coxx = t(coxx[,c(-1,-8)])
  colnames(coxx) = paste(title[i],terms[-c(1:2)],sep="_")
  cox_out[[i]] = coxx
  predictors =  data.frame(varImp(train_out,scale=T)$importance)
  predictors = predictors[order(predictors$ECA_UP,decreasing=T),]
  predictors$Genes = rownames(predictors)
  predictors = predictors[,-1]
  predictors = predictors$Genes[1:train_out$finalModel$mtry]
  predictors_out[[i]] = predictors
  write.table(predictors,paste(title[i],"_predictors.txt",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
}
dev.off()

predictors_out[[13]] = scan("BRCA_predictors.txt",what="")
predictors_out[[14]] = scan("antiPD1_predictors.txt",what="")
names(predictors_out)[13] = "BRCA"
names(predictors_out)[14] = "antiPD1"
pred_all = as.matrix(table(stack(predictors_out)))
pred_all = pred_all[do.call(order, c(decreasing=T,lapply(1:NCOL(pred_all), function(i) pred_all[, i]))), ]
pred_all_1 = pred_all[1:50,]; pred_all_2 = pred_all[51:100,]; pred_all_3 = pred_all[101:150,]; pred_all_4 = pred_all[151:214,]
pred_all_melt1 = melt(pred_all_1); pred_all_melt2 = melt(pred_all_2);pred_all_melt3 = melt(pred_all_3);pred_all_melt4 = melt(pred_all_4);
colnames(pred_all_melt1)[1:3] = colnames(pred_all_melt2)[1:3] = colnames(pred_all_melt3)[1:3] = colnames(pred_all_melt4)[1:3] = c("Genes","Tumor_Type","value")
pdf("gene_map.pdf",width=11,height=2)
ggplot(pred_all_melt1, aes(Genes,Tumor_Type)) + geom_tile(aes(fill = value), colour = "black") + scale_fill_gradient2(low ="blue", high ="black") + theme_tufte() + theme(axis.text.x=element_text(size=8,angle=90)) + theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())
ggplot(pred_all_melt2, aes(Genes,Tumor_Type)) + geom_tile(aes(fill = value), colour = "black") + scale_fill_gradient2(low ="blue", high ="black") + theme_tufte() + theme(axis.text.x=element_text(size=8,angle=90)) + theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())
ggplot(pred_all_melt3, aes(Genes,Tumor_Type)) + geom_tile(aes(fill = value), colour = "black") + scale_fill_gradient2(low ="blue", high ="black") + theme_tufte() + theme(axis.text.x=element_text(size=8,angle=90)) + theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())
ggplot(pred_all_melt4, aes(Genes,Tumor_Type)) + geom_tile(aes(fill = value), colour = "black") + scale_fill_gradient2(low ="blue", high ="black") + theme_tufte() + theme(axis.text.x=element_text(size=8,angle=90)) + theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())
dev.off()

cox_out1 = do.call("cbind",cox_out)
cox_out1 = round(cox_out1,3)
cox_out2 = rbind(apply(cox_out1,2,function(x) paste0(x[1]," (",x[3],"-",x[4],")")),cox_out1[2,])
rownames(cox_out2) = c("HR","pval")
write.table(cox_out2,"table1.txt",sep="\t",col.names=NA,quote=F)

save.image("all_validation.rds")
