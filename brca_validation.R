setwd("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/new_results")
load("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/new_results/finding_metagenes_final.rds")
load("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/new_results/validation_data.rds")

library(data.table); library(GSVA); library(plyr)

source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

################################################
## Random Forest Classifier based on metagenes
##
################################################

library(caret); library(survminer); library(survival); library(ggthemes); library(survutils); library(data.table)

rf_out = read.delim("rf_brca_downsampling.txt",header=T,row.names=1)
rf_sig = rf_out[rf_out$TCGA_BR_PFS<=0.05 & rf_out$TCGA_BR_OS<=0.05 & rf_out$METAB_BR_PFS<=0.05 & rf_out$METAB_BR_OS<=0.05,]

seed = 12345
pc_genes = pc5_genes
ntrees = 50
prop = 0.80
tr_ctrl = trainControl(method="cv",number=10,classProbs=TRUE,verboseIter=TRUE,sampling = "down", allowParallel=FALSE)

##################
# TCGA
##################

classifier_genes = Reduce("intersect",list(pc_genes, rownames(tcgaTN),rownames(metaB_TN)))
mtryGrid = expand.grid(mtry = if(length(classifier_genes) > 50){seq(50,min(length(classifier_genes),500),10)} else{seq(10,length(classifier_genes),10)})

### Training datasets
dat_cyt_train = dat_cyt[match(classifier_genes,rownames(dat_cyt)),]

### Test datasets
## TN
tcgaTN_cyt = tcgaTN[match(classifier_genes,rownames(tcgaTN)),]
metaB_TN_cyt = metaB_TN[match(classifier_genes,rownames(metaB_TN)),]
metaBV_TN_cyt = metaBV_TN[match(classifier_genes,rownames(metaBV_TN)),]
## ER
tcgaER_cyt = tcgaER[match(classifier_genes,rownames(tcgaER)),]
metaB_ER_cyt = metaB_ER[match(classifier_genes,rownames(metaB_ER)),]
metaBV_ER_cyt = metaBV_ER[match(classifier_genes,rownames(metaBV_ER)),]
## HER
tcgaHER_cyt = tcgaHER[match(classifier_genes,rownames(tcgaHER)),]
metaB_HER_cyt = metaB_HER[match(classifier_genes,rownames(metaB_HER)),]
metaBV_HER_cyt = metaBV_HER[match(classifier_genes,rownames(metaBV_HER)),]

set.seed(seed);
trainIndex = createDataPartition(pheno_cyt,prop,list=F,times=1)
trainData = t(dat_cyt_train[,trainIndex])
trainPheno = pheno_cyt[as.numeric(trainIndex)]

set.seed(seed);
train_out = train(trainData, trainPheno, method = "rf",importance=TRUE,metric="Accuracy",proximity=TRUE,trControl = tr_ctrl,tuneGrid = mtryGrid,ntree=ntrees)

## TCGA
## TN
clinTN$Prob = predict(train_out,apply(tcgaTN_cyt,1,scale),type="prob")[,1]
clinTN$Preds = predict(train_out,apply(tcgaTN_cyt,1,scale),type="raw")
tcga_brca_tn_pfs = as.numeric(pairwise_survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinTN,p.adjust="none")$p.value)
tcga_brca_tn_os = as.numeric(pairwise_survdiff(Surv(os_year,os_Ind) ~ Preds,clinTN,p.adjust="none")$p.value)
cox_TN_tcga =  suppressMessages(get_cox_res(clinTN,"pfs_year","pfs_Ind","Preds"))[,2:7]
clinTN = cbind(clinTN,apply(tcgaTN_gsva,1,scale))
## ER
clinER$Prob = predict(train_out,apply(tcgaER_cyt,1,scale),type="prob")[,1]
clinER$Preds = predict(train_out,apply(tcgaER_cyt,1,scale),type="raw")
tcga_brca_er_pfs = as.numeric(pairwise_survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinER,p.adjust="none")$p.value)
tcga_brca_er_os = as.numeric(pairwise_survdiff(Surv(os_year,os_Ind) ~ Preds,clinER,p.adjust="none")$p.value)
cox_ER_tcga =  suppressMessages(get_cox_res(clinER,"pfs_year","pfs_Ind","Preds"))[,2:7]
clinER = cbind(clinER,apply(tcgaER_gsva,1,scale))
## HER2
clinHER$Prob = predict(train_out,apply(tcgaHER_cyt,1,scale),type="prob")[,1]
clinHER$Preds = predict(train_out,apply(tcgaHER_cyt,1,scale),type="raw")
tcga_brca_her_pfs = as.numeric(pairwise_survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinHER,p.adjust="none")$p.value)
tcga_brca_her_os = as.numeric(pairwise_survdiff(Surv(os_year,os_Ind) ~ Preds,clinHER,p.adjust="none")$p.value)
cox_HER_tcga =  suppressMessages(get_cox_res(clinHER,"pfs_year","pfs_Ind","Preds"))[,2:7]
clinHER = cbind(clinHER,apply(tcgaHER_gsva,1,scale))

pdf("BRCA_tcga_surv_plots.pdf")
fitDSS_TCGA = survfit(Surv(pfs_year,pfs_Ind) ~ Preds,clinTN)
fitOS_TCGA = survfit(Surv(os_year,os_Ind) ~ Preds,clinTN)
ggsurvplot(fitDSS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinTN,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
ggsurvplot(fitOS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinTN,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_tcga_out = data.frame(matrix(0,nrow(tcgaTN_gsva)+1,8))
for(i in 1:nrow(cox_tcga_out)){
  print(i)
  cox_tcga_out[i,]=suppressMessages(get_cox_res(clinTN[,c(8,9,20:ncol(clinTN))],"pfs_year","pfs_Ind",colnames(clinTN[,-c(1:19)])[i]))
}
names(cox_tcga_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_tcga_out$term = as.factor(cox_tcga_out$term)
cox_tcga_out$term = factor(cox_tcga_out$term,levels=cox_tcga_out$term)
write.table(cox_tcga_out,"BRCA_TCGA_cox_results.txt",sep="\t",row.names=F,quote=F)
plot_cox_res(cox_tcga_out[c(1,5,9:11,18,20:21),]) + theme_classic()

fitDSS_TCGA = survfit(Surv(pfs_year,pfs_Ind) ~ Preds,clinER)
fitOS_TCGA = survfit(Surv(os_year,os_Ind) ~ Preds,clinER)
ggsurvplot(fitDSS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinER,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
ggsurvplot(fitOS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinER,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_tcga_out = data.frame(matrix(0,nrow(tcgaER_gsva)+1,8))
for(i in 1:nrow(cox_tcga_out)){
  print(i)
  cox_tcga_out[i,]=suppressMessages(get_cox_res(clinER[,c(8,9,20:ncol(clinER))],"pfs_year","pfs_Ind",colnames(clinER[,-c(1:19)])[i]))
}
names(cox_tcga_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_tcga_out$term = as.factor(cox_tcga_out$term)
cox_tcga_out$term = factor(cox_tcga_out$term,levels=cox_tcga_out$term)
plot_cox_res(cox_tcga_out[c(1,5,9:11,18,20:21),]) + theme_classic()

fitDSS_TCGA = survfit(Surv(pfs_year,pfs_Ind) ~ Preds,clinHER)
fitOS_TCGA = survfit(Surv(os_year,os_Ind) ~ Preds,clinHER)
ggsurvplot(fitDSS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinHER,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
ggsurvplot(fitOS_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=clinHER,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_tcga_out = data.frame(matrix(0,nrow(tcgaHER_gsva)+1,8))
for(i in 1:nrow(cox_tcga_out)){
  print(i)
  cox_tcga_out[i,]=suppressMessages(get_cox_res(clinHER[,c(8,9,20:ncol(clinTN))],"pfs_year","pfs_Ind",colnames(clinHER[,-c(1:19)])[i]))
}
names(cox_tcga_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_tcga_out$term = as.factor(cox_tcga_out$term)
cox_tcga_out$term = factor(cox_tcga_out$term,levels=cox_tcga_out$term)
plot_cox_res(cox_tcga_out[c(1,5,9:11,18,20:21),]) + theme_classic()
dev.off()

fit = survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinTN)
TN_p = 1 - pchisq(fit$chisq, length(fit$n) - 1)
TN_HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
TN_U95 = exp(log(TN_HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
TN_L95 = exp(log(TN_HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))

fit = survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinER)
ER_p = 1 - pchisq(fit$chisq, length(fit$n) - 1)
ER_HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
ER_U95 = exp(log(ER_HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
ER_L95 = exp(log(ER_HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))

fit = survdiff(Surv(pfs_year,pfs_Ind) ~ Preds,clinHER)
HER_p = 1 - pchisq(fit$chisq, length(fit$n) - 1)
HER_HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
HER_U95 = exp(log(HER_HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
HER_L95 = exp(log(HER_HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))


## Metabric
## TN
clinTN_metaB$Probs = predict(train_out,apply(metaB_TN_cyt,1,scale),type="prob")[,1]
clinTN_metaB$Preds = predict(train_out,apply(metaB_TN_cyt,1,scale))
clinTN_metaBV$Probs = predict(train_out,apply(metaBV_TN_cyt,1,scale),type="prob")[,1]
clinTN_metaBV$Preds = predict(train_out,apply(metaBV_TN_cyt,1,scale))
metaB_TN_valid = rbind(clinTN_metaB,clinTN_metaBV)
metab_TN_brca_pfs = as.numeric(pairwise_survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_TN_valid,p.adjust="none")$p.value)
metab_TN_brca_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_Ind) ~ Preds,metaB_TN_valid,p.adjust="none")$p.value)
cox_TN_metab = suppressMessages(get_cox_res(metaB_TN_valid,"DSS_days","DSS_event","Preds")[,2:7])
metaB_TN_valid = cbind(metaB_TN_valid,rbind(apply(metaB_TN_gsva,1,scale),apply(metaBV_TN_gsva,1,scale)))

## ER+
clinER_metaB$Probs = predict(train_out,apply(metaB_ER_cyt,1,scale),type="prob")[,1]
clinER_metaB$Preds = predict(train_out,apply(metaB_ER_cyt,1,scale))
clinER_metaBV$Probs = predict(train_out,apply(metaBV_ER_cyt,1,scale),type="prob")[,1]
clinER_metaBV$Preds = predict(train_out,apply(metaBV_ER_cyt,1,scale))
metaB_ER_valid = rbind(clinER_metaB,clinER_metaBV)
metab_ER_brca_pfs = as.numeric(pairwise_survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_ER_valid,p.adjust="none")$p.value)
metab_ER_brca_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_Ind) ~ Preds,metaB_ER_valid,p.adjust="none")$p.value)
cox_ER_metab = suppressMessages(get_cox_res(metaB_ER_valid,"DSS_days","DSS_event","Preds"))[,2:7]
metaB_ER_valid = cbind(metaB_ER_valid,rbind(apply(metaB_ER_gsva,1,scale),apply(metaBV_ER_gsva,1,scale)))

## HER2+
clinHER_metaB$Probs = predict(train_out,apply(metaB_HER_cyt,1,scale),type="prob")[,1]
clinHER_metaB$Preds = predict(train_out,apply(metaB_HER_cyt,1,scale))
clinHER_metaBV$Probs = predict(train_out,apply(metaBV_HER_cyt,1,scale),type="prob")[,1]
clinHER_metaBV$Preds = predict(train_out,apply(metaBV_HER_cyt,1,scale))
metaB_HER_valid = rbind(clinHER_metaB,clinHER_metaBV)
metab_HER_brca_pfs = as.numeric(pairwise_survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_HER_valid,p.adjust="none")$p.value)
metab_HER_brca_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_Ind) ~ Preds,metaB_HER_valid,p.adjust="none")$p.value)
cox_HER_metab = suppressMessages(get_cox_res(metaB_HER_valid,"DSS_days","DSS_event","Preds"))[,2:7]
metaB_HER_valid = cbind(metaB_HER_valid,rbind(apply(metaB_HER_gsva,1,scale),apply(metaBV_HER_gsva,1,scale)))

brca_rf_train = train_out

pdf("BRCA_metabric_surv_plots.pdf")
fitDSS_metaB = survfit(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_TN_valid)
ggsurvplot(fitDSS_metaB,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_TN_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
fitOS_metaB = survfit(Surv(OS_year,OS_Ind) ~ Preds,metaB_TN_valid)
ggsurvplot(fitOS_metaB,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_TN_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_metabric_out = data.frame(matrix(0,nrow(metaB_TN_gsva)+1,8))
for(i in 1:nrow(cox_metabric_out)){
  print(i)
  cox_metabric_out[i,]=suppressMessages(get_cox_res(metaB_TN_valid[,c(30:ncol(metaB_TN_valid))],"DSS_days","DSS_event",colnames(metaB_TN_valid[,-c(1:33)])[i]))
}
names(cox_metabric_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_metabric_out$term = as.factor(cox_metabric_out$term)
cox_metabric_out$term = factor(cox_metabric_out$term,levels=cox_metabric_out$term)
write.table(cox_metabric_out,"BRCA_metabric_cox_results.txt",sep="\t",row.names=F,quote=F)
plot_cox_res(cox_metabric_out[c(1,5,9:11,18,20:21),]) + theme_classic()

fitDSS_metaB = survfit(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_ER_valid)
ggsurvplot(fitDSS_metaB,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_ER_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
fitOS_metaB = survfit(Surv(OS_year,OS_Ind) ~ Preds,metaB_ER_valid)
ggsurvplot(fitOS_metaB,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_ER_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_metabric_out = data.frame(matrix(0,nrow(metaB_ER_gsva)+1,8))
for(i in 1:nrow(cox_metabric_out)){
  print(i)
  cox_metabric_out[i,]=suppressMessages(get_cox_res(metaB_ER_valid[,c(30:ncol(metaB_ER_valid))],"DSS_days","DSS_event",colnames(metaB_ER_valid[,-c(1:33)])[i]))
}
names(cox_metabric_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_metabric_out$term = as.factor(cox_metabric_out$term)
cox_metabric_out$term = factor(cox_metabric_out$term,levels=cox_metabric_out$term)
plot_cox_res(cox_metabric_out[c(1,5,9:11,18,20:21),]) + theme_classic()

fitDSS_metaB = survfit(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_ER_valid)
ggsurvplot(fitDSS_metaB,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_ER_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
fitOS_metaB = survfit(Surv(OS_year,OS_Ind) ~ Preds,metaB_ER_valid)
ggsurvplot(fitOS_metaB,palette = genespring.colors(2),risk.table=T,pval=T,data=metaB_ER_valid,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_metabric_out = data.frame(matrix(0,nrow(metaB_ER_gsva)+1,8))
for(i in 1:nrow(cox_metabric_out)){
  print(i)
  cox_metabric_out[i,]=suppressMessages(get_cox_res(metaB_ER_valid[,c(30:ncol(metaB_ER_valid))],"DSS_days","DSS_event",colnames(metaB_ER_valid[,-c(1:33)])[i]))
}
names(cox_metabric_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_metabric_out$term = as.factor(cox_metabric_out$term)
cox_metabric_out$term = factor(cox_metabric_out$term,levels=cox_metabric_out$term)
plot_cox_res(cox_metabric_out[c(1,5,9:11,18,20:21),]) + theme_classic()

dev.off()

fit = survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_TN_valid)
metab_TN_p = 1 - pchisq(fit$chisq, length(fit$n) - 1)
metab_TN_HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
metab_TN_U95 = exp(log(metab_TN_HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
metab_TN_L95 = exp(log(metab_TN_HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))

fit = survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_ER_valid)
metab_ER_p = 1 - pchisq(fit$chisq, length(fit$n) - 1)
metab_ER_HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
metab_ER_U95 = exp(log(metab_ER_HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
metab_ER_L95 = exp(log(metab_ER_HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))

fit = survdiff(Surv(DSS_days/365,DSS_event) ~ Preds,metaB_HER_valid)
metab_HER_p = 1 - pchisq(fit$chisq, length(fit$n) - 1)
metab_HER_HR = (fit$obs[2]/fit$exp[2])/(fit$obs[1]/fit$exp[1])
metab_HER_U95 = exp(log(metab_HER_HR) + qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))
metab_HER_L95 = exp(log(metab_HER_HR) - qnorm(0.975)*sqrt(1/fit$exp[2]+1/fit$exp[1]))

# Predictor
brca_predictors = data.frame(varImp(brca_rf_train,scale=T)$importance)
brca_predictors$Genes = rownames(brca_predictors)
brca_predictors = brca_predictors[order(brca_predictors$ECA_DOWN,decreasing=T),]
ggbarplot(brca_predictors[brca_rf_train$finalModel$mtry:1,],x="Genes",y="ECA_DOWN",fill="grey",ylab="Scaled Variable Importance",lab.size=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(1)))+theme(text = element_text(size=10))+coord_flip()
ggsave("BRCA_top_predictors.pdf",width=11.1,height=6.64,units="in")

datTCGA2 = datTCGA[match(brca_predictors$Genes[1:brca_rf_train$finalModel$mtry],rownames(datTCGA)),]
datTCGA2 = as.matrix(Matrix::t(scale(Matrix::t(datTCGA2))))
datTCGA2 = melt(datTCGA2)
datTCGA2$Subtype = clinical[match(datTCGA2$Var2,clinical$bcr_patient_barcode),"NEW_SUBTYPE"]
ggbarplot(datTCGA2,x="Subtype",y="value",add="mean_se",fill="grey") + stat_compare_means(method = "anova",label.y = 0.2)
ggsave("BRCA_tcga_top_predictors_enrichment.pdf")

clinical = clinical[match(colnames(datTCGA),clinical$bcr_patient_barcode),]
dim(clinical)

### Metabric
###
metaB = metaB[match(illumina_probes,rownames(metaB)),]
metaB$Genes = as.character(annot$ILMN_Gene)
metaB = setDT(metaB)[, lapply(.SD, mean), by = Genes]
metaB = as.data.frame(metaB)
rownames(metaB) = metaB[,1]
metaB = metaB[,-1]
metaB2 = metaB[match(brca_predictors$Genes[1:brca_rf_train$finalModel$mtry],rownames(metaB)),]
metaB2 = as.matrix(Matrix::t(scale(Matrix::t(metaB2))))

metaBV = metaBV[match(illumina_probes,rownames(metaBV)),]
metaBV$Genes = as.character(annot$ILMN_Gene)
metaBV = setDT(metaBV)[, lapply(.SD, mean), by = Genes]
metaBV = as.data.frame(metaBV)
rownames(metaBV) = metaBV[,1]
metaBV = metaBV[,-1]
metaBV2 = metaBV[match(brca_predictors$Genes[1:brca_rf_train$finalModel$mtry],rownames(metaBV)),]
metaBV2 = as.matrix(Matrix::t(scale(Matrix::t(metaBV2))))

metaB3 = cbind(metaB2,metaBV2)
metaB3 = melt(metaB3)
metaB_clin = rbind(metaB_clin,metaBV_clin)
metaB3$Subtype = metaB_clin[match(metaB3$Var2,metaB_clin$PATIENT_ID),"subtype"]
ggbarplot(metaB3,x="Subtype",y="value",add="mean_se",fill="grey") + stat_compare_means(method = "anova",label.y = 0.2)
ggsave("BRCA_METABRIC_top_predictors_enrichment.pdf")

# library(GSVA)
# brca_enrich = gsva(as.matrix(datTCGA),gset.idx.list=list(brca_predictors$Genes[1:train_out$finalModel$mtry]),method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
# datPRED = data.frame("ENRICHMENT"=as.numeric(brca_enrich),"Subtype"=clinical$NEW_SUBTYPE)
# rownames(datPRED)=colnames(datTCGA)
# datPRED$ENRICHMENT = scale(datPRED$ENRICHMENT)
# # ggboxplot(datPRED,x="Subtype",y="ENRICHMENT") + stat_compare_means(method = "anova",label.y = 3)
# ggbarplot(datPRED,x="Subtype",y="ENRICHMENT",add="mean_se",fill="grey") + stat_compare_means(method = "anova",label.y = 0.75)
# ggsave("BRCA_enrichment_score.pdf")

# tcga_prob = predict(train_out,apply(datTCGA[match(classifier_genes,rownames(datTCGA)),],1,scale),type="prob")[,2]
# datPRED = data.frame("ENRICHMENT"=as.numeric(brca_enrich),"PROB" = tcga_prob,"Subtype"=clinical$NEW_SUBTYPE)

library(limma); library(gplots); library(calibrate)

# 1. Training set
#
brca_pred_genes = brca_predictors$Genes[1:brca_rf_train$finalModel$mtry]
brca_train = brca_rf_train$trainingData
mod = model.matrix(~brca_train$.outcome)
colnames(mod)[2] = c("ECA_UP")
fit = eBayes(lmFit(t(brca_train[,-ncol(brca_train)]),mod))
top_gene_tab = topTable(fit,number=ncol(brca_train)-1,coef="ECA_UP")
top_gene_tab$Genes = rownames(top_gene_tab)

top_gene_tab1 = merge(brca_predictors,top_gene_tab,by.x=names(brca_predictors)[3],by.y=names(top_gene_tab)[7])
top_gene_tab1 = top_gene_tab1[,-3]
colnames(top_gene_tab1)[2]="Var_Importance"
top_gene_tab1 = top_gene_tab1[order(top_gene_tab1$Var_Importance,decreasing=T),]
top_gene_tab2 = subset(top_gene_tab1,top_gene_tab1$Genes %in% brca_pred_genes)
top_gene_tab2 = top_gene_tab2[order(top_gene_tab2$logFC,decreasing=T),]
top_gene_tab2$Genes = factor(top_gene_tab2$Genes,levels=top_gene_tab2$Genes)
top_gene_tab2$DIR = ifelse(top_gene_tab2$logFC<0,"Down Reg","Up Reg")

# top_sig = top_gene_tab2[top_gene_tab2$adj.P.Val<=0.05,]
# top_sig = top_sig[order(top_sig$logFC,decreasing=T),]
# top_up = head(top_sig,25)
# top_sig = top_sig[order(top_sig$logFC,decreasing=F),]
# top_down = head(top_sig,25)
# top_down=top_down[order(top_down$logFC,decreasing=T),]
# top_genes = rbind(top_up,top_down)
# top_genes$Genes = factor(top_genes$Genes,levels=top_genes$Genes)
# top_genes$DIR = ifelse(top_genes$logFC<0,"Down Reg","Up Reg")

ggplot(top_gene_tab2,aes(x=Genes,y=logFC,label=logFC))+geom_bar(stat='identity', aes(fill=DIR),width=.5)+coord_flip() + theme_classic()
ggsave("brca_train_top50_genes.pdf")

brca_train1 = brca_train[,-ncol(brca_train)]
top_gene_tab2 = top_gene_tab2[order(top_gene_tab2$logFC),]
brca_train1 = brca_train1[,match(top_gene_tab2$Genes,colnames(brca_train1))]
outcome = as.factor(as.character(brca_train$.outcome))

library(gplots)
col_breaks = c( seq(-1,0,length=100),               # for red
                seq(0.01,0.8,length=100),           # for yellow
                seq(0.81,1,length=100))             # for green
my_palette <- rev(colorRampPalette(c("red", "white", "blue"))(n = 299))

pdf("BRCA_gene_signature_heatmap.pdf",width=11.1,height=6.64)
# heatmap.2(t(brca_train1),scale="none",key=T,density.info="none",col=matlab.colors(100),dendrogram="none",margins=c(12,12),trace="none",cexRow=0.5,cexCol=0.75,symm=F,symkey=F,symbreaks=T,ColSideColors=genespring.colors(2)[as.numeric(outcome)],Colv=F,Rowv=F)
heatmap.2(t(brca_train1),main="",col=my_palette,notecol="black",breaks=col_breaks,scale="none",key=F,keysize=1.2,symkey=FALSE,density.info="none",cexRow=0.5,cexCol=0.8,labCol="",srtCol=90,dendrogram=c("none"),Rowv=NULL,Colv=NULL,trace="none",margins=c(9,15),ColSideColors=genespring.colors(2)[as.numeric(outcome)])
dev.off()

pdf("BRCA_DEgenes_volcanoplot.pdf")
with(top_gene_tab, plot(logFC, -log10(P.Value), pch=20, main="", xlim=c(-5,5),xlab="log fold change",ylab="-log10 p value",col="grey"))
# with(subset(top_gene_tab, adj.P.Val<0.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
# with(subset(top_gene_tab, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(top_gene_tab, adj.P.Val<0.05), points(logFC, -log10(P.Value), pch=20, col="black"))
with(subset(top_gene_tab, abs(logFC)>0.5), textxy(logFC, -log10(P.Value), labs=Genes, cex=.5))
dev.off()

#################
## GSE87517
################

gse26304 = read.delim("~/Research/DCIS/GSE26304/GSE26304_gene_level_new.txt",header=T,row.names=1,check.names=F)
gse26304_test = gse26304[match(brca_pred_genes,rownames(gse26304)),]
gse26304_test = as.matrix(Matrix::t(scale(Matrix::t(gse26304_test))))
gse26304_test = melt(gse26304_test)
gse26304_clin = read.delim("~/Research/DCIS/GSE26304/pheno_new.txt",header=T,row.names=1,check.names=F)

gse26304_test$Diagnosis = gse26304_clin[match(gse26304_test$Var2,rownames(gse26304_clin)),"Diagnosis"]
gse26304_test$Diagnosis = factor(gse26304_test$Diagnosis,levels = c("NORMAL","DCIS","MIXED","IDC"))
gse26304_test$PAM50 = gse26304_clin[match(gse26304_test$Var2,rownames(gse26304_clin)),"PAM50"]
gse26304_test$Subtype = gse26304_clin[match(gse26304_test$Var2,rownames(gse26304_clin)),"SUBTYPE"]
ggbarplot(gse26304_test,x="Diagnosis",y="value",add="mean_se",fill="grey") + stat_compare_means(method = "anova",label.y = 0.2)

gse26304_test2 = gse26304_test
gse26304_test2$Diagnosis = gsub("MIXED","IDC",gse26304_test2$Diagnosis)
gse26304_test2$Diagnosis = factor(gse26304_test2$Diagnosis,levels = c("NORMAL","DCIS","IDC"))
ggbarplot(gse26304_test2,x="Diagnosis",y="value",add="mean_se",fill="grey") + stat_compare_means(method = "anova",label.y = 0.2)
ggsave("gse26304_top_predictors_enrichment_overall.pdf")

#################
## GSE59246
################

gse59246 = read.delim("~/Research/DCIS/GSE59246/expression_dat/GSE59246_gene_level.txt",header=T,row.names=1,check.names=F)
geo_samp = read.delim("~/Research/DCIS/GSE59246/GEO_sample_info.txt",header=T,row.names=1)
colnames(gse59246) = geo_samp$SampleID
gse59246_clin = read.delim("~/Research/DCIS/GSE59246/mmc2.txt",header=T,row.names=1)
gse59246_clin = gse59246_clin[match(colnames(gse59246),rownames(gse59246_clin)),]
gse59246_clin = gse59246_clin[-grep("NA",rownames(gse59246_clin)),]
gse59246 = gse59246[,match(rownames(gse59246_clin),colnames(gse59246))]
gse59246_clin$Tissue_Type = as.factor(as.character(gse59246_clin$Tissue_Type))

gse59246_test = gse59246[match(brca_pred_genes,rownames(gse59246)),]
gse59246_test = as.matrix(Matrix::t(scale(Matrix::t(gse59246_test))))
gse59246_test = melt(gse59246_test)
gse59246_test$Tissue_Type = gse59246_clin[match(gse59246_test$Var2,rownames(gse59246_clin)),"Tissue_Type"]
gse59246_test$PAM50 = gse59246_clin[match(gse59246_test$Var2,rownames(gse59246_clin)),"PAM50"]

ggbarplot(gse59246_test,x="Tissue_Type",y="value",add="mean_se",fill="grey") + stat_compare_means(method = "anova",label.y = 0.2)
ggsave("gse59246_top_predictors_enrichment_overall.pdf")


save.image("brca_validation.rds")
