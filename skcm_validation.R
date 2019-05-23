##############################################
## Author: Chaitanya Acharya
## Code to validation melanoma datasets
## 1) TCGA primary melanoma
## 2) TCGA met. melanoma
## 3) GSE91061 (immune checkpoint dataset)
##############################################
load("rf_validations_downSampling.rds")
load("validation_data.rds")

# ggplotRegression <- function (fit){
#   require(ggplot2)
#   ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + geom_point() + stat_smooth(method = "lm", col = "red") +
#         labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
#         "Intercept =",signif(fit$coef[[1]],5 )," Slope =",signif(fit$coef[[2]], 5), " P =",signif(summary(fit)$coef[2,4], 5)))
# }

library(data.table); library(GSVA); library(caret); library(survminer); library(survival); library(ggthemes); library(survutils);
library(dplyr)

source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

imm = read.gmt.file("~/Research/pathways/immune_system.gmt")

tr_ctrl = trainControl(method="cv",number=10,classProbs=TRUE,verboseIter=TRUE,sampling = "down", allowParallel=FALSE)

#################
### SKCM primary
###
######################

tmp = subset(pc_fin_SKCM,PC=="5" & TCGA_SKCM_PFS<=0.05) %>% slice(which.min(TCGA_SKCM_PFS))
seed = 12345
ntrees = tmp$ntrees
pc_genes = pc5_genes
prop = tmp$train_prop/100

classifier_genes = Reduce("intersect",list(pc_genes,rownames(tcgaLUAD)))
### Training datasets
dat_cyt_train = dat_cyt[match(classifier_genes,rownames(dat_cyt)),]
### Test datasets
tcgaMELA_cyt = tcgaMELA[match(classifier_genes,rownames(tcgaMELA)),]

### RF classifier
mtryGrid = expand.grid(mtry = if(length(classifier_genes) > 50){seq(50,min(length(classifier_genes),500),10)} else{seq(10,length(classifier_genes),10)})
set.seed(seed)
trainIndex = createDataPartition(pheno_cyt,prop,list=F,times=1)
trainData = t(dat_cyt_train[,trainIndex])
trainPheno = pheno_cyt[as.numeric(trainIndex)]
set.seed(seed);
train_out = train(trainData, trainPheno, method = "rf",importance=TRUE,metric="Accuracy",proximity=TRUE,trControl = tr_ctrl,tuneGrid = mtryGrid,ntree=ntrees)
tcga_clinMELA$Probs = predict(train_out,apply(tcgaMELA_cyt,1,scale),type="prob")[,1]
tcga_clinMELA$Preds = predict(train_out,apply(tcgaMELA_cyt,1,scale),type="raw")
tcga_skcm_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinMELA,p.adjust="none")$p.value)
tcga_skcm_os = as.numeric(pairwise_survdiff(Surv(OS_year,OS_Ind) ~ Preds,tcga_clinMELA,p.adjust="none")$p.value)

tcga_clinMELA = cbind(tcga_clinMELA,apply(tcga_mela_gsva,1,scale))

pdf("SKCM_surv_plots.pdf")
fitDSS_TCGA = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinMELA)
ggsurvplot(fitDSS_TCGA,palette = genespring.colors(2),title="SKCM - MOD5",risk.table=T,pval=T,data=tcga_clinMELA,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
fitOS_TCGA = survfit(Surv(OS_year,OS_Ind) ~ Preds,tcga_clinMELA)
ggsurvplot(fitOS_TCGA,palette = genespring.colors(2),title="SKCM - MOD5",risk.table=T,pval=T,data=tcga_clinMELA,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_skcm_out = data.frame(matrix(0,nrow(tcga_mela_gsva)+1,8))
for(i in 1:nrow(cox_skcm_out)){
  print(i);
  cox_skcm_out[i,]=suppressMessages(get_cox_res(tcga_clinMELA[,c(2,3,17:ncol(tcga_clinMELA))],"PFS_year","PFS_Ind",colnames(tcga_clinMELA[,-c(1:16)])[i]))
}
names(cox_skcm_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_skcm_out$term = as.factor(cox_skcm_out$term)
cox_skcm_out$term = factor(cox_skcm_out$term,levels=cox_skcm_out$term)
write.table(cox_skcm_out,"SKCM_TCGA_cox_results.txt",sep="\t",row.names=F,quote=F)
plot_cox_res(cox_skcm_out[c(1,5,9:11,16:19),]) + theme_minimal()
skcm_predictors = data.frame(varImp(train_out,scale=T)$importance)
skcm_predictors$Genes = rownames(skcm_predictors)
skcm_predictors = skcm_predictors[order(skcm_predictors$ECA_DOWN,decreasing=T),]
ggbarplot(skcm_predictors[25:1,],x="Genes",y="ECA_DOWN",fill="grey",ylab="Scaled Variable Importance",lab.size=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(1)))+theme(text = element_text(size=10))+coord_flip()
dev.off()
skcm_rf_train = train_out

################
### SKCM mets
###
######################

tmp = subset(pc_fin_SKCMmet,PC=="5" & TCGA_SKCM_PFS<=0.05) %>% slice(which.min(TCGA_SKCM_PFS))
seed = 12345
ntrees = tmp$ntrees
pc_genes = pc5_genes
prop = tmp$train_prop/100

# classifier_genes = Reduce("intersect",list(pc_genes,rownames(tcgaMELA_met)))
classifier_genes = Reduce("intersect",list(pc_genes,rownames(tcgaLUAD)))
### Training datasets
dat_cyt_train = dat_cyt[match(classifier_genes,rownames(dat_cyt)),]
### Test datasets
tcgaMELA_met_cyt = tcgaMELA_met[match(classifier_genes,rownames(tcgaMELA_met)),]

mtryGrid = expand.grid(mtry = if(length(classifier_genes) > 50){seq(50,min(length(classifier_genes),500),10)} else{seq(10,length(classifier_genes),10)})
set.seed(seed)
trainIndex = createDataPartition(pheno_cyt,prop,list=F,times=1)
trainData = t(dat_cyt_train[,trainIndex])
trainPheno = pheno_cyt[as.numeric(trainIndex)]
set.seed(seed);
train_out = train(trainData, trainPheno, method = "rf",importance=TRUE,metric="Accuracy",proximity=TRUE,trControl = tr_ctrl,tuneGrid = mtryGrid,ntree=ntrees)
tcga_clinMELA_met$Probs = predict(train_out,apply(tcgaMELA_met_cyt,1,scale),type="prob")[,1]
tcga_clinMELA_met$Preds = predict(train_out,apply(tcgaMELA_met_cyt,1,scale),type="raw")
tcga_met_skcm_pfs = as.numeric(pairwise_survdiff(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinMELA_met,p.adjust="none")$p.value)
tcga_clinMELA_met = cbind(tcga_clinMELA_met,apply(tcga_mela_met_gsva,1,scale))
pdf("SKCM_met_surv_plots.pdf")
fitDSS_met_TCGA = survfit(Surv(PFS_year,PFS_Ind) ~ Preds,tcga_clinMELA_met)
ggsurvplot(fitDSS_met_TCGA,palette = genespring.colors(2),risk.table=T,pval=T,data=tcga_clinMELA_met,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
cox_skcm_met_out = data.frame(matrix(0,nrow(tcga_mela_met_gsva)+1,8))
for(i in 1:nrow(cox_skcm_met_out)){
  print(i);
  cox_skcm_met_out[i,]=suppressMessages(get_cox_res(tcga_clinMELA_met[,c(1:2,31:ncol(tcga_clinMELA_met))],"PFS_year","PFS_Ind",colnames(tcga_clinMELA_met[,c(32:ncol(tcga_clinMELA_met))])[i]))
}
names(cox_skcm_met_out)=c("term","estimate","std.error","statistic","p.value","conf.low","conf.high","test_type")
cox_skcm_met_out$term = as.factor(cox_skcm_met_out$term)
cox_skcm_met_out$term = factor(cox_skcm_met_out$term,levels=cox_skcm_met_out$term)
write.table(cox_skcm_met_out,"SKCM_TCGA_met_cox_results.txt",sep="\t",row.names=F,quote=F)
plot_cox_res(cox_skcm_met_out[c(1,5,9:11,18:22),]) + theme_classic()

skcm_met_predictors = data.frame(varImp(train_out,scale=T)$importance)
skcm_met_predictors$Genes = rownames(skcm_met_predictors)
skcm_met_predictors = skcm_met_predictors[order(skcm_met_predictors$ECA_DOWN,decreasing=T),]
ggbarplot(skcm_met_predictors[25:1,],x="Genes",y="ECA_DOWN",fill="grey",ylab="Scaled Variable Importance",lab.size=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(1)))+theme(text = element_text(size=10))+coord_flip()
dev.off()

###########
## GSE91061
## Immune checkpoint blockade therapy
##############################################

# load("GSE91061_validation.rds")

gse91061_pre = read.delim("~/Research/antiPD1_data/GSE91061/GSE91061_preTx_normalized.dat",header=T,row.names=1,check.names=F)
gse91061_CLIN_pre = read.delim("~/Research/antiPD1_data/GSE91061/GSE91061_preTx_clinical.dat",header=T,row.names=1,check.names=F)
gse91061_CLIN_pre$newResponse = ifelse(gse91061_CLIN_pre$myBOR=="SD" | gse91061_CLIN_pre$myBOR=="PD","No Response","Response")

library(vcd);library(ROCR); library(plotROC); library(pROC); library(limma); library(calibrate); library(gplots); library(ggthemes);

dat_cyt_new = dat_cyt[,ncol(dat_cyt):1]
pheno_cyt_new = rev(pheno_cyt)

dat_cyt1 = dat_cyt_new
pheno_cyt1 = pheno_cyt_new

classifier_genes = Reduce("intersect",list(pc_genes,rownames(gse91061_pre)))
mtryGrid = expand.grid(mtry = if(length(classifier_genes) > 50){seq(50,min(length(classifier_genes),500),10)} else{seq(10,length(classifier_genes),10)})

rf_out = read.delim("rf_pd_downSampling.txt",header=T,row.names=1)
rf_out = rf_out[rf_out$PFS<=0.05 & rf_out$FISHER<=0.05,]
ntrees = rf_out$ntrees
prop = rf_out$train_prop/100

dat_cyt_train = dat_cyt1[match(classifier_genes,rownames(dat_cyt1)),]
gse91061_cyt = gse91061_pre[match(classifier_genes,rownames(gse91061_pre)),]
set.seed(seed)
trainIndex = createDataPartition(pheno_cyt1,prop,list=F,times=1)
trainData = t(dat_cyt_train[,trainIndex])
trainPheno = pheno_cyt1[as.numeric(trainIndex)]
set.seed(seed)
train_out = train(trainData, trainPheno, method = "rf",importance=TRUE,metric="Accuracy",proximity=TRUE,trControl = tr_ctrl,tuneGrid = mtryGrid,ntree=ntrees)
gse91061_CLIN_pre$Prob = predict(train_out,apply(gse91061_cyt,1,scale),type="prob")[,1]
gse91061_CLIN_pre$Preds = predict(train_out,apply(gse91061_cyt,1,scale),type="raw")
gse91061_CLIN_pre_pfs = as.numeric(pairwise_survdiff(Surv(PFSWK,PFS_SOR) ~ Preds,gse91061_CLIN_pre,p.adjust="none")$p.value)
gse91061_CLIN_pre_os = as.numeric(pairwise_survdiff(Surv(OSWK,OS_SOR) ~ Preds,gse91061_CLIN_pre,p.adjust="none")$p.value)
gse91061_CLIN_pre_os = ifelse(length(gse91061_CLIN_pre_os)<1,NA,gse91061_CLIN_pre_os)
gse91061_CLIN_pre$pred_resp = as.character(gse91061_CLIN_pre$Preds)
gse91061_CLIN_pre$pred_resp = as.factor(ifelse(gse91061_CLIN_pre$pred_resp=="ECA_DOWN","No Response","Response"))
gse91061_CLIN_fisher1 = fisher.test(table(gse91061_CLIN_pre$Preds,gse91061_CLIN_pre$newResponse))$p.value
gse91061_CLIN_fisher2 = fisher.test(table(gse91061_CLIN_pre$Preds,gse91061_CLIN_pre$myBOR))$p.value
conf_mat = caret::confusionMatrix(table(gse91061_CLIN_pre$pred_resp,gse91061_CLIN_pre$newResponse),positive="Response")

pdf("antipd1_plots.pdf")

gse91061_rf_train = train_out
gse91061_predictors = data.frame(varImp(gse91061_rf_train,scale=T)$importance)
gse91061_predictors$Genes = rownames(gse91061_predictors)
gse91061_predictors = gse91061_predictors[order(gse91061_predictors$ECA_UP,decreasing=T),]
gse91061_predictors = gse91061_predictors[,-1]
colnames(gse91061_predictors)[1] = c("Overall")
ggbarplot(gse91061_predictors[1:60,],x="Genes",y="Overall",fill="grey",ylab="Scaled Variable Importance",lab.size=1)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(1)))+theme(text = element_text(size=8))+coord_flip()

# mosaic(structable(gse91061_CLIN_pre$Preds,gse91061_CLIN_pre$myBOR),shade=T,legend=T)
mosaic(structable(gse91061_CLIN_pre$Preds,gse91061_CLIN_pre$newResponse),shade=T,legend=T)
mosaic(structable(subset(gse91061_CLIN_pre,Cohort=="NIV3-PROG")$Preds,subset(gse91061_CLIN_pre,Cohort=="NIV3-PROG")$newResponse),shade=T,legend=T)
mosaic(structable(subset(gse91061_CLIN_pre,Cohort=="NIV3-NAIVE")$Preds,subset(gse91061_CLIN_pre,Cohort=="NIV3-NAIVE")$newResponse),shade=T,legend=T)
fitPFS = survfit(Surv(PFSWK,PFS_SOR) ~ Preds,gse91061_CLIN_pre)
ggsurvplot(fitPFS,palette = genespring.colors(2),risk.table=T,pval=T,data=gse91061_CLIN_pre,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
fitOS = survfit(Surv(OSWK,OS_SOR) ~ Preds,gse91061_CLIN_pre)
ggsurvplot(fitOS,palette = genespring.colors(2),risk.table=T,pval=T,data=gse91061_CLIN_pre,ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
# fitPFS_PROG = survfit(Surv(PFSWK,PFS_SOR) ~ Preds,subset(gse91061_CLIN_pre,Cohort=="NIV3-PROG"))
# ggsurvplot(fitPFS_PROG,palette = genespring.colors(2),risk.table=T,pval=T,data=subset(gse91061_CLIN_pre,Cohort=="NIV3-PROG"),ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")
# fitPFS_NAIVE = survfit(Surv(PFSWK,PFS_SOR) ~ Preds,subset(gse91061_CLIN_pre,Cohort=="NIV3-NAIVE"))
# ggsurvplot(fitPFS_NAIVE,palette = genespring.colors(2),risk.table=T,pval=T,data=subset(gse91061_CLIN_pre,Cohort=="NIV3-NAIVE"),ncensor.plot = FALSE,ncensor.plot.height = 0.25,risk.table.y.text.col = T,risk.table.fontsize = 3,risk.table.height=0.25,size = 0.5,ggtheme = theme_classic(),legend.title="")

fitER = survdiff(Surv(PFSWK,PFS_SOR) ~ Preds,gse91061_CLIN_pre)
ER_p = 1 - pchisq(fitER$chisq, length(fitER$n) - 1)
ER_HR = (fitER$obs[2]/fitER$exp[2])/(fitER$obs[1]/fitER$exp[1])
ER_U95 = exp(log(ER_HR) + qnorm(0.975)*sqrt(1/fitER$exp[2]+1/fitER$exp[1]))
ER_L95 = exp(log(ER_HR) - qnorm(0.975)*sqrt(1/fitER$exp[2]+1/fitER$exp[1]))

temp = read.delim("~/Research/antiPD1_data/GSE91061/archived_files/clinical_dat.txt",header=T,check.names=F)
temp = temp[match(gse91061_CLIN_pre$PatientID,temp$Patient),]
gse91061_CLIN_pre$Mutation_Load = temp$'Mutation Load'
gse91061_CLIN_pre$Neoantigen_Load = temp$'Neo-antigen Load'
gse91061_CLIN_pre$Neopeptide_Load = temp$'Neo-peptide Load'
gse91061_CLIN_pre$Cytolytic_score = temp$'Cytolytic Score'

sig_genes  = gse91061_predictors$Genes[1:gse91061_rf_train$finalModel$mtry]
gse91061_pre1 = gse91061_pre[match(sig_genes,rownames(gse91061_pre)),]
gse91061_pre1 = as.matrix(Matrix::t(scale(Matrix::t(gse91061_pre1))))
gse91061_pre1 = data.frame(t(gse91061_pre1),check.names=F)
gse91061_pre1$Response = gse91061_CLIN_pre$newResponse
gse91061_pre1$Cohort = gse91061_CLIN_pre$Cohort
gse91061_pre1_melt = melt(gse91061_pre1,id.vars=c("Response","Cohort"))
ggbarplot(gse91061_pre1_melt,x="Response",y="value",add="mean_se",fill="grey",title="Samples before Nivolumab") + stat_compare_means(method = "t.test",label.y = 0.1)
ggbarplot(gse91061_pre1_melt,x="Response",y="value",add="mean_se",fill="grey",title="Samples before Nivolumab",facet.by="Cohort") + stat_compare_means(method = "t.test",label.y = 0.1)


enrichment_scores = gsva(as.matrix(gse91061_pre),gset.idx.list=list(sig_genes),method="gsva",min.sz=1,max.sz=1000)
enrich_out = data.frame("ENRICH_SCORES"=as.numeric(scale(as.numeric(enrichment_scores))),"Cohort"=gse91061_CLIN_pre$Cohort,"PROB"=gse91061_CLIN_pre$Prob,
                        "PREDS"=gse91061_CLIN_pre$pred_resp,"Response"=gse91061_CLIN_pre$newResponse,"MutationLoad"=gse91061_CLIN_pre$Mutation_Load,"NeoAntigen_Load"=gse91061_CLIN_pre$Neoantigen_Load,"NeoPeptide_Load"=gse91061_CLIN_pre$Neopeptide_Load,"Cytolytic_Score"=gse91061_CLIN_pre$Cytolytic_score)
rownames(enrich_out) = colnames(gse91061_cyt)
enrich_out$PRED_CLASS = as.numeric(enrich_out$PREDS)-1
enrich_out$TRUE_CLASS = as.numeric(enrich_out$Response)-1
enrich_out = enrich_out[order(enrich_out$ENRICH_SCORES),]
enrich_out$Response = as.factor(as.character(enrich_out$Response))
ggbarplot(enrich_out,x="PREDS",y="ENRICH_SCORES",add="mean_se",fill="grey",title="Predicted Response") + stat_compare_means(method = "t.test",label.y = 0.75)
ggbarplot(enrich_out,x="PREDS",y="MutationLoad",add="mean_se",fill="grey",title="Predicted Response") + stat_compare_means(method = "t.test",label.y = 600)
ggbarplot(enrich_out,x="PREDS",y="NeoAntigen_Load",add="mean_se",fill="grey",title="Predicted Response") + stat_compare_means(method = "t.test",label.y = 600)
ggbarplot(enrich_out,x="PREDS",y="NeoPeptide_Load",add="mean_se",fill="grey",title="Predicted Response") + stat_compare_means(method = "t.test",label.y = 600)

ggbarplot(enrich_out,x="Response",y="ENRICH_SCORES",add="mean_se",fill="grey",title="True Response") + stat_compare_means(method = "t.test",label.y = 0.75)
ggbarplot(enrich_out,x="Response",y="PROB",add="mean_se",fill="grey",title="True Response") + stat_compare_means(method = "t.test",label.y = 0.75)
ggbarplot(enrich_out,x="Response",y="MutationLoad",add="mean_se",fill="grey",title="True Response") + stat_compare_means(method = "t.test",label.y = 600)
ggbarplot(enrich_out,x="Response",y="NeoAntigen_Load",add="mean_se",fill="grey",title="True Response") + stat_compare_means(method = "t.test",label.y = 600)
ggbarplot(enrich_out,x="Response",y="NeoPeptide_Load",add="mean_se",fill="grey",title="True Response") + stat_compare_means(method = "t.test",label.y = 600)

# preds = prediction(enrich_out$PRED_CLASS,enrich_out$TRUE_CLASS)
# auc = as.numeric(performance(preds, measure = "auc")@y.values)
# auc_prog = as.numeric(performance(prediction(subset(enrich_out,Cohort=="NIV3-PROG")$PRED_CLASS,subset(enrich_out,Cohort=="NIV3-PROG")$TRUE_CLASS), measure = "auc")@y.values)
# auc_naive = as.numeric(performance(prediction(subset(enrich_out,Cohort=="NIV3-NAIVE")$PRED_CLASS,subset(enrich_out,Cohort=="NIV3-NAIVE")$TRUE_CLASS), measure = "auc")@y.values)
#
# ggplot(enrich_out, aes(d = PRED_CLASS, m = TRUE_CLASS)) + geom_roc() + theme_classic() + annotate("text", x = .75, y = .25,label = paste("AUC = ", round(auc, 2)))+
#       scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .2)) + scale_y_continuous("Sensitivity", breaks = seq(0, 1, by = .2))
#
# ggplot(subset(enrich_out,Cohort=="NIV3-PROG"), aes(d = TRUE_CLASS, m = PROB)) + geom_roc() + theme_classic() + annotate("text", x = .75, y = .25,label = paste("AUC = ", round(auc_prog, 2)))+
#       scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .2)) + scale_y_continuous("Sensitivity", breaks = seq(0, 1, by = .2))
#
# ggplot(subset(enrich_out,Cohort=="NIV3-NAIVE"), aes(d = TRUE_CLASS, m = PROB)) + geom_roc() + theme_classic() + annotate("text", x = .75, y = .25,label = paste("AUC = ", round(auc_naive, 2)))+
#       scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .2)) + scale_y_continuous("Sensitivity", breaks = seq(0, 1, by = .2))

roc(enrich_out$TRUE_CLASS,enrich_out$PROB,ci=TRUE, boot.n=500, ci.alpha=0.95, stratified=FALSE,plot=TRUE, grid=F,print.auc=T, show.thres=T,smooth=F)
roc(subset(enrich_out,Cohort=="NIV3-PROG")$TRUE_CLASS,subset(enrich_out,Cohort=="NIV3-PROG")$PROB,ci=TRUE, boot.n=500, ci.alpha=0.95, stratified=FALSE,plot=TRUE, grid=F,print.auc=T, show.thres=T)
roc(subset(enrich_out,Cohort=="NIV3-NAIVE")$TRUE_CLASS,subset(enrich_out,Cohort=="NIV3-NAIVE")$PROB,ci=TRUE, boot.n=500, ci.alpha=0.95, stratified=FALSE,plot=TRUE, grid=F,print.auc=T, show.thres=T)

# 1. Training set
#
pd1_train = gse91061_rf_train$trainingData
mod = model.matrix(~pd1_train$.outcome)
colnames(mod)[2] = c("ECA_UP")
fit = eBayes(lmFit(t(pd1_train[,-ncol(pd1_train)]),mod))
top_gene_tab = topTable(fit,number=ncol(pd1_train)-1,coef="ECA_UP")
top_gene_tab$Genes = rownames(top_gene_tab)

top_gene_tab1 = merge(gse91061_predictors,top_gene_tab,by.x=names(gse91061_predictors)[2],by.y=names(top_gene_tab)[7])
colnames(top_gene_tab1)[2]="Var_Importance"
top_gene_tab1 = top_gene_tab1[order(top_gene_tab1$adj.P.Val),]
top_gene_tab2 = subset(top_gene_tab1,top_gene_tab1$Genes %in% sig_genes)
# top_up = top_gene_tab2[top_gene_tab2$logFC>0,]
# top_up = top_up[1:25,]
# top_up = top_up[order(top_up$logFC),]
# top_down = top_gene_tab2[top_gene_tab2$logFC<0,]
# top_down = top_down[order(top_down$logFC),]
# top_gene_tab3 = rbind(top_down,top_up)
top_gene_tab3 = top_gene_tab2[order(top_gene_tab2$logFC),]
top_gene_tab3$Genes = factor(top_gene_tab3$Genes,levels=top_gene_tab3$Genes)
top_gene_tab3$DIR = ifelse(top_gene_tab3$logFC<0,"Down Reg","Up Reg")
ggplot(top_gene_tab3,aes(x=Genes,y=logFC,label=logFC))+geom_bar(stat='identity', aes(fill=DIR),width=.5)+coord_flip() + theme_classic()

with(top_gene_tab2, plot(logFC, -log10(P.Value), pch=20, main="", xlim=c(-5,5),xlab="log fold change",ylab="-log10 p value",col="grey"))
with(subset(top_gene_tab2, adj.P.Val<0.05), points(logFC, -log10(P.Value), pch=20, col="black"))
with(subset(top_gene_tab2, abs(logFC)>0.75), textxy(logFC, -log10(P.Value), labs=Genes, cex=.5))

gse91061_valid = gse91061_cyt[match(sig_genes,rownames(gse91061_cyt)),]
gse91061_valid = data.frame(apply(gse91061_valid,1,scale),"Response"=as.factor(as.character(gse91061_CLIN_pre$newResponse)),"myBOR"=as.factor(as.character(gse91061_CLIN_pre$myBOR)),"Preds"=as.factor(as.character(gse91061_CLIN_pre$Preds)),"Cohort"=as.factor(as.character(gse91061_CLIN_pre$Cohort)),
                            "CytolyticScore"=as.numeric(gse91061_CLIN_pre$Cytolytic_score),"NeoantigenLoad"=as.numeric(gse91061_CLIN_pre$Neoantigen_Load),"NeopeptideLoad"=as.numeric(gse91061_CLIN_pre$Neopeptide_Load),"MutationLoad"=as.numeric(gse91061_CLIN_pre$Mutation_Load),check.names=F)
gse91061_valid = melt(gse91061_valid,id.vars=c("Response","myBOR","Preds","Cohort","CytolyticScore","NeoantigenLoad","NeopeptideLoad","MutationLoad"))
head(gse91061_valid)
ggbarplot(gse91061_valid,x="Response",y="value",add="mean_se",fill="grey",title="Full data") + stat_compare_means(method = "t.test",label.y = 0.3)
ggbarplot(subset(gse91061_valid,Cohort=="NIV3-PROG"),x="Response",y="value",add="mean_se",fill="grey",title="IPI EXPOSED COHORT") + stat_compare_means(method = "t.test",label.y = 0.6)
ggbarplot(subset(gse91061_valid,Cohort=="NIV3-NAIVE"),x="Response",y="value",add="mean_se",fill="grey",title="IPI NAIVE COHORT") + stat_compare_means(method = "t.test",label.y = 0.1)

### Breakdown by cohort
gse91061_CLIN_pre_naive = subset(gse91061_CLIN_pre,Cohort=="NIV3-NAIVE")
gse91061_pre_naive = gse91061_pre[,match(gse91061_CLIN_pre_naive$Sample,colnames(gse91061_pre))]
ctla4_pre_naive = as.numeric(scale(as.numeric(gse91061_pre_naive[grep("CTLA4",rownames(gse91061_pre_naive)),])))
pdcd1_pre_naive = as.numeric(scale(as.numeric(gse91061_pre_naive[grep("^PDCD1$",rownames(gse91061_pre_naive)),])))
gse91061_CLIN_pre_naive$ctla4 = ctla4_pre_naive
gse91061_CLIN_pre_naive$pd1 = pdcd1_pre_naive
ggbarplot(gse91061_CLIN_pre_naive,x="Preds",y="ctla4",add="mean_se",fill="grey",title="IPI NAIVE COHORT") + stat_compare_means(method = "t.test",label.y = 0.1)
ggbarplot(gse91061_CLIN_pre_naive,x="Preds",y="pd1",add="mean_se",fill="grey",title="IPI NAIVE COHORT") + stat_compare_means(method = "t.test",label.y = 0.8)

gse91061_CLIN_pre_prog = subset(gse91061_CLIN_pre,Cohort=="NIV3-PROG")
gse91061_pre_prog = gse91061_pre[,match(gse91061_CLIN_pre_prog$Sample,colnames(gse91061_pre))]
ctla4_pre_prog = as.numeric(scale(as.numeric(gse91061_pre_prog[grep("CTLA4",rownames(gse91061_pre_prog)),])))
pdcd1_pre_prog = as.numeric(scale(as.numeric(gse91061_pre_prog[grep("^PDCD1$",rownames(gse91061_pre_prog)),])))
gse91061_CLIN_pre_prog$ctla4 = ctla4_pre_prog
gse91061_CLIN_pre_prog$pd1 = pdcd1_pre_prog
ggbarplot(gse91061_CLIN_pre_prog,x="Preds",y="ctla4",add="mean_se",fill="grey",title="IPI EXPOSED COHORT") + stat_compare_means(method = "t.test",label.y = 0.1)
ggbarplot(gse91061_CLIN_pre_prog,x="Preds",y="pd1",add="mean_se",fill="grey",title="IPI EXPOSED COHORT") + stat_compare_means(method = "t.test",label.y = 0.1)

## Before and after Tx
library(qvalue)
gse91061_matched = read.delim("~/Research/antiPD1_data/GSE91061/GSE91061_matched_normalized.dat",header=T,row.names=1)
gse91061_matched_clin = read.delim("~/Research/antiPD1_data/GSE91061/GSE91061_matched_clinical.dat",header=T,row.names=1)
tmp = gse91061_matched_clin
tmp$Sample = gsub("_Pre","_On",tmp$Sample)
tmp$SampleType = "On"
gse91061_matched_clin = rbind(gse91061_matched_clin,tmp)
gse91061_matched_clin = gse91061_matched_clin[match(colnames(gse91061_matched),gse91061_matched_clin$Sample),]
gse91061_matched_sig = gse91061_matched[match(sig_genes,rownames(gse91061_matched)),]
pair = as.factor(as.numeric(gse91061_matched_clin$Patient))
on_pre = as.factor(as.numeric(gse91061_matched_clin$SampleType))
# out = apply(gse91061_matched,1,function(x) t.test(x ~ gse91061_matched_clin$SampleType,paired=T)$p.value)
# out_qval = qvalue(out)$qval
mod = model.matrix(~pair+on_pre)
top_diff_expr = topTable(eBayes(lmFit(gse91061_matched_sig, mod)), coef="on_pre2",number=nrow(gse91061_matched_sig))
top_diff_expr$qval = qvalue(top_diff_expr$P.Value)$qval
gse91061_matched_new = data.frame(apply(gse91061_matched_sig,1,scale),"Response"=gse91061_matched_clin$myBOR,"Treatment"=gse91061_matched_clin$SampleType,"Cohort"=gse91061_matched_clin$Cohort)
gse91061_matched_new$Response = ifelse(gse91061_matched_new$Response=="PD" | gse91061_matched_new$Response=="SD","No Response","Response")
gse91061_matched_melt = melt(gse91061_matched_new,id.vars=c("Response","Treatment","Cohort"))
ggbarplot(gse91061_matched_melt,x="Response",y="value",add="mean_se",fill="grey",title="Matched Samples") + stat_compare_means(method = "t.test",label.y = 0.1)
ggbarplot(gse91061_matched_melt,x="Response",y="value",add="mean_se",fill="grey",title="Matched Samples",facet.by="Treatment") + stat_compare_means(method = "t.test",label.y = 0.3)
ggbarplot(gse91061_matched_melt,x="Treatment",y="value",add="mean_se",fill="grey",facet.by="Cohort",title="Matched Samples") + stat_compare_means(method = "t.test",label.y = 0.25)
ggbarplot(subset(gse91061_matched_melt,Cohort=="NIV3-PROG"),x="Response",y="value",add="mean_se",fill="grey",facet.by="Treatment",title="Matched Samples:IPI EXPOSED COHORT") + stat_compare_means(method = "t.test",label.y = 0.5)
ggbarplot(subset(gse91061_matched_melt,Cohort=="NIV3-NAIVE"),x="Response",y="value",add="mean_se",fill="grey",facet.by="Treatment",title="Matched Samples:IPI NAIVE COHORT") + stat_compare_means(method = "t.test",label.y = 0.4)

## Before and after Tx
gse91061_genes_before_after = aggregate(value ~ variable + Treatment + Cohort,gse91061_matched_melt,FUN="mean")
ggplot(gse91061_genes_before_after, aes(x = variable, y = Treatment, fill = value)) + geom_tile(color = "white")+scale_fill_gradientn(limits = c(-0.25,0.25),colours=genespring.colors(3))+theme_classic()+facet_grid(~Cohort,drop=T,space="free",scales="free")+coord_equal()+theme(axis.text=element_text(size=5))+coord_flip()
ggplot(gse91061_genes_before_after, aes(x = variable, y = Treatment, fill = value)) + geom_tile(color = "white")+scale_fill_gradientn(limits = c(-0.25,0.25),colours=genespring.colors(3))+theme_classic()+coord_equal()+theme(axis.text=element_text(size=7))+coord_flip()
ggbarplot(gse91061_matched_melt,x="Treatment",y="value",add="mean_se",fill="grey",facet.by=c("Cohort","Response"),title="Matched Samples") + stat_compare_means(method = "t.test",label.y = 0.6)

gse91061_genes_before_after = aggregate(value ~ variable + Treatment,gse91061_matched_melt,FUN="mean")
gse91061_genes_before_after = gse91061_genes_before_after[order(gse91061_genes_before_after$value,decreasing=T),]
gse91061_genes_before_after$variable = factor(gse91061_genes_before_after$variable,levels=split(gse91061_genes_before_after,gse91061_genes_before_after$Treatment)$On$variable)
ggplot(gse91061_genes_before_after,aes(x=variable,y=value,fill=Treatment)) + geom_bar(stat="identity",width=0.6)+coord_flip()+theme_classic()+scale_fill_manual(values=genespring.colors(2))
ggplot(gse91061_genes_before_after, aes(x = variable, y = Treatment, fill = value)) + geom_tile(color = "white")+scale_fill_gradientn(limits = c(-0.25,0.25),colours=genespring.colors(3))+theme_classic()+coord_equal()+theme(axis.text=element_text(size=7))+coord_flip()

### Test BRCA mod genes
temp_clin = gse91061_CLIN_pre
mod = model.matrix(~newResponse,temp_clin)
colnames(mod)[2]=c("Response")
fit = eBayes(lmFit(gse91061_pre,mod))
top_gene_tab = topTable(fit,number=nrow(gse91061_pre),coef="Response")
barcodeplot(top_gene_tab$logFC,index = which(rownames(top_gene_tab) %in% sig_genes))
barcodeplot(top_gene_tab$logFC,index = which(rownames(top_gene_tab) %in% sig_genes),gene.weights=top_gene_tab$logFC[which(rownames(top_gene_tab) %in% sig_genes)],weights.label="logFC")

dev.off()

save.image("skcm_validation.rds")
