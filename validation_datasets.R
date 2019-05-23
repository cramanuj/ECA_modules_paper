#!/usr/bin/env Rscript

### Preparing files for validation
### Author: Chaitanya R. Acharya

library(data.table); library(GSVA); library(plyr)

source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

imm = read.gmt.file("~/Research/pathways/immune_system.gmt")

#########
## BREAST
#########
print("preprocessing breast (TCGA)")
## TCGA data

clinical = fread("~/Research/BRCA/tcga_brca_clinical_data.txt",header=TRUE,data.table=F)
dim(clinical)
clinical = clinical[clinical$history_of_neoadjuvant_treatment!=1,]
clinical = clinical[clinical$history_other_malignancy!=1,]
clinical$subtype = as.character(clinical$NEW_SUBTYPE)
clinical$os_year = as.numeric(as.character(clinical$os_year))
clinical$os_Ind = as.numeric(as.character(clinical$os_Ind))
clinical$rfs_year = as.numeric(as.character(clinical$rfs_year))
clinical$rfs_Ind = as.numeric(as.character(clinical$rfs_Ind))
clinical$pfs_year = as.numeric(as.character(clinical$pfs_year))
clinical$pfs_Ind = as.numeric(as.character(clinical$pfs_Ind))
clinical$Age = as.numeric(as.character(clinical$Age))
clinical = clinical[,-c(18:19)]
clinical = clinical[!is.na(clinical$NEW_SUBTYPE),]
samples_clinical = as.character(clinical$bcr_patient_barcode)

datTCGA = fread("~/Research/BRCA/breast_immune/TCGA_BR_Exp.txt",header=T,data.table=F);
rownames(datTCGA) = datTCGA[,1];
datTCGA = datTCGA[,-1];
rownames(clinical)=clinical$bcr_patient_barcode

# TN
clinTN = clinical[clinical$subtype=="TN",]
clinTN = clinTN[clinTN$PR != 1,]
clinTN = clinTN[!is.na(clinTN$subtype),]
keep_samples = intersect(clinTN$bcr_patient_barcode,colnames(datTCGA))
tcgaTN = datTCGA[,match(keep_samples,colnames(datTCGA))]
clinTN = clinTN[match(keep_samples,clinTN$bcr_patient_barcode),]

# ER
clinER = clinical[clinical$subtype=="ER+",]
keep_ER = intersect(clinER$bcr_patient_barcode,colnames(datTCGA))
tcgaER = datTCGA[,match(keep_ER,colnames(datTCGA))]

# HER
clinHER = clinical[clinical$subtype=="HER+",]
keep_HER = intersect(clinHER$bcr_patient_barcode,colnames(datTCGA))
tcgaHER = datTCGA[,match(keep_HER,colnames(datTCGA))]

## GSVA/ssGSEA
tcgaTN_gsva = gsva(as.matrix(tcgaTN),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcgaER_gsva = gsva(as.matrix(tcgaER),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcgaHER_gsva = gsva(as.matrix(tcgaHER),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)


###############################
## METABRIC datasets
## 1) Discovery 2) Validation
###############################
print("preprocessing breast (Metabric)")

library(impute)

# 1) MB Discovery
metaB = fread("~/Research/Metabric/_ega-box-04_discovery_ExpressionMatrix.txt",header=F,data.table=F)
dim(metaB)
rownames(metaB) = metaB[,1]; metaB = metaB[,-1]
system("head -1 ~/Research/Metabric/_ega-box-04_discovery_ExpressionMatrix.txt > metaB_samples.txt")
metaB_samples = scan("metaB_samples.txt",what="")
colnames(metaB) = metaB_samples

temp = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_patient.txt",header=T)
temp = temp[match(colnames(metaB),temp$PATIENT_ID),]
temp1 = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_sample.txt",header=T)
temp1 = temp1[match(colnames(metaB),temp1$PATIENT_ID),]
temp1$ER_STATUS = as.character(temp1$ER_STATUS)
temp1$HER2_STATUS = as.character(temp1$HER2_STATUS)
temp1$PR_STATUS = as.character(temp1$PR_STATUS)
temp1$TN = laply(1:nrow(temp1),function(i) length(grep("-",temp1[i,6:8])))
subtype = ifelse(temp1$ER_STATUS=="+" & temp1$HER2_STATUS=="-","ER+",ifelse(temp1$ER_STATUS=="+" | temp1$ER_STATUS=="-" & temp1$HER2_STATUS=="+","HER2+","TN"))
temp1$subtype = subtype
temp2 = cbind(temp,temp1[,-c(1:2)])
metaB_clin = temp2

## TN
clinTN_metaB = metaB_clin[metaB_clin$TN=="3",]
clinTN_metaB$OS_year = clinTN_metaB$OS_MONTHS/12
clinTN_metaB$OS_Ind = ifelse(clinTN_metaB$OS_STATUS=="DECEASED",1,0)
metaB_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaB_dss1 = metaB_dss[match(clinTN_metaB$PATIENT_ID,metaB_dss$SampleID),]
clinTN_metaB$DSS_days = metaB_dss1$time
clinTN_metaB$DSS_event = metaB_dss1$status
metaB_TN = metaB[,match(clinTN_metaB$PATIENT_ID,colnames(metaB))]
annot = read.delim("~/Illumina/HumanHT-12_V4_0_R2_15002873_B.txt",header=T)
annot = annot[,c(14,5,12)]
annot$Symbol = as.character(annot$Symbol)
illumina_probes = intersect(rownames(metaB_TN),annot$Probe_Id)
annot = annot[match(illumina_probes,annot$Probe_Id),]
metaB_TN = metaB_TN[match(illumina_probes,rownames(metaB_TN)),]
metaB_TN$Genes = as.character(annot$ILMN_Gene)
metaB_TN = setDT(metaB_TN)[, lapply(.SD, mean), by = Genes]
metaB_TN = as.data.frame(metaB_TN)
rownames(metaB_TN) = metaB_TN[,1]
metaB_TN = metaB_TN[,-1]

## ER
clinER_metaB = metaB_clin[metaB_clin$subtype=="ER+",]
clinER_metaB$OS_year = clinER_metaB$OS_MONTHS/12
clinER_metaB$OS_Ind = ifelse(clinER_metaB$OS_STATUS=="DECEASED",1,0)
metaB_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaB_dss1 = metaB_dss[match(clinER_metaB$PATIENT_ID,metaB_dss$SampleID),]
clinER_metaB$DSS_days = metaB_dss1$time
clinER_metaB$DSS_event = metaB_dss1$status
metaB_ER = metaB[,match(clinER_metaB$PATIENT_ID,colnames(metaB))]
annot = read.delim("~/Illumina/HumanHT-12_V4_0_R2_15002873_B.txt",header=T)
annot = annot[,c(14,5,12)]
annot$Symbol = as.character(annot$Symbol)
illumina_probes = intersect(rownames(metaB_ER),annot$Probe_Id)
annot = annot[match(illumina_probes,annot$Probe_Id),]
metaB_ER = metaB_ER[match(illumina_probes,rownames(metaB_ER)),]
metaB_ER$Genes = as.character(annot$ILMN_Gene)
metaB_ER = setDT(metaB_ER)[, lapply(.SD, mean), by = Genes]
metaB_ER = as.data.frame(metaB_ER)
rownames(metaB_ER) = metaB_ER[,1]
metaB_ER = metaB_ER[,-1]

## HER2
clinHER_metaB = metaB_clin[metaB_clin$subtype=="HER2+",]
clinHER_metaB$OS_year = clinHER_metaB$OS_MONTHS/12
clinHER_metaB$OS_Ind = ifelse(clinHER_metaB$OS_STATUS=="DECEASED",1,0)
metaB_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaB_dss1 = metaB_dss[match(clinHER_metaB$PATIENT_ID,metaB_dss$SampleID),]
clinHER_metaB$DSS_days = metaB_dss1$time
clinHER_metaB$DSS_event = metaB_dss1$status
metaB_HER = metaB[,match(clinHER_metaB$PATIENT_ID,colnames(metaB))]
annot = read.delim("~/Illumina/HumanHT-12_V4_0_R2_15002873_B.txt",header=T)
annot = annot[,c(14,5,12)]
annot$Symbol = as.character(annot$Symbol)
illumina_probes = intersect(rownames(metaB_HER),annot$Probe_Id)
annot = annot[match(illumina_probes,annot$Probe_Id),]
metaB_HER = metaB_HER[match(illumina_probes,rownames(metaB_HER)),]
metaB_HER$Genes = as.character(annot$ILMN_Gene)
metaB_HER = setDT(metaB_HER)[, lapply(.SD, mean), by = Genes]
metaB_HER = as.data.frame(metaB_HER)
rownames(metaB_HER) = metaB_HER[,1]
metaB_HER = metaB_HER[,-1]

## GSVA/ssGSEA
metaB_TN_gsva = gsva(as.matrix(metaB_TN),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
metaB_ER_gsva = gsva(as.matrix(metaB_ER),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
metaB_HER_gsva = gsva(as.matrix(metaB_HER),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)

### 2) MB Validation
metaBV = fread("~/Research/Metabric/_ega-box-04_validation_ExpressionMatrix.txt",header=F,sep=" ",data.table=F)
rownames(metaBV) = metaBV[,1]; metaBV = metaBV[,-1]
system("head -1 ~/Research/Metabric/_ega-box-04_validation_ExpressionMatrix.txt > metaBV_samples.txt")
metaBV_samples = scan("metaBV_samples.txt",what="")
colnames(metaBV) = metaBV_samples

temp = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_patient.txt",header=T)
keep_samples = intersect(colnames(metaBV),temp$PATIENT_ID)
temp = temp[match(keep_samples,temp$PATIENT_ID),]
metaBV = metaBV[,match(keep_samples,colnames(metaBV))]
temp1 = read.delim("~/Research/Metabric/brca_metabric/data_clinical_supp_sample.txt",header=T)
temp1 = temp1[match(colnames(metaBV),temp1$PATIENT_ID),]
temp1$ER_STATUS = as.character(temp1$ER_STATUS)
temp1$HER2_STATUS = as.character(temp1$HER2_STATUS)
temp1$PR_STATUS = as.character(temp1$PR_STATUS)
temp1$TN = laply(1:nrow(temp1),function(i) length(grep("-",temp1[i,6:8])))
subtype = ifelse(temp1$ER_STATUS=="+" & temp1$HER2_STATUS=="-","ER+",ifelse(temp1$ER_STATUS=="+" | temp1$ER_STATUS=="-" & temp1$HER2_STATUS=="+","HER2+","TN"))
temp1$subtype = subtype
temp2 = cbind(temp,temp1[,-c(1:2)])
metaBV_clin = temp2

## TN
clinTN_metaBV = metaBV_clin[metaBV_clin$TN=="3",]
clinTN_metaBV$OS_year = clinTN_metaBV$OS_MONTHS/12
clinTN_metaBV$OS_Ind = ifelse(clinTN_metaBV$OS_STATUS=="DECEASED",1,0)
metaBV_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaBV_dss1 = metaBV_dss[match(clinTN_metaBV$PATIENT_ID,metaBV_dss$SampleID),]
clinTN_metaBV$DSS_days = metaBV_dss1$time
clinTN_metaBV$DSS_event = metaBV_dss1$status
metaBV_TN = metaBV[,match(clinTN_metaBV$PATIENT_ID,colnames(metaBV))]
metaBV_TN = metaBV_TN[match(illumina_probes,rownames(metaBV_TN)),]
metaBV_TN$Genes = as.character(annot$ILMN_Gene)
metaBV_TN = setDT(metaBV_TN)[, lapply(.SD, mean),  by = Genes]
metaBV_TN = as.data.frame(metaBV_TN)
rownames(metaBV_TN) = metaBV_TN[,1]
metaBV_TN = metaBV_TN[,-1]
temp_imp = impute.knn(as.matrix(metaBV_TN))
metaBV_TN = temp_imp$data

# ER
clinER_metaBV = metaBV_clin[metaBV_clin$subtype=="ER+",]
clinER_metaBV$OS_year = clinER_metaBV$OS_MONTHS/12
clinER_metaBV$OS_Ind = ifelse(clinER_metaBV$OS_STATUS=="DECEASED",1,0)
metaBV_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaBV_dss1 = metaBV_dss[match(clinER_metaBV$PATIENT_ID,metaBV_dss$SampleID),]
clinER_metaBV$DSS_days = metaBV_dss1$time
clinER_metaBV$DSS_event = metaBV_dss1$status
metaBV_ER = metaBV[,match(clinER_metaBV$PATIENT_ID,colnames(metaBV))]
metaBV_ER = metaBV_ER[match(illumina_probes,rownames(metaBV_ER)),]
metaBV_ER$Genes = as.character(annot$ILMN_Gene)
metaBV_ER = setDT(metaBV_ER)[, lapply(.SD, mean),  by = Genes]
metaBV_ER = as.data.frame(metaBV_ER)
rownames(metaBV_ER) = metaBV_ER[,1]
metaBV_ER = metaBV_ER[,-1]
temp_imp = impute.knn(as.matrix(metaBV_ER))
metaBV_ER = temp_imp$data

# HER
clinHER_metaBV = metaBV_clin[metaBV_clin$subtype=="HER2+",]
clinHER_metaBV$OS_year = clinHER_metaBV$OS_MONTHS/12
clinHER_metaBV$OS_Ind = ifelse(clinHER_metaBV$OS_STATUS=="DECEASED",1,0)
metaBV_dss = read.csv("~/Research/Metabric/brca_metabric/Clinical_Overall_Survival_Data_from_METABRIC.txt",header=T)
metaBV_dss1 = metaBV_dss[match(clinHER_metaBV$PATIENT_ID,metaBV_dss$SampleID),]
clinHER_metaBV$DSS_days = metaBV_dss1$time
clinHER_metaBV$DSS_event = metaBV_dss1$status
metaBV_HER = metaBV[,match(clinHER_metaBV$PATIENT_ID,colnames(metaBV))]
metaBV_HER = metaBV_HER[match(illumina_probes,rownames(metaBV_HER)),]
metaBV_HER$Genes = as.character(annot$ILMN_Gene)
metaBV_HER = setDT(metaBV_HER)[, lapply(.SD, mean),  by = Genes]
metaBV_HER = as.data.frame(metaBV_HER)
rownames(metaBV_HER) = metaBV_HER[,1]
metaBV_HER = metaBV_HER[,-1]
temp_imp = impute.knn(as.matrix(metaBV_HER))
metaBV_HER = temp_imp$data

### ssGSEA
metaBV_TN_gsva = gsva(as.matrix(metaBV_TN),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
metaBV_ER_gsva = gsva(as.matrix(metaBV_ER),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
metaBV_HER_gsva = gsva(as.matrix(metaBV_HER),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)


###########################
### SKCM (Melanoma)
###
#########

print("Preprocessing SKCM (primary)")

### PRIMARY
tcgaMELA = read.delim("~/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/SKCM_primary_norm.txt",header=T,row.names=1,check.names=F)
tcga_mela_gsva = gsva(as.matrix(tcgaMELA),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinMELA = read.delim("/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/SKCM_clinical_primary.txt",header=T,row.names=1)
dim(tcga_clinMELA)


### METASTATIC

print("Preprocessing SKCM (metastatic)")
tcgaMELA_met = fread("~/Research/SKCM/SKCM_met_norm.txt",header=T,data.table=F)
rownames(tcgaMELA_met) = tcgaMELA_met[,1]
tcgaMELA_met = tcgaMELA_met[,-1]
tcga_mela_met_gsva = gsva(as.matrix(tcgaMELA_met),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinMELA_met = read.delim("~/Research/SKCM/SKCM_clinical_met.txt",header=T,row.names=1)
dim(tcga_clinMELA_met)

#################################
#### Immune checkpoint dataset
#### GSE91061
######################

# load("GSE91061_validation.rds")

#################################
## LUAD
#############
print("Preprocessing LUAD")

tcgaLUAD = fread("~/Research/LUAD/LUAD_normalized.txt",header=T,data.table=F)
rownames(tcgaLUAD) = tcgaLUAD[,1]
tcgaLUAD = tcgaLUAD[,-1]
tcga_luad_gsva = gsva(as.matrix(tcgaLUAD),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinLUAD = read.delim("~/Research/LUAD/LUAD_clinical.txt",header=T,row.names=1)
dim(tcga_clinLUAD)

###############################
# LUSC
################
print("Preprocessing LUSC")
tcgaLUSC = fread("~/Research/LUSC/LUSC_normalized.txt",header=T,data.table=F)
rownames(tcgaLUSC) = tcgaLUSC[,1]
tcgaLUSC = tcgaLUSC[,-1]
tcga_lusc_gsva = gsva(as.matrix(tcgaLUSC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinLUSC = read.delim("~/Research/LUSC/LUSC_clinical.txt",header=T,row.names=1)
dim(tcga_clinLUSC)

###############################
# HNSC
################
print("Preprocessing HNSC")
tcgaHNSC = fread("~/Research/HNSC/HNSC_normalized.txt",header=T,data.table=F)
rownames(tcgaHNSC) = tcgaHNSC[,1]
tcgaHNSC = tcgaHNSC[,-1]
tcga_hnsc_gsva = gsva(as.matrix(tcgaHNSC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinHNSC = read.delim("~/Research/HNSC/HNSC_clinical.txt",header=T,row.names=1)
dim(tcga_clinHNSC)

###############################
# COAD
################
print("Preprocessing COAD")
tcgaCOAD = fread("~/Research/COAD/COAD_normalized.txt",header=T,data.table=F)
rownames(tcgaCOAD) = tcgaCOAD[,1]
tcgaCOAD = tcgaCOAD[,-1]
tcga_COAD_gsva = gsva(as.matrix(tcgaCOAD),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinCOAD = read.delim("~/Research/COAD/COAD_clinical.txt",header=T,row.names=1)
dim(tcga_clinCOAD)

###############################
# KIRC
###############
print("Preprocessing KIRC")
tcgaKIRC = fread("~/Research/KIRC/KIRC_normalized.txt",header=T,data.table=F)
rownames(tcgaKIRC) = tcgaKIRC[,1]
tcgaKIRC = tcgaKIRC[,-1]
tcga_KIRC_gsva = gsva(as.matrix(tcgaKIRC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinKIRC = read.delim("~/Research/KIRC/KIRC_clinical.txt",header=T,row.names=1)
dim(tcga_clinKIRC)

###############################
# KIRP
################
print("Preprocessing KIRP")
tcgaKIRP = fread("~/Research/KIRP/KIRP_normalized.txt",header=T,data.table=F)
rownames(tcgaKIRP) = tcgaKIRP[,1]
tcgaKIRP = tcgaKIRP[,-1]
tcga_KIRP_gsva = gsva(as.matrix(tcgaKIRP),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinKIRP = read.delim("~/Research/KIRP/KIRP_clinical.txt",header=T,row.names=1)
dim(tcga_clinKIRP)

###############################
# BLCA
################
print("Preprocessing BLCA")
tcgaBLCA = fread("~/Research/BLCA/BLCA_normalized.txt",header=T,data.table=F)
rownames(tcgaBLCA) = tcgaBLCA[,1]
tcgaBLCA = tcgaBLCA[,-1]
tcga_BLCA_gsva = gsva(as.matrix(tcgaBLCA),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinBLCA = read.delim("~/Research/BLCA/BLCA_clinical.txt",header=T,row.names=1)
dim(tcga_clinBLCA)

###############################
# LIHC
################
print("Preprocessing LIHC")
tcgaLIHC = fread("~/Research/LIHC/LIHC_normalized.txt",header=T,data.table=F)
rownames(tcgaLIHC) = tcgaLIHC[,1]
tcgaLIHC = tcgaLIHC[,-1]
tcga_LIHC_gsva = gsva(as.matrix(tcgaLIHC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinLIHC = read.delim("~/Research/LIHC/LIHC_clinical.txt",header=T,row.names=1)
dim(tcga_clinLIHC)

###############################
# OV
################
print("Preprocessing OV")
tcgaOV = fread("~/Research/OV/OV_normalized.txt",header=T,data.table=F)
rownames(tcgaOV) = tcgaOV[,1]
tcgaOV = tcgaOV[,-1]
tcga_OV_gsva = gsva(as.matrix(tcgaOV),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinOV = read.delim("~/Research/OV/OV_clinical.txt",header=T,row.names=1)
dim(tcga_clinOV)

###############################
# CESC
################
print("Preprocessing CESC")
tcgaCESC = fread("~/Research/CESC/CESC_normalized.txt",header=T,data.table=F)
rownames(tcgaCESC) = tcgaCESC[,1]
tcgaCESC = tcgaCESC[,-1]
tcga_CESC_gsva = gsva(as.matrix(tcgaCESC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinCESC = read.delim("~/Research/CESC/CESC_clinical.txt",header=T,row.names=1)
dim(tcga_clinCESC)

###############################
# GBM
################
print("Preprocessing GBM")
tcgaGBM = fread("~/Research/GBM/GBM_normalized.txt",header=T,data.table=F)
rownames(tcgaGBM) = tcgaGBM[,1]
tcgaGBM = tcgaGBM[,-1]
tcga_GBM_gsva = gsva(as.matrix(tcgaGBM),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinGBM = read.delim("~/Research/GBM/GBM_clinical.txt",header=T,row.names=1)
dim(tcga_clinGBM)

###############################
# UCEC
################
print("Preprocessing UCEC")
tcgaUCEC = fread("~/Research/UCEC/UCEC_normalized.txt",header=T,data.table=F)
rownames(tcgaUCEC) = tcgaUCEC[,1]
tcgaUCEC = tcgaUCEC[,-1]
tcga_UCEC_gsva = gsva(as.matrix(tcgaUCEC),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinUCEC = read.delim("~/Research/UCEC/UCEC_clinical.txt",header=T,row.names=1)
dim(tcga_clinUCEC)

###############################
# PRAD
################
print("Preprocessing PRAD")
tcgaPRAD = fread("~/Research/PRAD/PRAD_normalized.txt",header=T,data.table=F)
rownames(tcgaPRAD) = tcgaPRAD[,1]
tcgaPRAD = tcgaPRAD[,-1]
tcga_PRAD_gsva = gsva(as.matrix(tcgaPRAD),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinPRAD = read.delim("~/Research/PRAD/PRAD_clinical.txt",header=T,row.names=1)
dim(tcga_clinPRAD)

###############################
## LGG
##################
print("Preprocessing LGG")
tcgaLGG = fread("~/Research/LGG/LGG_normalized.txt",header=T,data.table=F)
rownames(tcgaLGG) = tcgaLGG[,1]
tcgaLGG = tcgaLGG[,-1]
tcga_LGG_gsva = gsva(as.matrix(tcgaLGG),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinLGG = read.delim("~/Research/LGG/LGG_clinical.txt",header=T,row.names=1)
dim(tcga_clinLGG)

###############################
# STAD
################
print("Preprocessing STAD")
tcgaSTAD = fread("~/Research/STAD/STAD_normalized.txt",header=T,data.table=F)
rownames(tcgaSTAD) = tcgaSTAD[,1]
tcgaSTAD = tcgaSTAD[,-1]
tcga_STAD_gsva = gsva(as.matrix(tcgaSTAD),gset.idx.list=imm$genesets,method="ssgsea",kcdf="Gaussian",min.sz=1,max.sz=1000)
tcga_clinSTAD = read.delim("~/Research/STAD/STAD_clinical.txt",header=T,row.names=1)
dim(tcga_clinSTAD)

save.image("validation_data.rds")
