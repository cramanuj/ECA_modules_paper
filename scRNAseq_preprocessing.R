## Jan 19, 2019
## Chaitanya Acharya
## Data preprocessing

library(cellrangerRkit);library(data.table)
source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/scripts/r_scripts/useful_functions.R")

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
fwrite(seu_exprs,"BRCA_scRNAseq_all_patients.dat",sep="\t",col.names=T,row.names=T,quote=F)

## Impute

library(devtools)
install_github("Vivianstats/scImpute")
library(scImpute)
seu_imp = scimpute(count_path="/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/BRCA_scRNAseq_all_patients.dat",
                    infile="txt",outfile="txt",out_dir="/Users/ca31/Research/BRCA/Loi_BCTCells_NatMed_071018/metagene_based_final/",labeled = FALSE,drop_thre = 0.5,Kcluster=2,ncores=2)
seu_imp = fread("scimpute_count.txt",header=F,data.table=F)
rownames(seu_imp) = seu_imp[,1]
seu_imp = seu_imp[,-1]
colnames(seu_imp) = colnames(seu_exprs)
seu_imp = seu_imp[match(rownames(seu_exprs),rownames(seu_imp)),]
pdata = data.frame(barcode=colnames(seu_imp),patientid = c(rep(1,5174),rep(2,1137)))
rownames(pdata) = pdata$barcode
seu_imp$Genes = as.factor(as.character(probes_exprs$Genes))
dat = setDT(seu_imp)[, lapply(.SD, max), by = Genes]
dat = as.data.frame(dat)
rownames(dat) = dat[,1]
dat = dat[,-1]
fwrite(dat,"BRCA_scRNAseq_all_patients_gene_level.dat",sep="\t",col.names=T,row.names=T,quote=F)

save.image("brca_preprocessing.rds")
