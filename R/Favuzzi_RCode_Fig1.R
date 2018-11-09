### Favuzzi et al.: FIGURE 1
#Packages
source("https://bioconductor.org/biocLite.R")
biocLite("cummeRbund")
library(cummeRbund)
library(ggplot2)
library(biomaRt)
library(matrixStats) #rowMaxs function
library(gplots) #heatmap.2 function
#Create database connection
cuff <- readCufflinks(dir = getwd(), 
                      dbFile = "cuffData.db", gtfFile = NULL,
                      runInfoFile  = "run.info", 
                      repTableFile = "read_groups.info",
                      geneFPKM     = "genes.fpkm_tracking",
                      geneDiff     = "gene_exp.diff", 
                      geneCount    = "genes.count_tracking", 
                      geneRep      = "genes.read_group_tracking",
                      isoformFPKM  = "isoforms.fpkm_tracking", 
                      isoformDiff  = "isoform_exp.diff", 
                      isoformCount = "isoforms.count_tracking",
                      isoformRep   = "isoform.read_group_tracking",
                      TSSFPKM      = "tss_groups.fpkm_tracking", 
                      TSSDiff      = "tss_group_exp.diff", 
                      TSSCount     = "tss_groups.count_tracking", 
                      TSSRep       = "tss_groups.read_group_tracking",
                      CDSFPKM      = "cds.fpkm_tracking", 
                      CDSExpDiff   = "cds_exp.diff",
                      CDSCount     = "cds.count_tracking",
                      CDSRep       = "cds.read_group_tracking",
                      CDSDiff      = "cds.diff",
                      promoterFile = "promoters.diff",
                      splicingFile = "splicing.diff",
                      varModelFile = "var_model.info",
                      driver = "SQLite",      #Driver for backend database.
                      genome = "mm9",         #indicates to which genome build the .gtf annotations belong.
                      rebuild = TRUE, verbose = FALSE)
#To reopen the database connection
cuff <- readCufflinks(rebuild=FALSE)
#Matrix of FPKM values
gene_mat_fpkm<-fpkmMatrix(genes(cuff))
#Matrix of differential expression
gene_dif_ch10_pv10<-diffData(genes(cuff),"ch10","pv10")
gene_dif_ch10_ss10<-diffData(genes(cuff),"ch10","ss10")
gene_dif_ss10_pv10<-diffData(genes(cuff),"ss10","pv10")
gene_dif_ch8_pv5<-diffData(genes(cuff),"ch8","pv5")
gene_dif_ch8_ss5<-diffData(genes(cuff),"ch8","ss5")
gene_dif_ss5_pv5<-diffData(genes(cuff),"ss5","pv5")
gene_dif_ch8_ch10<-diffData(genes(cuff),"ch8","ch10")
gene_dif_pv5_pv10<-diffData(genes(cuff),"pv5","pv10")
gene_dif_ss5_ss10<-diffData(genes(cuff),"ss5","ss10")
gene_dif_pv5_ch10<-diffData(genes(cuff),"pv5","ch10")
gene_dif_ss5_ch10<-diffData(genes(cuff),"ss5","ch10")
gene_dif_ch8_pv10<-diffData(genes(cuff),"ch8","pv10")
gene_dif_ss5_pv10<-diffData(genes(cuff),"ss5","pv10")
gene_dif_ch8_ss10<-diffData(genes(cuff),"ch8","ss10")
gene_dif_pv5_ss10<-diffData(genes(cuff),"pv5","ss10")
gene_dif_inp0_ch10<-diffData(genes(cuff),"inp0","ch10")
gene_dif_inp0_pv10<-diffData(genes(cuff),"inp0","pv10")
gene_dif_inp0_ss10<-diffData(genes(cuff),"inp0","ss10")
gene_dif_pyr_ch10<-diffData(genes(cuff),"pyr","ch10")
gene_dif_pyr_pv10<-diffData(genes(cuff),"pyr","pv10")
gene_dif_pyr_ss10<-diffData(genes(cuff),"pyr","ss10")
gene_dif_plp_ch10<-diffData(genes(cuff),"plp","ch10")
gene_dif_plp_pv10<-diffData(genes(cuff),"plp","pv10")
gene_dif_plp_ss10<-diffData(genes(cuff),"plp","ss10")
### HEATMAP Figure 1B
##Calculation of log2 (specificity ratio) matrix
gene_mat_logsr_p10 <-data.frame(
  rownames(gene_mat_fpkm),
  logsp_ch10 = log2(gene_mat_fpkm[c(2)]/rowMaxs(as.matrix(gene_mat_fpkm[c(3,6)]))),
  logsp_pv10 = log2(gene_mat_fpkm[c(3)]/rowMaxs(as.matrix(gene_mat_fpkm[c(2,6)]))),
  logsp_ss10 = log2(gene_mat_fpkm[c(6)]/rowMaxs(as.matrix(gene_mat_fpkm[c(2,3)])))
)
names(gene_mat_logsr_p10) <- c("X","logsp_ch10","logsp_pv10","logsp_ss10")
gene_mat_logsr_p10[is.na(gene_mat_logsr_p10)] <- 0
gene_mat_spec <-gene_mat_logsr_p10
##Threshold values
v_fpkm = 5
v_fdr = 0.05
##Create Mart query (mm9)
listMarts(host='may2012.archive.ensembl.org')
ensembl67 = useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', 
                    dataset='mmusculus_gene_ensembl')
ensembl67_filters <-listFilters(ensembl67)
ensembl67_attributes <-listAttributes(ensembl67)
ensembl67_datasets <-listDatasets(ensembl67)
ensembl67_marts <-listMarts(ensembl67)
##Chandelier cell population at P10: subsetting genes
#fpkm threshold
gene_cel_fpkm <-gene_mat_fpkm[(gene_mat_fpkm$ch10>=v_fpkm),]
gene_cel_fpkm_ <-gene_mat_spec[which(gene_mat_spec$X %in% rownames(gene_cel_fpkm)),]
#qvalue ch10_pv10
gene_dif_ch10_pv10_qv <-gene_dif_ch10_pv10[(gene_dif_ch10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_[
  which(gene_cel_fpkm_$X %in% gene_dif_ch10_pv10_qv$gene_id),]
#qvalue ch10_ss10
gene_dif_ch10_ss10_qv <-gene_dif_ch10_ss10[(gene_dif_ch10_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ch10_ss10_qv$gene_id),]
#qvalue pyr_ch10 
gene_dif_pyr_ch10_qv <-gene_dif_pyr_ch10[(gene_dif_pyr_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pyr_ch10_qv$gene_id),]
#qvalue plp_ch10 
gene_dif_plp_ch10_qv <-gene_dif_plp_ch10[(gene_dif_plp_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_plp_ch10_qv$gene_id),]
#log2specificity threshold
gene_cel_spec <-gene_mat_spec[(gene_mat_spec$logsp_ch10>=0),]
gene_cel_fpkm_qv_spec <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_cel_spec$X),]
gene_cel_deg<-gene_cel_fpkm_qv_spec
#assigning ensembl names
deg_cel <-gene_cel_deg[c(1,2)]
ens_cel <-deg_cel$X
id_cel <-getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","mgi_symbol"),
               values = ens_cel,
               mart=ensembl67)
deg_ch_names <-merge(deg_cel, id_cel,
                     by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE)
deg_ch_names <- deg_ch_names[order(-deg_ch_names[c(2)]),]
#order
deg_ch<-deg_ch_names
deg_ch$order_id <- 1:nrow(deg_ch)
#top 100 genes
deg_ch<-head(deg_ch,100)
##Basket cell population at P10: subsetting genes
#fpkm threshold
gene_cel_fpkm <-gene_mat_fpkm[(gene_mat_fpkm$pv10>=v_fpkm),]
gene_cel_fpkm_ <-gene_mat_spec[which(gene_mat_spec$X %in% rownames(gene_cel_fpkm)),]
#qvalue ch10_pv10
gene_dif_ch10_pv10_qv <-gene_dif_ch10_pv10[(gene_dif_ch10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_[
  which(gene_cel_fpkm_$X %in% gene_dif_ch10_pv10_qv$gene_id),]
#qvalue ss10_pv10
gene_dif_ss10_pv10_qv <-gene_dif_ss10_pv10[(gene_dif_ss10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ss10_pv10_qv$gene_id),]
#qvalue pyr_pv10 
gene_dif_pyr_pv10_qv <-gene_dif_pyr_pv10[(gene_dif_pyr_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pyr_pv10_qv$gene_id),]
#qvalue plp_pv10 
gene_dif_plp_pv10_qv <-gene_dif_plp_pv10[(gene_dif_plp_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_plp_pv10_qv$gene_id),]
#log2specificity threshold
gene_cel_spec <-gene_mat_spec[(gene_mat_spec$logsp_pv10>=0),]
gene_cel_fpkm_qv_spec <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_cel_spec$X),]
gene_cel_deg<-gene_cel_fpkm_qv_spec
#assigning ensembl names
deg_cel <-gene_cel_deg[c(1,3)]
ens_cel <-deg_cel$X
id_cel <-getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","mgi_symbol"),
               values = ens_cel,
               mart=ensembl67)
deg_pv_names <-merge(deg_cel, id_cel,
                     by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE)
deg_pv_names <- deg_pv_names[order(-deg_pv_names[c(2)]),]
#order
deg_pv<-deg_pv_names
deg_pv$order_id <- (nrow(deg_ch)+1):((nrow(deg_ch))+nrow(deg_pv))
#top 100 genes
deg_pv<-head(deg_pv,100)
##Somatostatin cell population at P10: subsetting genes
#fpkm threshold
gene_cel_fpkm <-gene_mat_fpkm[(gene_mat_fpkm$ss10>=v_fpkm),]
gene_cel_fpkm_ <-gene_mat_spec[which(gene_mat_spec$X %in% rownames(gene_cel_fpkm)),]
#qvalue ch10_ss10
gene_dif_ch10_ss10_qv <-gene_dif_ch10_ss10[(gene_dif_ch10_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_[
  which(gene_cel_fpkm_$X %in% gene_dif_ch10_ss10_qv$gene_id),]
#qvalue ss10_pv10
gene_dif_ss10_pv10_qv <-gene_dif_ss10_pv10[(gene_dif_ss10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ss10_pv10_qv$gene_id),]
#qvalue pyr_ss10 
gene_dif_pyr_ss10_qv <-gene_dif_pyr_ss10[(gene_dif_pyr_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pyr_ss10_qv$gene_id),]
#qvalue plp_ss10 
gene_dif_plp_ss10_qv <-gene_dif_plp_ss10[(gene_dif_plp_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_plp_ss10_qv$gene_id),]
#log2specificity threshold
gene_cel_spec <-gene_mat_spec[(gene_mat_spec$logsp_ss10>=0),]
gene_cel_fpkm_qv_spec <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_cel_spec$X),]
gene_cel_deg<-gene_cel_fpkm_qv_spec
#assigning ensembl names
deg_cel <-gene_cel_deg[c(1,4)]
ens_cel <-deg_cel$X
id_cel <-getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","mgi_symbol"),
               values = ens_cel,
               mart=ensembl67)
deg_ss_names <-merge(deg_cel, id_cel,
                     by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE)
deg_ss_names <- deg_ss_names[order(-deg_ss_names[c(2)]),]
#order
deg_ss<-deg_ss_names
deg_ss$order_id <- (nrow(deg_ch)+nrow(deg_pv)+1):((nrow(deg_ch)+nrow(deg_pv))+nrow(deg_ss))
#top 100 genes
deg_ss<-head(deg_ss,100)
##Merge dataframes
deg_merge<-merge(deg_ch,deg_pv,by=c("X","mgi_symbol","order_id"),all=TRUE)
deg_merge<-merge(deg_merge,deg_ss,by=c("X","mgi_symbol","order_id"),all=TRUE)
deg_merge<-deg_merge[c(1:3)]
#intersection
deg_3in <-merge(deg_merge,gene_mat_spec,by="X",all.x=TRUE)
#sort
deg_3in <- deg_3in[order(deg_3in$order_id),] 
#maximum and minimum
funct_inf2max <- function(x)
{
  for (i in 2:ncol(x)){
    x[,i][is.infinite(x[,i]) & x[,i] > 0] = max(x[,i][is.finite(x[,i])])
  }
  return(x)
}
funct_inf2min <- function(x)
{
  for (i in 2:ncol(x)){
    x[,i][is.infinite(x[,i]) & x[,i] < 0] = min(x[,i][is.finite(x[,i])])
  }
  return(x)
}
deg_3in<-funct_inf2max(deg_3in)
deg_3in<-funct_inf2min(deg_3in)
##Matrix for heatmap
hmap<-deg_3in
row.names(hmap) <- hmap$X
hmap<-hmap[c(6,5,4)]
hmap_matrix<-data.matrix(hmap)
break_hm <-c(seq(-6,6,by=0.1))
color_hm <-colorRampPalette(c("grey","white","maroon"))(120)
#heatmap
heatmap.2(hmap_matrix, Colv=FALSE, Rowv=TRUE, dendrogram='row', labRow=FALSE,
          col=color_hm, breaks=break_hm, scale="none", trace="none", 
          key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="log2sr")

### GO ANALYSIS Figure 1C
#(Refer to script for Figure S3D)

### HEATMAPS Figure 1E
##Calculation of log2 (specificity ratio) matrix
gene_mat_logsr_p10 <-data.frame(
  rownames(gene_mat_fpkm),
  logsp_ch10 = log2(gene_mat_fpkm[c(2)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-2,-1)]))),
  logsp_pv10 = log2(gene_mat_fpkm[c(3)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-3,-4)]))),
  logsp_ss10 = log2(gene_mat_fpkm[c(6)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-6,-5)])))
)
names(gene_mat_logsr_p10) <- c("X","logsp_ch10","logsp_pv10","logsp_ss10")
gene_mat_logsr_p10[is.na(gene_mat_logsr_p10)] <- 0
gene_mat_spec <-gene_mat_logsr_p10
##Calculation of upregulation score
gene_mat_upreg <-data.frame(
  rownames(gene_mat_fpkm),
  ch8_ch10 = log2(gene_mat_fpkm[c(2)]/gene_mat_fpkm[c(1)]),
  pv5_pv10 = log2(gene_mat_fpkm[c(3)]/gene_mat_fpkm[c(4)]),
  ss5_ss10 = log2(gene_mat_fpkm[c(6)]/gene_mat_fpkm[c(5)])
)
names(gene_mat_upreg) <- c("X","ch8_ch10","pv5_pv10","ss5_ss10")
gene_mat_upreg[is.na(gene_mat_upreg)] <- 0
##Threshold values
v_fpkm = 10
v_fdr = 0.05
##Somatostatin cell population: subsetting genes
#fpkm threshold
gene_cel_fpkm <-gene_mat_fpkm[(gene_mat_fpkm$ss10>=v_fpkm),]
gene_cel_fpkm_ <-gene_mat_spec[which(gene_mat_spec$X %in% rownames(gene_cel_fpkm)),]
#qvalue ch10_ss10
gene_dif_ch10_ss10_qv <-gene_dif_ch10_ss10[(gene_dif_ch10_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_[
  which(gene_cel_fpkm_$X %in% gene_dif_ch10_ss10_qv$gene_id),]
#qvalue ss10_pv10
gene_dif_ss10_pv10_qv <-gene_dif_ss10_pv10[(gene_dif_ss10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ss10_pv10_qv$gene_id),]
#qvalue ss5_ss10
gene_dif_ss5_ss10_qv <-gene_dif_ss5_ss10[(gene_dif_ss5_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ss5_ss10_qv$gene_id),]
#qvalue ch8_ss10
gene_dif_ch8_ss10_qv <-gene_dif_ch8_ss10[(gene_dif_ch8_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ch8_ss10_qv$gene_id),]
#qvalue pv5_ss10
gene_dif_pv5_ss10_qv <-gene_dif_pv5_ss10[(gene_dif_pv5_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pv5_ss10_qv$gene_id),]
#qvalue inp0_ss10 
gene_dif_inp0_ss10_qv <-gene_dif_inp0_ss10[(gene_dif_inp0_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_inp0_ss10_qv$gene_id),]
#qvalue pyr_ss10 
gene_dif_pyr_ss10_qv <-gene_dif_pyr_ss10[(gene_dif_pyr_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pyr_ss10_qv$gene_id),]
#qvalue plp_ss10 
gene_dif_plp_ss10_qv <-gene_dif_plp_ss10[(gene_dif_plp_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_plp_ss10_qv$gene_id),]
#log2specificity threshold
gene_cel_spec <-gene_mat_spec[(gene_mat_spec$logsp_ss10>=0),]
gene_cel_fpkm_qv_spec <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_cel_spec$X),]
#upregulation P5-P10 threshold
gene_upreg_cel <-gene_mat_upreg[(gene_mat_upreg$ss5_ss10>=0),]
gene_cel_deg <-gene_cel_fpkm_qv_spec[
  which(gene_cel_fpkm_qv_spec$X %in% gene_upreg_cel$X),]
#Import database from Hinojosa et al., 2018: "AH_nex10_nkx10"
#threshold values
AH_nex10_nkx10_ <-AH_nex10_nkx10[(AH_nex10_nkx10$R.fold<1.5),]
AH_nex10_nkx10_subset <-AH_nex10_nkx10_[(AH_nex10_nkx10_$q.value<0.05),]
#substract genes
gene_cel_deg <-gene_cel_deg[
  !(gene_cel_deg$X %in% AH_nex10_nkx10_subset$probeset.ID),]
gene_ss_deg <- gene_cel_deg
##Basket cell population: subsetting genes
#fpkm threshold
gene_cel_fpkm <-gene_mat_fpkm[(gene_mat_fpkm$pv10>=v_fpkm),]
gene_cel_fpkm_ <-gene_mat_spec[which(gene_mat_spec$X %in% rownames(gene_cel_fpkm)),]
#qvalue ch10_pv10
gene_dif_ch10_pv10_qv <-gene_dif_ch10_pv10[(gene_dif_ch10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_[
  which(gene_cel_fpkm_$X %in% gene_dif_ch10_pv10_qv$gene_id),]
#qvalue ss10_pv10
gene_dif_ss10_pv10_qv <-gene_dif_ss10_pv10[(gene_dif_ss10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ss10_pv10_qv$gene_id),]
#qvalue pv5_pv10
gene_dif_pv5_pv10_qv <-gene_dif_pv5_pv10[(gene_dif_pv5_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pv5_pv10_qv$gene_id),]
#qvalue ch8_pv10
gene_dif_ch8_pv10_qv <-gene_dif_ch8_pv10[(gene_dif_ch8_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ch8_pv10_qv$gene_id),]
#qvalue ss5_pv10
gene_dif_ss5_pv10_qv <-gene_dif_ss5_pv10[(gene_dif_ss5_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ss5_pv10_qv$gene_id),]
#qvalue inp0_pv10 
gene_dif_inp0_pv10_qv <-gene_dif_inp0_pv10[(gene_dif_inp0_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_inp0_pv10_qv$gene_id),]
#qvalue pyr_pv10 
gene_dif_pyr_pv10_qv <-gene_dif_pyr_pv10[(gene_dif_pyr_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pyr_pv10_qv$gene_id),]
#qvalue plp_pv10 
gene_dif_plp_pv10_qv <-gene_dif_plp_pv10[(gene_dif_plp_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_plp_pv10_qv$gene_id),]
#log2specificity threshold
gene_cel_spec <-gene_mat_spec[(gene_mat_spec$logsp_pv10>=0),]
gene_cel_fpkm_qv_spec <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_cel_spec$X),]
#upregulation P5-P10 threshold
gene_upreg_cel <-gene_mat_upreg[(gene_mat_upreg$pv5_pv10>=0),]
gene_cel_deg <-gene_cel_fpkm_qv_spec[
  which(gene_cel_fpkm_qv_spec$X %in% gene_upreg_cel$X),]
gene_pv_deg <- gene_cel_deg
##Chandelier cell population: subsetting genes
#fpkm threshold
gene_cel_fpkm <-gene_mat_fpkm[(gene_mat_fpkm$ch10>=v_fpkm),]
gene_cel_fpkm_ <-gene_mat_spec[which(gene_mat_spec$X %in% rownames(gene_cel_fpkm)),]
#qvalue ch10_pv10
gene_dif_ch10_pv10_qv <-gene_dif_ch10_pv10[(gene_dif_ch10_pv10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_[
  which(gene_cel_fpkm_$X %in% gene_dif_ch10_pv10_qv$gene_id),]
#qvalue ch10_ss10
gene_dif_ch10_ss10_qv <-gene_dif_ch10_ss10[(gene_dif_ch10_ss10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ch10_ss10_qv$gene_id),]
#qvalue ch8_ch10
gene_dif_ch8_ch10_qv <-gene_dif_ch8_ch10[(gene_dif_ch8_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ch8_ch10_qv$gene_id),]
#qvalue pv5_ch10
gene_dif_pv5_ch10_qv <-gene_dif_pv5_ch10[(gene_dif_pv5_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pv5_ch10_qv$gene_id),]
#qvalue ss5_ch10
gene_dif_ss5_ch10_qv <-gene_dif_ch8_ch10[(gene_dif_ch8_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_ss5_ch10_qv$gene_id),]
#qvalue inp0_ch10
gene_dif_inp0_ch10_qv <-gene_dif_inp0_ch10[(gene_dif_inp0_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_inp0_ch10_qv$gene_id),]
#qvalue pyr_ch10 
gene_dif_pyr_ch10_qv <-gene_dif_pyr_ch10[(gene_dif_pyr_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_pyr_ch10_qv$gene_id),]
#qvalue plp_ch10 
gene_dif_plp_ch10_qv <-gene_dif_plp_ch10[(gene_dif_plp_ch10$q_value<v_fdr),]
gene_cel_fpkm_qv <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_dif_plp_ch10_qv$gene_id),]
#log2specificity threshold
gene_cel_spec <-gene_mat_spec[(gene_mat_spec$logsp_ch10>=0),]
gene_cel_fpkm_qv_spec <-gene_cel_fpkm_qv[
  which(gene_cel_fpkm_qv$X %in% gene_cel_spec$X),]
#upregulation P5-P10 threshold
gene_upreg_cel <-gene_mat_upreg[(gene_mat_upreg$ch8_ch10>=0),]
gene_cel_deg <-gene_cel_fpkm_qv_spec[
  which(gene_cel_fpkm_qv_spec$X %in% gene_upreg_cel$X),]
gene_ch_deg <- gene_cel_deg
##Dataframe for heatmaps
gene_hm_logsr_p10 <-data.frame(
  rownames(gene_mat_fpkm),
  logsp_in = log2(gene_mat_fpkm[c(7)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-7)]))),
  logsp_ch8 = log2(gene_mat_fpkm[c(1)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-2,-1)]))),
  logsp_ch10 = log2(gene_mat_fpkm[c(2)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-2,-1)]))),
  logsp_pv5 = log2(gene_mat_fpkm[c(4)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-3,-4)]))),
  logsp_pv10 = log2(gene_mat_fpkm[c(3)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-3,-4)]))),
  logsp_ss5 = log2(gene_mat_fpkm[c(5)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-6,-5)]))),
  logsp_ss10 = log2(gene_mat_fpkm[c(6)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-6,-5)]))),
  logsp_pyr = log2(gene_mat_fpkm[c(8)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-8)]))),
  logsp_plp = log2(gene_mat_fpkm[c(9)]/rowMaxs(as.matrix(gene_mat_fpkm[c(-9)])))
)
names(gene_hm_logsr_p10) <- c("X","logsp_inp0","logsp_ch8","logsp_ch10",
                              "logsp_pv5","logsp_pv10","logsp_ss5","logsp_ss10","logsp_pyr","logsp_plp")
gene_hm_logsr_p10[is.na(gene_hm_logsr_p10)] <- 0
gene_hm_spec <-gene_hm_logsr_p10
##Create Mart query (mm9)
listMarts(host='may2012.archive.ensembl.org')
ensembl67 = useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', 
                    dataset='mmusculus_gene_ensembl')
ensembl67_filters <-listFilters(ensembl67)
ensembl67_attributes <-listAttributes(ensembl67)
ensembl67_datasets <-listDatasets(ensembl67)
ensembl67_marts <-listMarts(ensembl67)
##Assign gene names
#somatostatin population
deg_ss <-gene_hm_spec[which(gene_hm_spec$X %in% gene_ss_deg$X),]
ens_ss <-deg_ss$X
id_ss <-getBM(filters = "ensembl_gene_id",
              attributes = c("ensembl_gene_id","mgi_symbol"),
              values = ens_ss,
              mart=ensembl67)
gene_ss_deg_names <-merge(deg_ss, id_ss,
                          by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE)
gene_ss_deg_names <- gene_ss_deg_names[order(-gene_ss_deg_names[c(8)]),]
gene_ss_deg_hm <-gene_ss_deg_names[c(2:10)]
rownames(gene_ss_deg_hm)<-gene_ss_deg_names$mgi_symbol
#basket population
deg_pv <-gene_hm_spec[which(gene_hm_spec$X %in% gene_pv_deg$X),]
ens_pv <-deg_pv$X
id_pv <-getBM(filters = "ensembl_gene_id",
              attributes = c("ensembl_gene_id","mgi_symbol"),
              values = ens_pv,
              mart=ensembl67)
gene_pv_deg_names <-merge(deg_pv, id_pv,
                          by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE)
gene_pv_deg_names <- gene_pv_deg_names[order(-gene_pv_deg_names[c(6)]),]
gene_pv_deg_hm <-gene_pv_deg_names[c(2:10)]
rownames(gene_pv_deg_hm)<-gene_pv_deg_names$mgi_symbol
#chandelier population
deg_ch <-gene_hm_spec[which(gene_hm_spec$X %in% gene_ch_deg$X),]
ens_ch <-deg_ch$X
id_ch <-getBM(filters = "ensembl_gene_id",
              attributes = c("ensembl_gene_id","mgi_symbol"),
              values = ens_ch,
              mart=ensembl67)
gene_ch_deg_names <-merge(deg_ch, id_ch,
                          by.x = "X", by.y = "ensembl_gene_id", all.x = TRUE)
gene_ch_deg_names <- gene_ch_deg_names[order(-gene_ch_deg_names[c(4)]),]
gene_ch_deg_hm <-gene_ch_deg_names[c(2:10)]
rownames(gene_ch_deg_hm)<-gene_ch_deg_names$mgi_symbol
##Matrix for heatmap
#ss
gene_deg_hm<-gene_ss_deg_hm[1:25,]
break_hm <-c(seq(-4,4,by=0.1))
color_hm <-colorRampPalette(c("grey25","white","#FF8C5A"))(80)
#pv
gene_deg_hm<-gene_pv_deg_hm[1:25,]
break_hm <-c(seq(-4,4,by=0.1))
color_hm <-colorRampPalette(c("grey25","white","#189FC9"))(80)
#ch
gene_deg_hm<-gene_ch_deg_hm[1:25,]
break_hm <-c(seq(-6,6,by=0.1))
color_hm <-colorRampPalette(c("grey25","white","#E0141A"))(120)
#matrix
gene_deg_hm_matrix<-data.matrix(gene_deg_hm[c(1,6,7,4,5,2,3,8,9)])
#heatmaps
heatmap.2(gene_deg_hm_matrix, Colv=FALSE, Rowv=FALSE, dendrogram='none',
          col=color_hm, breaks=break_hm, scale="none", trace="none",
          key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=45, cexCol=0.6, cexRow=0.6,
          colsep=1:ncol(gene_deg_hm_matrix),rowsep=1:nrow(gene_deg_hm_matrix),
          sepcolor="white")
####