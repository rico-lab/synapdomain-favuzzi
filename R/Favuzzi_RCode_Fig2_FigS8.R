### Favuzzi et al.: FIGURE 2 & SUPPLEMENTARY FIGURE S8
#Packages
library(cummeRbund)
library(ggplot2)
library(gplots) #heatmap.2 function
library(matrixStats) #rowMaxs function
library(biomaRt)
library(ggtern) #ternary diagrams
#Input matrix
cuff <- readCufflinks(rebuild=FALSE)
gene.matrix<-fpkmMatrix(genes(cuff))
### HEATMAPS
#Figure 2A: Heatmap of SST+ interneuron-specific genes
#Cbln4=ENSMUSG00000067578; Igsf21=ENSMUSG00000040972; Cd59a=ENSMUSG00000032679; Ptpru=ENSMUSG00000028909; Pcdh18=ENSMUSG00000037892;
specific_genes <- c("ENSMUSG00000067578", "ENSMUSG00000040972", "ENSMUSG00000032679", "ENSMUSG00000028909", "ENSMUSG00000037892")
heatmap_genes <- subset(gene.matrix, rownames(gene.matrix) %in% specific_genes)
heatmap_genes<-heatmap_genes[c(5,4,2,1,3),]
rownames(heatmap_genes) <- c("Cbln4","Igsf21","Cd59a","Ptpru", "Pcdh18")
heatmap_genes_df <- data.frame(heatmap_genes[,c(5,6,4,3,1,2,7,8,9)])
heatmap_genes_df$expr <-log2(heatmap_genes_df$ss10+1)
heatmap_genes_df$upreg <-log2(heatmap_genes_df$ss10/heatmap_genes_df$ss5)
heatmap_genes_df[c(12)] <-log2(heatmap_genes_df[c(2)]/rowMaxs(as.matrix(heatmap_genes_df[c(3:9)])))
names(heatmap_genes_df)[(12)]<-"spec"
#expression level
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(10,10)])
break_hm <-c(seq(0,8,by=0.1))
color_hm <-colorRampPalette(c("white","#BA946E"))(80)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="Expression level")
#upregulation score
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(11,11)])
break_hm <-c(seq(0,3,by=0.1))
color_hm <-colorRampPalette(c("white","#0B9B75"))(30)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="Upregulation score")
#specificity ratio (ss)
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(12,12)])
break_hm <-c(seq(0,4,by=0.1))
color_hm <-colorRampPalette(c("white","#FF8C5A"))(40)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="log2(specificity ratio)")
#Figure 2A: Heatmap of PV+ basket cell-specific genes
#Lgals1=ENSMUSG00000068220; Lgi2=ENSMUSG00000039252; Tmem91=ENSMUSG00000061702; Thy1=ENSMUSG00000032011; Cgref1=ENSMUSG00000029161;
specific_genes <- c("ENSMUSG00000068220", "ENSMUSG00000039252", "ENSMUSG00000061702", "ENSMUSG00000032011", "ENSMUSG00000029161")
heatmap_genes <- subset(gene.matrix, rownames(gene.matrix) %in% specific_genes)
heatmap_genes<-heatmap_genes[c(5,3,4,2,1),]
rownames(heatmap_genes) <- c("Lgals1","Lgi2","Tmem91","Thy1","Cgref1")
heatmap_genes_df <- data.frame(heatmap_genes[,c(5,6,4,3,1,2,7,8,9)])
heatmap_genes_df$expr <-log2(heatmap_genes_df$pv10+1)
heatmap_genes_df$upreg <-log2(heatmap_genes_df$pv10/heatmap_genes_df$pv5)
heatmap_genes_df[c(12)] <-log2(heatmap_genes_df[c(4)]/rowMaxs(as.matrix(heatmap_genes_df[c(1:2,5:9)])))
names(heatmap_genes_df)[(12)]<-"spec"
#expression level
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(10,10)])
break_hm <-c(seq(0,8,by=0.1))
color_hm <-colorRampPalette(c("white","#BA946E"))(80)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="Expression level")
#upregegulation score
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(11,11)])
break_hm <-c(seq(0,3,by=0.1))
color_hm <-colorRampPalette(c("white","#0B9B75"))(30)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="Upregulation score")
#specificity ratio (pv)
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(12,12)])
break_hm <-c(seq(0,4,by=0.1))
color_hm <-colorRampPalette(c("white","#189FC9"))(40)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="log2(specificity ratio)")
#Figure 2A: Heatmap of chandelier cell-specific genes
#Hapln1=ENSMUSG00000021613; Thsd7a=ENSMUSG00000032625; Fgf13=ENSMUSG00000031137; Sema3c=ENSMUSG00000028780; Mme=ENSMUSG00000027820;
specific_genes <- c("ENSMUSG00000021613", "ENSMUSG00000032625", "ENSMUSG00000031137", "ENSMUSG00000028780", "ENSMUSG00000027820")
heatmap_genes <- subset(gene.matrix, rownames(gene.matrix) %in% specific_genes)
heatmap_genes<-heatmap_genes[c(1,5,4,3,2),]
rownames(heatmap_genes) <- c("Hapln1","Thsd7a","Fgf13","Sema3c","Mme")
heatmap_genes_df <- data.frame(heatmap_genes[,c(5,6,4,3,1,2,7,8,9)])
heatmap_genes_df$expr <-log2(heatmap_genes_df$ch10+1)
heatmap_genes_df$upreg <-log2(heatmap_genes_df$ch10/heatmap_genes_df$ch8)
heatmap_genes_df[c(12)] <-log2(heatmap_genes_df[c(6)]/rowMaxs(as.matrix(heatmap_genes_df[c(1:4,7:9)])))
names(heatmap_genes_df)[(12)]<-"spec"
#expression level
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(10,10)])
break_hm <-c(seq(0,8,by=0.1))
color_hm <-colorRampPalette(c("white","#BA946E"))(80)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="Expression level")
#upregulation score
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(11,11)])
break_hm <-c(seq(0,3,by=0.1))
color_hm <-colorRampPalette(c("white","#0B9B75"))(30)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="Upregulation score")
#specificity ratio (ch)
heatmap_genes_m<-data.matrix(heatmap_genes_df[c(12,12)])
break_hm <-c(seq(0,6,by=0.1))
color_hm <-colorRampPalette(c("white","#E0141A"))(60)
heatmap.2(heatmap_genes_m, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="log2(specificity ratio)")
### TERNARY DIAGRAMS
##Calculate specficity score for ternary diagrams
#create Mart Query
listMarts(host='may2012.archive.ensembl.org')
ensembl67 = useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', 
                    dataset='mmusculus_gene_ensembl')
ensembl67_filters <-listFilters(ensembl67)
ensembl67_attributes <-listAttributes(ensembl67)
ensembl67_datasets <-listDatasets(ensembl67)
ensembl67_marts <-listMarts(ensembl67)
#assigning ensembl-mgi gene names
gene_matrix_names <-rownames(gene.matrix)
gene_matrix_ensembl <-getBM(filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","mgi_symbol"),
                  values = gene_matrix_names,
                  mart=ensembl67)
gene_matrix_ids <-merge(gene.matrix, gene_matrix_ensembl,
                     by.x = 0, by.y = "ensembl_gene_id", all.x = TRUE)
#aggregate duplicates
gene_matrix_mgi <- gene_matrix_ids[c(11,2:10)]
gene_matrix_mgi <- unique(gene_matrix_mgi)
gene_matrix_mgi_sum <-aggregate(gene_matrix_mgi[, -c(1)],by=list(gene_matrix_mgi$mgi_symbol),FUN="sum",na.rm=TRUE)
gene_matrix_mgi_sum<-unique(gene_matrix_mgi_sum)
gene_matrix_mgi_sum<-gene_matrix_mgi_sum[-1,]
#specificity score function
tissuespecificity<-function(fpkms,relative=FALSE,...){
  fpkms<-t(makeprobs(t(fpkms)))
  d<-diag(ncol(fpkms))
  res<-apply(d,MARGIN=1,function(q){
    JSdistFromP(fpkms,q)
  })
  colnames(res)<-paste(colnames(fpkms),"_spec",sep="")
  if(relative){
    res<-res/max(res)
  }
  1-res
}
#specificity score for P10 samples
genes_P10 <- gene_matrix_mgi_sum[c(1,7,4,3)] #drop other cells and ages
genes_P10 <- unique(genes_P10)
rownames(genes_P10) <- genes_P10$Group.1
genes_P10 <- genes_P10[c(2:4)]
genes_P10_specifscore <- tissuespecificity(genes_P10)
genes_P10_specifscore_sum<-as.data.frame(genes_P10_specifscore)
genes_P10_specifscore_sum$sum <- rowSums(genes_P10_specifscore_sum[c(1:3)])
genes_P10_specifscore_sum_rel<-genes_P10_specifscore_sum
genes_P10_specifscore_sum_rel_max<-genes_P10_specifscore_sum_rel/apply(genes_P10_specifscore_sum_rel,1,max)
genes_P10_specifscore_sum_rel_max<-genes_P10_specifscore_sum_rel_max[c(1:3)]
#specificity score for P5/8 samples
genes_P5 <- gene_matrix_mgi_sum[c(1,6,5,2)] #drop other cells and ages
genes_P5 <- unique(genes_P5)
rownames(genes_P5) <- genes_P5$Group.1
genes_P5 <- genes_P5[c(2:4)]
genes_P5_specifscore <- tissuespecificity(genes_P5)
genes_P5_specifscore_sum<-as.data.frame(genes_P5_specifscore)
genes_P5_specifscore_sum$sum <- rowSums(genes_P5_specifscore_sum[c(1:3)])
genes_P5_specifscore_sum_rel<-genes_P5_specifscore_sum
genes_P5_specifscore_sum_rel_max<-genes_P5_specifscore_sum_rel/apply(genes_P5_specifscore_sum_rel,1,max)
genes_P5_specifscore_sum_rel_max<-genes_P5_specifscore_sum_rel_max[c(1:3)]
#combine dataframes
gene_matrix_specif <- cbind(genes_P10,genes_P10_specifscore_sum_rel_max,genes_P5,genes_P5_specifscore_sum_rel_max)
##Figure 2B: Ternary diagrams
gene_matrix_specif_fpkm5 <-subset(gene_matrix_specif,ss10>=5 | pv10>=5 | ch10>=5)
df_tern <-gene_matrix_specif_fpkm5
df_tern[c(13)] <-log2(rowMaxs(as.matrix(df_tern[c(1:3)]))+1)
df_tern[c(14)] <-log2(rowMaxs(as.matrix(df_tern[c(7:9)]))+1)
names(df_tern)[(13)]<-"log2size_p10"
names(df_tern)[(14)]<-"log2size_p5"
#Cbln protein family
cbln_family <- c("Cbln1", "Cbln2", "Cbln3", "Cbln4")
df_tern_cbln <- subset(df_tern, rownames(df_tern) %in% cbln_family)
df_all <-df_tern_cbln
#Lgi protein family
lgi_family <- c("Lgi1", "Lgi2", "Lgi3")
df_tern_lgi <- subset(df_tern, rownames(df_tern) %in% lgi_family)
df_all <-df_tern_lgi
#Fgf protein family
fgf_family <- c("Fgf1", "Fgf10", "Fgf11", "Fgf12", "Fgf13", "Fgf14", "Fgf15", "Fgf16", "Fgf17", "Fgf18", 
                "Fgf2", "Fgf20", "Fgf21", "Fgf22", "Fgf23", "Fgf3", "Fgf4", "Fgf5", "Fgf6", "Fgf7", "Fgf8", "Fgf9")
df_tern_fgf <- subset(df_tern, rownames(df_tern) %in% fgf_family)
df_all <-df_tern_fgf
#Ternary diagrams
ggtern()+
  geom_point(data=df_all,aes(pv10_spec,ss10_spec,ch10_spec),size=log2size_p10,shape=21,fill="gray",alpha=0.7)+
  geom_point(data=df_all,aes(pv5_spec,ss5_spec,ch8_spec),size=log2size_p5,shape=21,alpha=0.7)+
  geom_text(data=df_all,aes(pv10_spec,ss10_spec,ch10_spec,label=rownames(df_all)),hjust=0, vjust=0)+
  geom_text(data=df_all,aes(pv5_spec,ss5_spec,ch8_spec,label=rownames(df_all)),hjust=0, vjust=0)+
  theme_bw()
##Figure S8: Ternary diagrams
#Cntn protein family
cntn_family <- c("Cntn1", "Cntn2", "Cntn3", "Cntn4", "Cntn5", "Cntn6")
df_tern_cntn <- subset(df_tern, rownames(df_tern) %in% cntn_family)
df_all <-df_tern_cntn
#Cntnap protein family
cntnap_family <- c("Cntnap1", "Cntnap2", "Cntnap3", "Cntnap4", "Cntnap5a", "Cntnap5b", "Cntnap5c")
df_tern_cntnap <- subset(df_tern, rownames(df_tern) %in% cntnap_family)
df_all <-df_tern_cntnap
#Igsf protein family
igsf_family <- c("Igsf1", "Igsf10", "Igsf11", "Igsf21", "Igsf3", "Igsf5", "Igsf6", "Igsf8", "Igsf9", "Igsf9b")
df_tern_igsf <- subset(df_tern, rownames(df_tern) %in% igsf_family)
df_all <-df_tern_igsf
#Nxph protein family
nxph_family <- c("Nxph1", "Nxph2", "Nxph3", "Nxph4")
df_tern_nxph <- subset(df_tern, rownames(df_tern) %in% nxph_family)
df_all <-df_tern_nxph
#Cdh protein family (subgroup 1)
cdh_family1 <- c("Cdh1", "Cdh3", "Cdh4", "Cdh5", "Cdh7", "Cdh10", "Cdh12", "Cdh13", "Cdh18", "Cdh19", "Cdh20")
df_tern_cdh1 <- subset(df_tern, rownames(df_tern) %in% cdh_family1)
df_all <-df_tern_cdh1
#Cdh protein family (subgroup 1)
cdh_family2 <- c("Cdh2", "Cdh6", "Cdh8", "Cdh9", "Cdh11", "Cdh15", "Cdh16", "Cdh17", "Cdh22", "Cdh23", "Cdh24", "Cdh26")
df_tern_cdh2 <- subset(df_tern, rownames(df_tern) %in% cdh_family2)
df_all <-df_tern_cdh2
#Pcdh protein family (subgroup 1)
pcdh_family1 <- c("Pcdh1", "Pcdh8", "Pcdh10", "Pcdh12", "Pcdh15", "Pcdh18","Pcdh19")
df_tern_pcdh1 <- subset(df_tern, rownames(df_tern) %in% pcdh_family1)
df_all <-df_tern_pcdh1
#Pcdh protein family (subgroup 2)
pcdh_family2 <- c("Pcdh7", "Pcdh9", "Pcdh11x", "Pcdh12", "Pcdh17","Pcdh20")
df_tern_pcdh2 <- subset(df_tern, rownames(df_tern) %in% pcdh_family2)
df_all <-df_tern_pcdh2
#IgCAM protein family (subgroup 1)
igcam_family1 <- c("Dscam", "Kirrel", "Kirrel2", "Kirrel3", "Kit", "Kitl", "L1cam", "Nrcam", "Sdk1", "Sdk2")
df_tern_igcam1 <- subset(df_tern, rownames(df_tern) %in% igcam_family1)
df_all <-df_tern_igcam1
#IgCAM protein family (subgroup 2)
igcam_family2 <- c("Cttnbp2", "Cttnbp2nl", "Dscaml1", "Jam2", "Jam3", "Nfasc", "Icam1", "Icam2", "Icam4", "Icam5")
df_tern_igcam2 <- subset(df_tern, rownames(df_tern) %in% igcam_family2)
df_all <-df_tern_igcam2
#Sema protein family (subgroup 1)
sema_family1 <- c("Sema3a", "Sema3b", "Sema3d", "Sema4a", "Sema4d", "Sema4f", "Sema5a", "Sema6b", "Sema6c", "Sema7a")
df_tern_sema1 <- subset(df_tern, rownames(df_tern) %in% sema_family1)
df_all <-df_tern_sema1
#Sema protein family (subgroup 2)
sema_family2 <- c("Sema3c", "Sema3e", "Sema3f", "Sema3g", "Sema4b", "Sema4c", "Sema4g", "Sema5b", "Sema6a", "Sema6d")
df_tern_sema2 <- subset(df_tern, rownames(df_tern) %in% sema_family2)
df_all <-df_tern_sema2
#Ptpr protein family (subgroup 1)
ptpr_family1 <- c("Ptprg", "Ptprk", "Ptprm", "Ptpru", "Ptprz1")
df_tern_ptpr1 <- subset(df_tern, rownames(df_tern) %in% ptpr_family1)
df_all <-df_tern_ptpr1
#Ptpr protein family (subgroup 2)
ptpr_family2 <- c("Ptpre", "Ptprn2", "Ptpro", "Ptprt", "Ptprv")
df_tern_ptpr2 <- subset(df_tern, rownames(df_tern) %in% ptpr_family2)
df_all <-df_tern_ptpr2
#Ptpr protein family (subgroup 3)
ptpr_family3 <- c("Ptprb", "Ptprc", "Ptprcap", "Ptprf", "Ptprn", "Ptprr")
df_tern_ptpr3 <- subset(df_tern, rownames(df_tern) %in% ptpr_family3)
df_all <-df_tern_ptpr3
#Ptpr protein family (subgroup 4)
ptpr_family4 <- c("Ptpra", "Ptprj", "Ptprh", "Ptprq", "Ptprs")
df_tern_ptpr4 <- subset(df_tern, rownames(df_tern) %in% ptpr_family4)
df_all <-df_tern_ptpr4
#Ternary diagrams
ggtern()+
  geom_point(data=df_all,aes(pv10_spec,ss10_spec,ch10_spec),size=log2size_p10,shape=21,fill="gray",alpha=0.7)+
  geom_point(data=df_all,aes(pv5_spec,ss5_spec,ch8_spec),size=log2size_p5,shape=21,alpha=0.7)+
  geom_text(data=df_all,aes(pv10_spec,ss10_spec,ch10_spec,label=rownames(df_all)),hjust=0, vjust=0)+
  geom_text(data=df_all,aes(pv5_spec,ss5_spec,ch8_spec,label=rownames(df_all)),hjust=0, vjust=0)+
  theme_bw()
### Note that the individual dots in the ternary diagrams were coloured according to 
### the colour code explained in the corresponding Figure Legends (Fig.2B and Fig.S8A).
####