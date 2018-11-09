### Favuzzi et al.: SUPPLEMENTARY FIGURE S3
#Packages
source("https://bioconductor.org/biocLite.R")
biocLite("cummeRbund")
library(cummeRbund)
library(ggplot2)
library(gplots) #heatmap.2 function
library(matrixStats) #rowMaxs function
library(biomaRt)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GOstats)
### LOAD MATRICES (from Figure S2)
#Reopen the database connection (from Figure S2)
cuff <- readCufflinks(rebuild=FALSE)
#Matrix of FPKM values
gene_mat_fpkm<-fpkmMatrix(genes(cuff))
#Matrix of differential expression
gene_dif_ch10_pv10<-diffData(genes(cuff),"ch10","pv10")
gene_dif_ch10_ss10<-diffData(genes(cuff),"ch10","ss10")
gene_dif_ss10_pv10<-diffData(genes(cuff),"ss10","pv10")
gene_dif_pyr_ch10<-diffData(genes(cuff),"pyr","ch10")
gene_dif_pyr_pv10<-diffData(genes(cuff),"pyr","pv10")
gene_dif_pyr_ss10<-diffData(genes(cuff),"pyr","ss10")
gene_dif_plp_ch10<-diffData(genes(cuff),"plp","ch10")
gene_dif_plp_pv10<-diffData(genes(cuff),"plp","pv10")
gene_dif_plp_ss10<-diffData(genes(cuff),"plp","ss10")
### HEATMAPS
#Figure S3A: Heatmap of well-known cell type-specific genes
#Elfn1=ENSMUSG00000048988; Syt2=ENSMUSG00000026452; Pthlh=ENSMUSG00000048776;
specific_genes <- c("ENSMUSG00000048988", "ENSMUSG00000026452", "ENSMUSG00000048776")
heatmap_genes <- subset(gene_mat_fpkm, rownames(gene_mat_fpkm) %in% specific_genes)
rownames(heatmap_genes) <- c("Syt2","Pthlh","Elfn1")
heatmap_genes_m <- data.matrix(heatmap_genes[,c(5,6,4,3,1,2)])
heatmap_genes_m_log <-log2(heatmap_genes_m+1)
break_hm <-c(seq(-0,4,by=0.1))
color_hm <-colorRampPalette(c("grey","white","maroon"))(40)
heatmap.2(heatmap_genes_m_log, 
          Colv=FALSE, Rowv=TRUE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="log2(FPKM+1)")
#Figure S3B: Heatmap of well-known synaptogenic genes
#Syn2=ENSMUSG00000009394; Syp=ENSMUSG00000031144; Syt1=ENSMUSG00000035864;
specific_genes <- c("ENSMUSG00000009394", "ENSMUSG00000031144", "ENSMUSG00000035864")
heatmap_genes <- subset(gene_mat_fpkm, rownames(gene_mat_fpkm) %in% specific_genes)
rownames(heatmap_genes) <- c("Syn2","Syp","Syt1")
heatmap_genes_m <- data.matrix(heatmap_genes[,c(5,6,4,3,1,2)])
heatmap_genes_m_log <-log2(heatmap_genes_m+1)
break_hm <-c(seq(-0,10,by=0.1))
color_hm <-colorRampPalette(c("grey","white","maroon"))(100)
heatmap.2(heatmap_genes_m_log, 
          Colv=FALSE, Rowv=FALSE, dendrogram='none', col=color_hm, breaks=break_hm,
          scale="none", trace="none", key=TRUE, key.title=NA, key.ylab=NA, density.info=c("none"),
          srtCol=90, cexCol=0.6, cexRow=0.6, main="log2(FPKM+1)")
#Figure S3C: Heatmap of differentially-expressed genes
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
mySigMat+
  theme(panel.border = element_rect(colour="black"),plot.title=element_blank())+
  scale_fill_gradient(low="white", high="#9b1f51")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))
### GENE ONTOLOGY (GO) ANALYSIS
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
#combine the interneuron target gene list & ensembl-entrez database (to convert to gene IDs for GO test)
df_eg_ensembl<-as.data.frame(org.Mm.egENSEMBL)
deg_ss_entrez<-merge(deg_ss_names,df_eg_ensembl,by.x="X",by.y="ensembl_id",all=FALSE)
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
#combine the interneuron target gene list & ensembl-entrez database (to convert to gene IDs for GO test)
df_eg_ensembl<-as.data.frame(org.Mm.egENSEMBL)
deg_pv_entrez<-merge(deg_pv_names,df_eg_ensembl,by.x="X",by.y="ensembl_id",all=FALSE)
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
#combine the interneuron target gene list & ensembl-entrez database (to convert to gene IDs for GO test)
df_eg_ensembl<-as.data.frame(org.Mm.egENSEMBL)
deg_ch_entrez<-merge(deg_ch_names,df_eg_ensembl,by.x="X",by.y="ensembl_id",all=FALSE)
##Gene universe
#subset genes expressed in interneurons at P10 (threshold FPKM>5)
gene_in_fpkm <-gene_mat_fpkm[c(2,3,6)]
gene_in_fpkm5 <-gene_in_fpkm[apply(gene_in_fpkm, MARGIN = 1, function(x) any(x >= 5)), ]
#call the "org.Mm.egENSEMBL" database
df_eg_ensembl<-as.data.frame(org.Mm.egENSEMBL)
#combine the universe gene list & ensembl-entrez database (to convert gene names)
gene_in_fpkm5_entrez<-merge(gene_in_fpkm5,df_eg_ensembl,by.x="row.names",by.y="ensembl_id",all=FALSE)
#Get EntrezGene names (entrez gene = EG)
entrez_universe<-as.character(unique(gene_in_fpkm5_entrez$gene_id))
##Gene target
#get EntrezGene names (entrez gene = EG) for each interneuron population
entrez_target<-as.character(unique(deg_ss_entrez$gene_id))
entrez_target<-as.character(unique(deg_pv_entrez$gene_id))
entrez_target<-as.character(unique(deg_ch_entrez$gene_id))
#top genes after entrez gene name conversion
entrez_target<-head(entrez_target,100)
##GO test
#gene ontology category
go_category="CC"
#p-value cut-off
cutoff_value=0.05
#create "paramater object"
params <-new("GOHyperGParams",
             geneIds=entrez_target,
             universeGeneIds=entrez_universe,
             annotation="org.Mm.eg.db",
             ontology=go_category,
             pvalueCutoff=cutoff_value,
             testDirection="over")
#run the GO hypergeometric test
hgover<-hyperGTest(params)
##Create summary of data
#dataframe with hypergeometric tests and pvalues for GO terms
df_go_hg<-summary(hgover)
#dataframe with genes in each category that are represented in the target lists
hgover_geneids<-geneIdsByCategory(hgover)
df_go_geneids<-data.frame(go_term=rep(names(hgover_geneids),lapply(hgover_geneids,length)),
                          gene_id=unlist(hgover_geneids))
##Graphs for GO analysis
#Figure S3D: Bar plots for the three interneuron populations at P10
#select GO terms (gene set size minimum<25, maximum>2500, & only top 20 terms ranked by pvalue)
df_go<-df_go_hg
df_go<-df_go[ -which(df_go$Size<25),]
df_go<-df_go[ -which(df_go$Size>2500),]
df_go<-head(df_go,20)
#calculate  -log10(pvalue) for GO terms
df_go_plot<-df_go
df_go_plot$logpval<- -log10(df_go_plot$Pvalue)
#reorder by -log10(pvalue)
df_go_plot <- df_go_plot[order(df_go_plot$logpval),] 
#create new combined name "Term + GO:number"
df_go_plot$Term_go <-paste0(df_go_plot$Term," ","(",df_go_plot$GOCCID,")")
df_go_plot$Term_go_ <- factor(df_go_plot$Term_go, levels=unique(df_go_plot$Term_go))
#ss cells
color_go<-c("white","#FF8C5A")
breaks_go<-c(0,10,20,30,40)
limits_go<-c(0,40)
scale_go<-c(1,2,3,4,5)
#pv cells
color_go<-c("white","#04A4DD")
breaks_go<-c(0,10,20,30)
limits_go<-c(0,30)
scale_go<-c(1,2,3,4,5)
#ch cells
color_go<-c("white","#E0141A")
breaks_go<-c(0,10,20,30,40)
limits_go<-c(0,40)
scale_go<-c(1,2,3,4,5)
#plot
ggplot()+
  geom_bar(data=df_go_plot, aes(x=Term_go_,y=Count,fill=logpval), size=1.0,width=0.7,alpha=1,stat="identity")+
  coord_flip()+
  scale_fill_gradientn(colours=color_go,limits=c(1,max(df_go_plot$logpval)),breaks=scale_go)+
  scale_y_continuous(expand=c(0,0),limits=limits_go,breaks=breaks_go)+
  theme_bw()+
  theme(axis.title.y=element_blank())+
  theme(legend.position = "bottom")+
  labs(x="GO process", y="Gene count")
####