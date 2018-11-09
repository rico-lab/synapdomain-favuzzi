### Favuzzi et al.: SUPPLEMENTARY FIGURE S2
#Packages
source("https://bioconductor.org/biocLite.R")
biocLite("cummeRbund")
library(cummeRbund)
library(ggplot2)
library(cowplot)
library(reshape2)
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
#Define colors
cols <- c("ch8" = "#F46277", "ch10" = "#E0141A",
          "pv5" = "#7FD8E2", "pv10" = "#059FDB",
          "ss5" = "#FFA992", "ss10" = "#FF8C5A",
          "inp0" = "#6EB348", "pyr" = "#714091", "plp" = "#878787")
#Figure S2B: FPKM dispersion box plots
b<-csBoxplot(genes(cuff),logMode=TRUE,pseudocount=1,
             replicates=FALSE)
b$data$condition = factor(b$data$condition, 
                          levels=c('ss5','ss10','pv5','pv10',
                                   'ch8','ch10','inp0','pyr','plp'))
b+
  theme_bw()+
  theme(axis.text.y=element_text(size=11),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=11))+
  scale_fill_manual(values = cols)
#Figure S2C: Count dispersion plots
disp<-dispersionPlot(genes(cuff))
disp$data$sample_name_ = factor(disp$data$sample_name, 
                                levels=c('ss5','ss10','inp0',
                                         'pv5','pv10','pyr','ch8','ch10','plp'))
disp+
  theme_bw()+
  facet_wrap( ~ sample_name_, ncol=3)+
  theme(axis.text.y=element_text(size=11),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=11))+
  scale_colour_manual(values = cols)
#Figure S2D: Density plots
dens<-csDensity(genes(cuff))
dens+
  theme_bw()+
  theme(axis.text.y=element_text(size=11),axis.text.x=element_text(size=11),plot.title=element_blank())+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)
#Figure S2E: CV2 plots
genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv+
  theme_bw()+
  theme(axis.text.y=element_text(size=11),axis.text.x=element_text(size=11),plot.title=element_blank())+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)
#Figure S2F: Clustering based on Jensen-Shannon distance for replicates
dend.rep <-csDendro(genes(cuff),replicates=T)
#Figure S2G: Principal component analysis
genes.PCA<-PCAplot(genes(cuff),x="PC2",y="PC3")   #scale=TRUE(default) #showPoints=TRUE(default)
genes.PCA$layers[[1]]<-geom_point(alpha=0.1,color="gray55")
genes.PCA$layers[[2]]<-geom_hline(yintercept = 0, linetype="dashed", color="gray25")
genes.PCA$layers[[3]]<-geom_vline(xintercept = 0, linetype="dashed", color="gray25")
genes.PCA+
  scale_colour_manual(values = cols)+
  theme_bw()
## Create expression matrices
#Matrix of Read counts
gene.counts<-countMatrix(genes(cuff))
gene.rep.matrix<-repCountMatrix(genes(cuff))
#Matrix of FPKM values
gene.matrix<-fpkmMatrix(genes(cuff))
gene.rep.matrix<-repFpkmMatrix(genes(cuff))
## Matrices of Differential Expression analysis
#all
gene.diff<-diffData(genes(cuff))
#ch10
gene.diff_ch8.ch10<-diffData(genes(cuff),"ch8","ch10")
gene.diff_ch10.pv10<-diffData(genes(cuff),"ch10","pv10")
gene.diff_ch10.ss10<-diffData(genes(cuff),"ch10","ss10")
gene.diff_pv5.ch10<-diffData(genes(cuff),"pv5","ch10")
gene.diff_ss5.ch10<-diffData(genes(cuff),"ss5","ch10")
#pv10
gene.diff_pv5.pv10<-diffData(genes(cuff),"pv5","pv10")
gene.diff_ch10.pv10<-diffData(genes(cuff),"ch10","pv10")
gene.diff_ss10.pv10<-diffData(genes(cuff),"ss10","pv10")
gene.diff_ch8.pv10<-diffData(genes(cuff),"ch8","pv10")
gene.diff_ss5.pv10<-diffData(genes(cuff),"ss5","pv10")
#ss10
gene.diff_ss5.ss10<-diffData(genes(cuff),"ss5","ss10")
gene.diff_ch10.ss10<-diffData(genes(cuff),"ch10","ss10")
gene.diff_ss10.pv10<-diffData(genes(cuff),"ss10","pv10")
gene.diff_ch8.ss10<-diffData(genes(cuff),"ch8","ss10")
gene.diff_pv5.ss10<-diffData(genes(cuff),"pv5","ss10")
#ch8
gene.diff_ch8.ch10<-diffData(genes(cuff),"ch8","ch10")
gene.diff_ch8.pv5<-diffData(genes(cuff),"ch8","pv5")
gene.diff_ch8.ss5<-diffData(genes(cuff),"ch8","ss5")
#pv5
gene.diff_pv5.pv10<-diffData(genes(cuff),"pv5","pv10")
gene.diff_ch8.pv5<-diffData(genes(cuff),"ch8","pv5")
gene.diff_ss5.pv5<-diffData(genes(cuff),"ss5","pv5")
#ss5
gene.diff_ss5.ss10<-diffData(genes(cuff),"ss5","ss10")
gene.diff_ch8.ss5<-diffData(genes(cuff),"ch8","ss5")
gene.diff_ss5.pv5<-diffData(genes(cuff),"ss5","pv5")
#inp0
gene.diff_inp0.ch8<-diffData(genes(cuff),"inp0","ch8")
gene.diff_inp0.pv5<-diffData(genes(cuff),"inp0","pv5")
gene.diff_inp0.ss5<-diffData(genes(cuff),"inp0","ss5")
gene.diff_inp0.ch10<-diffData(genes(cuff),"inp0","ch10")
gene.diff_inp0.pv10<-diffData(genes(cuff),"inp0","pv10")
gene.diff_inp0.ss10<-diffData(genes(cuff),"inp0","ss10")
#pyr
gene.diff_pyr.ch10<-diffData(genes(cuff),"pyr","ch10")
gene.diff_pyr.pv10<-diffData(genes(cuff),"pyr","pv10")
gene.diff_pyr.ss10<-diffData(genes(cuff),"pyr","ss10")
gene.diff_pyr.ch8<-diffData(genes(cuff),"pyr","ch8")
gene.diff_pyr.pv5<-diffData(genes(cuff),"pyr","pv5")
gene.diff_pyr.ss5<-diffData(genes(cuff),"pyr","ss5")
gene.diff_pyr.inp0<-diffData(genes(cuff),"pyr","inp0")
#plp
gene.diff_plp.ch10<-diffData(genes(cuff),"plp","ch10")
gene.diff_plp.pv10<-diffData(genes(cuff),"plp","pv10")
gene.diff_plp.ss10<-diffData(genes(cuff),"plp","ss10")
gene.diff_plp.pyr<-diffData(genes(cuff),"plp","pyr")
gene.diff_plp.ch8<-diffData(genes(cuff),"plp","ch8")
gene.diff_plp.pv5<-diffData(genes(cuff),"plp","pv5")
gene.diff_plp.ss5<-diffData(genes(cuff),"plp","ss5")
gene.diff_plp.inp0<-diffData(genes(cuff),"plp","inp0")
####