### Favuzzi et al.: SUPPLEMENTARY FIGURES S4, S5, & S6
#Packages
library(ggplot2)
## Single cell RNA sequencing databases from the adult brain:
#Tasic et al., 2016, Nat. Neurosci.: Figure S4
#Paul et al., 2017, Cell: Figure S5
#Zeisel et al., 2015, Science: Figure S6
## Function used for violin plots:
ggplot()+
  geom_violin(data=dataframe_scrnaseq,aes(x=Cell_id,y=Expression_level,fill=Cell_id),scale="width")+
  facet_grid(Gene_id ~ ., scales="free_y", switch="y")+
  scale_y_continuous(position = "right", expand=c(0,0))+
  scale_x_discrete(position = "top")+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90), strip.text.y=element_text(angle=180),
        axis.title.x = element_blank(), axis.ticks.y=element_blank(),
        strip.background = element_blank())
####