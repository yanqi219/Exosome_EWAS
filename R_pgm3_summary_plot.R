rm(list=ls())
options(stringsAsFactor=F)
setwd('G:/Dropbox/JUno_work/HD/BloodMethylationEnrollHD_2018/')
library(WGCNA)
require(openxlsx)
library(writexl)
library(gridExtra)
library(ggpubr)
library(grid)
#
input0=read.csv('mygreatoutput.csv')
out0.png='myplot_DATA.png'
#
#input0$class.old=input0$class
#input0$class=ifelse(input0$class.old=='pos','Hypermethylation','Hypomethylation')#I added a column "class" on the rGREAT output for pos/neg
#input=subset(input0,ObsGenes>=5)#I sometimes only use this QC after I browse the results 
input=subset(input,TotalGenes<=3000 & TotalGenes >=10)#the number of genes in annotation (ex:Go term)
dim(input)
input$BinomRank.old=input$BinomRank
input$BinomRank<-NULL
input$HyperRank.old=input$HyperRank
input$HyperRank<-NULL
#
input$DatabaseClass=paste(input$Database,input$class)
input=input[order(input$DatabaseClass,input$BinomP),]
input$BinomRank<-unlist(sapply(split(input,input$DatabaseClass),function(x){BinomRank=rank(x$BinomP,ties.method='min')}))
summary(input$BinomRank.old-input$BinomRank)
#
input=input[order(input$DatabaseClass,input$HyperP),]
input$HyperRank<-unlist(sapply(split(input,input$DatabaseClass),function(x){HyperRank=rank(x$HyperP,ties.method='min')}))
summary(input$HyperRank.old-input$HyperRank)
#
input$log10P= -log10(input$BinomP)
input2.all=subset(input,BinomRank<=10)

#
input2.all=input2.all[order(input2.all$BinomP),]
unique(input2.all$Database)
#
library(ggpubr)
library(dplyr)
#
id1=substr(input2.all$Database,1,2)=='GO'
id2=input2.all$Database%in%c("MSigDB Predicted Promoter Motifs",  "MSigDB miRNA Motifs")
id3=input2.all$Database%in%c("MSigDB Immunologic Signatures"   , "MSigDB Cancer Neighborhood" ,"MSigDB Oncogenic Signatures" ,
                             "MSigDB Perturbation")
id4=input2.all$Database%in%c("MGI Expression: Detected", "Mouse Phenotype" )
id5=input2.all$Database%in%c("PANTHER Pathway","BioCyc Pathway")
id6=input2.all$Database%in%c("Disease Ontology","Human Phenotype","HGNC Gene Families","InterPro","Placenta Disorders" )
input2.all$MegaDatabase=NA
input2.all$MegaDatabase[id1]='GO'
input2.all$MegaDatabase[id2]='MSigDB Motif'
input2.all$MegaDatabase[id3]='MsigDB'
input2.all$MegaDatabase[id4]='Mouse'
input2.all$MegaDatabase[id5]='PANTHER& Metabolic pathway'
input2.all$MegaDatabase[id6]='Human Phenotype etc.'
table(input2.all$MegaDatabase,input2.all$Database)
table(input2.all$Database[is.na(input2.all$MegaDatabase)])

#
mega=unique(input2.all$MegaDatabase)
p1.list=list()
for(k in 1:length(mega)){
input2=subset(input2.all,MegaDatabase==mega[k])
out.png=gsub(out0.png,pattern='DATA',rep=mega[[k]])
input2$Desc=substr(as.character(input2$Desc),1,50)

input2$Desc=factor(as.character(input2$Desc),levels=unique(input2$Desc[order(input2$BinomP)]),ordered=T)
p1<-ggplot(input2,aes(x=log10P,y=Desc,size=TotalGenes,colour=-log10(BinomP)))+
  geom_point()+scale_color_gradient(low='blue',high='red')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=20),
        axis.text.y = element_text(size=12),
        axis.title=   element_text(size=20),
        plot.margin = margin(0.5, 1, 0, 1, "cm")
  )+  ylab(label = mega[k])

p1<-p1+theme(strip.background=element_rect(fill=c('lightblue')),
          strip.text.x   =element_text(color='black',size=20,face='bold'))+
          theme(strip.text.x = element_text(margin = margin(1.5,0,1.5,0, "cm")))
p1<-p1+facet_grid(.~class)
#
p1.list[[k]]=p1
png(out.png,width=16,height=12,units='in',res=300)
p <- grid.arrange( arrangeGrob(p1.list[[k]]))
dev.off()
}


         


