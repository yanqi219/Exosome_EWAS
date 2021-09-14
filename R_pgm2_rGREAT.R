rm(list=ls())
options(stringsAsFactor=F)
setwd('G:/Dropbox/JUno_work/HD/BloodMethylationEnrollHD_2018/')
library(WGCNA)
library(rGREAT)
pos.txt='mypos.txt'
neg.txt='myneg.txt'
background.txt='mybackground.txt'
out.csv='myoutput.csv'
#
#
pos1=read.delim(pos.txt)
neg1=read.delim(neg.txt)
backgroud=read.delim(backgroud.txt)
#
input=vector(len=2,mode='list')
input[[1]]=pos1
input[[2]]=neg1
names(input)[1:2]=c('pos','neg')
#
#(1)Binomial test: I use 50kb for a buff of gene regulation region. The default is 1000kb
#
output12={}
for(i in 1:length(input)){
  output.all={}
  job = submitGreatJob(input[[i]], 
                       species               = "hg19",
                       includeCuratedRegDoms = TRUE,
                       rule                  = c("basalPlusExt"),
                       adv_upstream          = 5.0,
                       adv_downstream        = 1.0,
                       adv_span              = 50,
                       request_interval = 300,
                       version="3",
                       max_tries = 10)
  #
  ontology.all=availableOntologies(job)
  print(ontology.all)
  
  for(k in 1:length(ontology.all)){
    print(ontology.all[k])
    out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k], download_by = "tsv"),error=function(e){NULL})
    if(!is.null(out0.list)){
      db0.list=as.list(names(out0.list))
      output<-Map(cbind,Database=db0.list,out0.list)
      output<-do.call('rbind',output)
      output.all=rbind(output.all,output)
    }
  }
  table(output.all$Database)
  output.all=data.frame(class=names(input)[i],extend='50kb',output.all)
  output12=rbind(output12,output.all)
}#end i
print(tail(output12[,1:10]))
#}#end TISSUE.id
write.table(output12,out.csv,sep=',',row.names = F)
#
#(2)Hypergeometic genomic region based analysis
#
job = submitGreatJob(input[[i]], 
                     bg=backgroud,
                     species               = "hg19",
                     includeCuratedRegDoms = TRUE,
                     rule                  = c("basalPlusExt"),
                     adv_upstream          = 5.0,
                     adv_downstream        = 1.0,
                     adv_span              = 50,
                     request_interval = 300,
                     version="3",
                     max_tries = 10)
