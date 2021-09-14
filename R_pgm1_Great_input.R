#
#(1) prepare input files for hg19 assembly
#
#
rm(list=ls())
options(stringsAsFactor=F)
setwd('H:/JUno_work2/HD/BloodMethylationEnrollHD_2018/')
library(nlme)
library(WGCNA)
library(qqman)
library(calibrate)
library(R.utils)

library(WGCNA)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)

new0='level3/RandomSlope/EnrollHD2020/Metal/DNAmlongitudinal/output/summary/Metal_pgm3_combine_TRAIT_progression_two_enroll_last_reg_model4a_TEST_EnrollHD2020_final.gz'
pos.txt='level3/RandomSlope/EnrollHD2020/GREAT/TRAIT_progression_pos_top1000'
neg.txt='level3/RandomSlope/EnrollHD2020/GREAT/TRAIT_progression_neg_top1000'
#
X.test=c("motscore")
k=1
new.file=gsub(new0,pattern='TRAIT',rep=X.test[k])
new=read.delim(gzfile(new.file))
new=new[order(new$Meta.P),]
dim(new)
head(new)
new$Gene=as.character(new$Gene)
new$CHR1=paste('chr',new$CHR,sep='')
new$bp1=new$bp+1
id=new$CHR1=='chr23'
new$CHR1[id]='chrX'
id=new$CHR1=='chr24'
new$CHR1[id]='chrY'

pos=subset(new,Meta.cor>0)
neg=subset(new,Meta.cor<0)
pos=pos[1:1000,]
neg=new[1:1000,]
#
head(subset(pos,select=c(CHR1,bp,bp1,CpG)))
#CHR1       bp      bp1        CpG
#9    chr5 36877216 36877217 cg16453206
#60  chr12 57985290 57985291 cg13672103
#72  chr16 85866348 85866349 cg00709979
#93   chr8 10697017 10697018 cg01619129
#95  chr19 36980503 36980504 cg21775899
#100  chr6 30656577 30656578 cg08743794

write.table(subset(pos,select=c(CHR1,bp,bp1,CpG)),pos.txt,sep=' ',row.names = F,quote=F,col.names=F)
write.table(subset(neg,select=c(CHR1,bp,bp1,CpG)),neg.txt,sep=' ',row.names = F,quote=F,col.names=F)
#
#prepare background: example 450k
#
#The backgroud file has exactly the same format that you simply list all the CpGs in the 450k. If your EWAS involves both 450k and EPIC then use the overlap
#CpGs to build the background for hypergeometric test
