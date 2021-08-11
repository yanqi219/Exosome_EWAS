rm(list=ls())
options(stringsAsFactors = F)
setwd('G:/Steve Horvath Lab Dropbox/Ake Lu/')
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)
#
info=read.csv('JamesClement_hUCP/ResultsClockClement.csv')
#info=read.csv('JUno_work/JamesClement_hUCP')

#
yvars    =c('AgeAccelerationResidual','AgeAccelerationResidualHannum','DNAmAgeSkinBloodClockAdjAge','AgeAccelPheno','AgeAccelGrim')
yvars.txt=c('AgeAccelHorvath','AgeAccelHannum','AgeAccelSkinClock'    ,'AgeAccelPhenoAge','AgeAccelGrimAge')
info=info[order(info$BH.ID),]
id1=which(info$Time=='BaselineDraw')
id2=which(info$Time=='SecondDraw')
info$Age.new=info$Age
info$Age.new[id2]=info$Age[id2]+7
head(subset(info,select=c(BH.ID,Age,Age.new,Time)))
info1=subset(info,Time=='BaselineDraw')
info2=subset(info,Time!='BaselineDraw')
summary(info1$Age)
summary(info2$Age.new-info1$Age)
#
p1.list = list()
out.all1={}
for(k in 1:length(yvars)){
info$y=info[,yvars[k]]
info1$y=info1[,yvars[k]]
info2$y=info2[,yvars[k]]
info2$y.correct=info2$y+(info1$Age-info2$Age)
p1=t.test(info2$y,info1$y,paired=T)$p.value
diff=round(mean(info2$y-info1$y),digits=2)
p.plot<-ggplot(data=info,aes(x=Age.new,y=y))+geom_point()+geom_smooth(method='lm',se=F)+
  aes(colour=factor(BH.ID))+xlab('Age (baseline) to Age+7 (0.2yr)')+ylab(yvars.txt[k])+
  ggtitle(paste(LETTERS[k],'     ',yvars.txt[k],'\n','drop=',diff,'p=',
                format(p1,scientific=T,digits=2)))+labs(col='BH.ID')
p1.list[[k]] = p.plot
out0=data.frame(yvars[k],diff,p1)
out.all1=rbind(out.all1,out0)
}
#names(out.all1)=c('DNAm','FollowUp-baseline','T','FollowUp-baseline with age correction','T')
#
png('JUno_work/JamesClement_hUCP/SpaghettiPlot_DNAmClocks.png',width=20,height=9,units='in',res=300)
p <- grid.arrange(              arrangeGrob(p1.list[[1]]),
                                arrangeGrob(p1.list[[2]]),
                                arrangeGrob(p1.list[[3]]),
                                arrangeGrob(p1.list[[4]]),
                                arrangeGrob(p1.list[[5]]),
                                ncol=3,nrow=2)
dev.off()

#
#pro
#DNAmADM	DNAmB2M	DNAmCystatinC	DNAmGDF15	DNAmLeptin	DNAmPACKYRS	DNAmPAI1	DNAmTIMP1

yvars    =c('DNAmADMAdjAge',	'DNAmB2MAdjAge',	'DNAmCystatinCAdjAge',	
            'DNAmGDF15AdjAge',	'DNAmLeptinAdjAge',	'DNAmPACKYRSAdjAge'	,
            'DNAmPAI1AdjAge',	'DNAmTIMP1AdjAge')
yvars.txt=yvars
#
p1.list = list()
np=length(yvars)
rm(p1,diff)
for(k in 1:np){
  info$y=info[,yvars[k]]
  info1$y=info1[,yvars[k]]
  info2$y=info2[,yvars[k]]
  #info2$y.correct=info2$y+(info1$Age-info2$Age)
  diff=round(mean(info2$y-info1$y),digits=2)
  p1=t.test(info2$y,info1$y,paired=T)$p.value
 
  if(k<np){
  p.plot<-ggplot(data=info,aes(x=Age.new,y=y))+geom_point()+geom_smooth(method='lm',se=F)+
    aes(colour=factor(BH.ID))+xlab('Age (baseline) to Age+7 (0.2yr)')+ylab(yvars.txt[k])+
    ggtitle(paste(LETTERS[k],'     ',yvars.txt[k],'\n','drop=',diff,'p=',format(p1,scientific=T,digits=1)))+
    labs(col='BH.ID')+ theme(legend.position = "none")
  }else{
    p.plot<-ggplot(data=info,aes(x=Age.new,y=y))+geom_point()+geom_smooth(method='lm',se=F)+
      aes(colour=factor(BH.ID))+xlab('Age (baseline) to Age+7 (0.2yr)')+ylab(yvars.txt[k])+
      ggtitle(paste(LETTERS[k],'     ',yvars.txt[k],'\n','drop=',diff,'p=',
                    format(p1,scientific=T,digits=1)))+labs(col='BH.ID')
  }
  p1.list[[k]] = p.plot
  out0=data.frame(yvars[k],diff,p1)
  out.all1=rbind(out.all1,out0)
  
  
}
#names(out.all1)=c('DNAm','FollowUp-baseline based on AgeAccel','T')

#
png('JUno_work/JamesClement_hUCP/SpaghettiPlot_DNAmProtiens.png',width=18,height=9,units='in',res=300)
p <- grid.arrange(              arrangeGrob(p1.list[[1]]),
                                arrangeGrob(p1.list[[2]]),
                                arrangeGrob(p1.list[[3]]),
                                arrangeGrob(p1.list[[4]]),
                                arrangeGrob(p1.list[[5]]),
                                arrangeGrob(p1.list[[6]]),
                                arrangeGrob(p1.list[[7]]),
                                arrangeGrob(p1.list[[8]]),
                                ncol=3,nrow=3)
dev.off()
#
yvars.txt<-yvars<-c('CD8.naive','CD4.naive','CD8pCD28nCD45RAn','PlasmaBlast','CD4T','NK','Mono','Gran')
#
p1.list = list()
np=length(yvars)
rm(p1,diff)
for(k in 1:np){
  info$y=info[,yvars[k]]
  info1$y=info1[,yvars[k]]
  info2$y=info2[,yvars[k]]
  #info2$y.correct=info2$y+(info1$Age-info2$Age)
  diff=round(mean(info2$y-info1$y),digits=2)
  p1=t.test(info2$y,info1$y,paired=T)$p.value
  
  if(k<np){
    p.plot<-ggplot(data=info,aes(x=Age.new,y=y))+geom_point()+geom_smooth(method='lm',se=F)+
      aes(colour=factor(BH.ID))+xlab('Age (baseline) to Age+7 (0.2yr)')+ylab(yvars.txt[k])+
      ggtitle(paste(LETTERS[k],'     ',yvars.txt[k],'\n','FU-baseline=',diff,'p=',format(p1,scientific=T,digits=1)))+
      labs(col='BH.ID')+ theme(legend.position = "none")
  }else{
    p.plot<-ggplot(data=info,aes(x=Age.new,y=y))+geom_point()+geom_smooth(method='lm',se=F)+
      aes(colour=factor(BH.ID))+xlab('Age (baseline) to Age+7 (0.2yr)')+ylab(yvars.txt[k])+
      ggtitle(paste(LETTERS[k],'     ',yvars.txt[k],'\n','FU-baseline=',diff,'p=',
                    format(p1,scientific=T,digits=1)))+labs(col='BH.ID')
  }
  p1.list[[k]] = p.plot
  out0=data.frame(yvars[k],diff,p1)
  out.all1=rbind(out.all1,out0)
}
#
png('JUno_work/JamesClement_hUCP/SpaghettiPlot_CellCounts.png',width=18,height=9,units='in',res=300)
p <- grid.arrange(grobs=p1.list,ncol=3,nrow=3)
dev.off()

write.table(out.all1,'JUno_work/JamesClement_hUCP/paired  t-test of age acceleration.csv',sep=',',row.names = F,quote=F)