library(minfi) 
library(RPMM) 
library(lumi) 
library(wateRmelon)
library(IlluminaHumanMethylation450kmanifest) 
sheetf = "/u/flashscratch/j/joshuazh/NMR/LBC_meth/idats/SampleSheetLBC36.csv" 
targets = read.csv(sheetf) baseDir="/u/flashscratch/j/joshuazh/NMR/LBC_meth/idats/idats/" 
targets$Basename = paste0(baseDir, targets$Basename) #For the read.metharray.exp() function: the base directory is used if you set targets = NULL, 
#if recursive = T then the function searches the base directory and all subdirectories for .idat file pairs (green and red) 
#otherwise it looks in the provided vector at targets = <vector> for file paths. 
#Depending on the size of methylation array and how much RAM you have (such as if it's for 450K array), 
#I would suggest splitting it into ~200 .idat file pairs (green and red) per run or else you will run out of RAM if you have around 16gb or 32gb.
#If force = T it just means that it ignores the usual error given if the .idat files are not all of the same size
RGset = read.metharray.exp(base = baseDir, targets = targets$Basename, recursive=T, force=F) 
#I also find it helpful to save the files after every step (like the RGset file) because my session crashes a lot if I don't allocate enough RAM 
Mset <- preprocessNoob(RGset) print("created Methylset") 
#this step takes the longest which is why I print something out to see if it's done gset = mapToGenome(Mset) grset = ratioConvert(gset, what = "both") 
#str(grset) 
#I like to just double check this to make sure nothing went wrong beta = getBeta(grset) finalnorm = data.frame(ID = row.names(beta),beta) saveRDS(finalnorm, file = yourDirectoryHere)
