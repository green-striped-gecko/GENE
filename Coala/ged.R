
library(coala)
library(ape)
library(abc)
library(dartR)
#new multi-D-sfs function
source("./Coala/gl.msfs.R")


#load ged data
#takes a minute for the full data set, therefore we save afterwards
#ged <- gl.read.dart("./data/SNP.csv", ind.metafile = "./data/GED.csv")


#please do not push to github so ignore when saved
#write as RDS for next time
#saveRDS(ged, "./data/ged.rda")
#download the rds file from google drive and store in the data folder...
ged <- readRDS("./data/ged.rda")
#check ge

table(pop(ged))


#filter the data
#for funwe use only JerraEast or West

jerra <- gl.keep.pop(ged, pop.list="Jerrabomberra East")


nInd(jerra)
nLoc(jerra)

jerra2 <- gl.filter.callrate(jerra, method="loc", threshold = 1)
nInd(jerra2)
nLoc(jerra2)

#filter for minor allele frequency

#jerra3 <- gl.filter.maf(jerra2, threshold = 0.05)

#down to 4727 loci, and 49 individuals
#!no missing data when running sfs function
#run gl.impute if necessary
sfs <- gl.msfs(jerra2)
barplot(sfs)






