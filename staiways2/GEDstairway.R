library(dartR)
library(furrr)

### idea: 

setwd("d:/bernd/Honours/ryan/data/")

#install.packages("dartR")
library(dartR)

tympo <- gl.read.dart(filename = "Report_DTym20-5504_3_moreOrders_SNP_mapping_2.csv", ind.metafile = "GedDartR2020_4.csv")

#filter for CGED only
lineata <- gl.keep.pop(tympo, pop.list = c("Cookanalla", "Aerial Paddock (lower)", "Aerial Paddock (upper)", "Airport/Tidbinbilla", "Bonshaw", "Campbell Park", "Jerra East", "Jerra West", "Majura", "Poplars", "Queanbeyan Nature Reserve" ))



#setwd("d:/bernd/projects/anna")
#gl <- gl.read.dart("Report-DRef16-2230_AM.csv", covfilename = "foxcov_anna_3.csv")



#get rid of skulls
#glf <- gl[!(indNames(gl)=="Tas_skull" | indNames(gl)=="Tas_skull_b")]

#gl<- gl.filter.callrate(glf, threshold = 0.99)



#no missing values
gl<- gl.filter.callrate(lineata, threshold = 1)



source("d:/bernd/r/dartr/r/gl2stairwayplot2.r")
gl2stairwayplot2(gl,simfolder = "foxes", verbose = T,outpath = "d:/programms/Stairway2", mu = 1.2e-8, gentime = 1, outfile="fox.blueprint",run=TRUE, nreps = 200, parallel=7)
