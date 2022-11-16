#simulate via fsc and run gl2stairway
#source("~/GENE/Coala/gl.msfs.R")
setwd("d:/bernd/r/GENE/")
source("./Coala/gl.msfs.R")
source("./staiways2/gl2stairwayplot2.r")
library(strataG)
library(dartR)
library(parallel)
library(future)
library(furrr)



setwd("d:/programms/fsc/")
demes <- fscSettingsDemes(fscDeme(deme.size = 200, sample.size = 50, growth=0.05))
genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = 10, mut.rate = 1e-6), num.chrom = 100000)
events = fscSettingsEvents(
  #fscEvent(event.time = 20, source = 0, sink = 0, prop.migrants = 0, new.size=10),
  fscEvent(event.time = 40, source = 0, sink = 0, prop.migrants = 0, new.growth=0)
  
  )


p <- fscWrite(demes = demes, genetics = genetics, label = "shrink333", events = events)

prun <- fscRun(p, all.sites = F, dna.to.snp =TRUE, num.sims = 1, no.arl.output = F, sfs.type = "maf", exec = "fsc2709")
snp.df <- fscReadArp(prun, coded.snps = TRUE, marker = "snp")
gl <- new("genlight", snp.df[,-c(1:2)], ploidy=2, pop=snp.df[,2], ind.names=snp.df[,1])

#can be done to have only the first SNP on a sequence
#filter only first snp on sequence
#ll <- nchar(gl@loc.names)
#nn <- substr(locNames(gl),ll-2,ll)
#index <- (nn=="SNP" | nn=="_L1")
#glf <- gl[, index]
#table(as.matrix(glf))
#glf <- gl.compliance.check(glf)

#stairways has to have no D0!!!!!
#sfss <-gl.msfs(jerra2, minbinsize = 1) #now in the function via minbinsize
sfss <-gl.msfs(gl, minbinsize = 0)
barplot(sfss)
#setwd("~/GENE/staiways2/")

setwd("d:/programms/Stairway2/")
#Sys.setenv(DISPLAY="localhost:10.0")
gl2stairwayplot2(gl,simfolder = "simjerra", verbose = T,outpath = "d:/programms/Stairway2/", mu = 1e-6, gentime = 1,run=TRUE, nreps = 20, parallel=8, L=1e6, minbinsize =1, cleanup = FALSE, pct_training = 0.9)

