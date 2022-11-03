#simulate via fsc and run gl2stairway
source("~/GENE/Coala/gl.msfs.R")
library(strataG)
library(dartR)
library(parallel)
library(future)
library(furrr)
setwd("~/fsc/fsc27_linux64/")

mean.sim <- function(p) {
  ns <- p$run.params$num.sims       
  ni <- p$settings$demes$sample.size        
  out.mean <- matrix(nrow = ns, ncol = ni+1)
  for(i in 1:ns) {
    snp <- fscReadSFS(p, sim = i)
    out.mean[i, ] <- snp$sfs$marginal[[1]][1:(ni+1)]
  }
  return(colMeans(out.mean))
}


demes <- fscSettingsDemes(fscDeme(deme.size = 100, sample.size = 100))
genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = 10, mut.rate = 2.5e-8), num.chrom = 100000)
events = fscSettingsEvents(
  fscEvent(event.time = 20, source = 0, sink = 0, prop.migrants = 0, new.size=0.1),
  fscEvent(event.time = 30, source = 0, sink = 0, prop.migrants = 0, new.size=10)
  
  )


p <- fscWrite(demes = demes, genetics = genetics, label = "shrink333", events = events)

prun <- fscRun(p, all.sites = F, dna.to.snp =TRUE, num.sims = 1, no.arl.output = F, sfs.type = "maf")


#sfss <- mean.sim(prun) #do not forget to take away the zero bin!!
snp.df <- fscReadArp(prun, coded.snps = TRUE, marker = "snp")
gl <- new("genlight", snp.df[,-c(1:2)], ploidy=2, pop=snp.df[,2], ind.names=snp.df[,1])

#obs.sfs <- fscReadSFS(prun)
#barplot(obs.sfs$sfs$marginal[[1]][2:40])

#filter only first snp on sequence
#ll <- nchar(gl@loc.names)
#nn <- substr(locNames(gl),ll-2,ll)
#index <- (nn=="SNP" | nn=="_L1")

#glf <- gl[, index]
#table(as.matrix(glf))

#glf <- gl.compliance.check(glf)
#stairways has to have no D0!!!!!
#sfss <-gl.msfs(jerra2, minbinsize = 1)
sfss <-gl.msfs(gl, minbinsize = 0)
barplot(sfss)
setwd("~/GENE/staiways2/")
source("gl2stairwayplot2.r")
Sys.setenv(DISPLAY="localhost:10.0")
gl2stairwayplot2(gl,simfolder = "simjerra", verbose = T,outpath = "~/GENE/staiways2", mu = 2.5e-8, gentime = 1,run=TRUE, nreps = 200, parallel=30, L=1.8e9, minbinsize =1)

