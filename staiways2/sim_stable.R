#simulate via fsc and run gl2stairway
source("~/GENE/Coala/gl.msfs.R")
library(strataG)
library(dartR)
library(parallel)
library(future)
library(furrr)
setwd("~/fsc/fsc27_linux64/")
demes <- fscSettingsDemes(fscDeme(deme.size = 1000, sample.size = 50))
genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = 10, mut.rate = 1e-7), num.chrom = 100000)
events = fscSettingsEvents(
  fscEvent(event.time = 50, source = 0, sink = 0, prop.migrants = 0, new.size = 5))


p <- fscWrite(demes = demes, genetics = genetics, label = "stable333", events = events)

p <- fscRun(p, all.sites = F, dna.to.snp =FALSE, num.sims = 1, no.arl.output = F)
snp.df <- fscReadArp(p, coded.snps = TRUE, marker = "snp")
gl <- new("genlight", snp.df[,-c(1:2)], ploidy=2, pop=snp.df[,2], ind.names=snp.df[,1])

#filter only first snp on sequence
ll <- nchar(gl@loc.names)
nn <- substr(locNames(gl),ll-2,ll)
index <- (nn=="SNP" | nn=="_L1")

glf <- gl[, index]
table(as.matrix(glf))

glf <- gl.compliance.check(glf)
sfs <-gl.msfs(glf)
barplot(sfs)

setwd("~/GENE/staiways2/")
source("gl2stairwayplot2.r")
Sys.setenv(DISPLAY="localhost:10.0")
gl2stairwayplot2(glf,simfolder = "stable_100", verbose = T,outpath = "~/GENE/staiways2", mu = 1e-7, gentime = 1, outfile="stable100.blueprint",run=TRUE, nreps = 200, parallel=30, L=1e6)



