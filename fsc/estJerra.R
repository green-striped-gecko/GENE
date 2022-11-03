
# function to get mean of num.sims
mean.sim <- function(p) {
  ns <- p$run.params$num.sims       
  ni <- p$settings$demes$sample.size        
  out.mean <- matrix(nrow = ns, ncol = ni)
  for(i in 1:ns) {
    snp <- fscReadSFS(p, sim = i)
    out.mean[i, ] <- snp$sfs$marginal[[1]][2:(ni+1)]
  }
  return(colMeans(out.mean))
}


# specify scaled mutation rate theta = 4*N*mu
N0 <- 70
mut.rate <- 2.5e-7
ss <- nInd(jerra2)
nloc <- 10000
num.sims <- 20

# fastsimcoal sfs
demes <- fscSettingsDemes(fscDeme(deme.size = N0, sample.size = ss))
events <- fscSettingsEvents(fscEvent(event.time = 17, new.growth = 0))
genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = 1, mut.rate = mut.rate), num.chrom = nloc)

p <- fscWrite(demes = demes, genetics = genetics, events = events, label = "ex2.snps.1k2", use.wd = T)
system.time(
  p <- fscRun(p, all.sites = FALSE, dna.to.snp = TRUE, sfs.type = "maf", num.sims = 10, num.cores = 1)
)

snp <- fscReadSFS(p)




sfs.fsc <- mean.sim(p = p)
sfs.fsc

#plot(target.sfs, type = "h", ylim = c(0, 450))
#lines(sfs.fsc, col = "red")

sum(sfs.fsc)




snp.df <- fscReadArp(p, coded.snps = TRUE, marker = "snp")
gl <- new("genlight", snp.df[,-c(1:2)], ploidy=2, pop=snp.df[,2], ind.names=snp.df[,1])
sfs.fscb <- gl.msfs(gl)


sum(sfs.fsc)
sum(sfs.fscb)

par(mfrow=c(2,1))
barplot(sfs.fsc)
barplot(sfs.fscb)

#demes <- fscSettingsDemes(fscDeme(deme.size = N0, sample.size = ss, growth = 0.5))
#events <- fscSettingsEvents(fscEvent(event.time = 17, new.growth = 0))
#genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = 1, mut.rate = mut.rate), num.chrom = nloc)

demes <- fscSettingsDemes(fscDeme(deme.size = N0, sample.size = ss, growth="GROWTH"))
events <- fscSettingsEvents(
  fscEvent(event.time = "TDEC", 0, 0, 0, new.growth = 0))

est <- fscSettingsEst(
  #fscEstParam("NCUR", is.int = TRUE, distr = "unif", min = 10, max = 1000),
  # default for is.int = TRUE and distr = "unif" 
   fscEstParam("TDEC", min = 1, max = 100),
   fscEstParam("GROWTH",is.int=FALSE,min=-0.4, max=1),
  obs.sfs = snp$sfs$marginal[[1]]
)
est

est.p <- fscWrite(
  demes = demes,
  events = events,
  genetics = fscSettingsGenetics(fscBlock_freq(mut.rate)),
  est = est,
  label = "estTDEC"
)
system.time(est.p <- fscRun(est.p, num.sims = 1e6, num.cores = 20))


param.est <- fscReadParamEst(est.p )
cat(readLines("/tmp/Rtmp5ftkV3/estTDEC/estTDEC.bestlhoods"))
#readLines(file.path(est.p$folders, est.p$))


