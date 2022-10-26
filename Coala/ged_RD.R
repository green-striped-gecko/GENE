

library(coala)
library(ape)
library(abc)
library(dartR)


gl.msfs<- function(x) {
  
  if (sum(is.na(as.matrix(x)))>0) cat("Your data contains missing data, better to use gl.impute to fill those gaps meaningful!\n")
  pp <- seppop(x)
  ml <- nLoc(x)
  mix <- max(table(pop(x)))
  sfsl <- list()
  for (i in 1:length(pp))
  {
    #mm <- sapply(pp, function(x) table((0.5-abs(0.5-gl.alf(x)[,1]))*nInd(x)*2), simplify = FALSE)
    mi <- nInd(pp[[i]])
    cs <- colSums(as.matrix(pp[[i]]), na.rm=TRUE)
    sfs0 <- table(mi-(abs(mi-cs)))
    sfsf <- rep(0,mix+1)
    sfsf[as.numeric(names(sfs0))+1] <- sfs0
    #delete monomorphs
    sfsf <- sfsf[-1]
    names(sfsf)<- paste0("d",1:mix)
    sfsl[[i]]<- sfsf
  }
  ll <- sapply(sfsl,length)
  if (all(ll==max(ll))) sfsl <- do.call(rbind,sfsl)
  return(sfsl)
}

#load ged data
#takes a minute for the full data set, therefore we save afterwards
#ged <- gl.read.dart("./data/SNP.csv", ind.metafile = "./data/GED.csv")


#please do not push to github so ignore when saved
#write as RDS for next time
#saveRDS(ged, "./data/ged.rda")
#download the rds file from google drive and store in the data folder...
ged <- readRDS("ged.rda")

summary(ged)

class(ged) <- "genlight"
summary(ged)

table(pop(ged))


#filter the data
#for fun we use only JerraEast or West

jerra <- gl.keep.pop(ged, pop.list="Jerrabomberra East")


nInd(jerra)
nLoc(jerra)

jerra2 <- gl.filter.callrate(jerra, method="loc", threshold = 1)
nInd(jerra2)
nLoc(jerra2)

# filter for minor allele frequency
# jerra2 <- gl.filter.maf(jerra2, threshold = 0.01)
# nInd(jerra2)
# nLoc(jerra2)

#!no missing data when running sfs function
#run gl.impute if necessary
target.sfs <- gl.msfs(jerra2)
target.sfs <- as.vector(target.sfs)
plot(target.sfs, type = "h")

# example fit
mod <- coal_model(sample_size = nInd(jerra2), loci_number = nLoc(jerra2), loci_length = 1) +
  feat_mutation(0.2) +
  sumstat_sfs()

# simulate
system.time(
  mod.sim <- simulate(mod, nsim = 1)
)
mod.sim$sfs

lines(mod.sim$sfs, col = "red")

# with decline
mod <- coal_model(sample_size = nInd(jerra2) + 1, loci_number = nLoc(jerra2), loci_length = 1) +
  feat_mutation(0.04) +
  feat_growth(rate = -2.15, time = 0) +
  feat_growth(rate = 0, time = 1) +
  sumstat_sfs()

# simulate
system.time(
  sim_data <- simulate(mod, nsim = 30, cores = 30)
)
sim_param <- create_abc_param(sim_data, mod)
sim_sumstat <- create_abc_sumstat(sim_data, mod)

plot(target.sfs, type = "h")
for(i in 1:nrow(sim_sumstat)) {
  lines(sim_sumstat[i, ], col = "red")
}

#-----------------------------------------------------------------------------------------------
# ABC example 
target.sfs

# set up the model to simulate 
sim.mod <- coal_model(sample_size = nInd(jerra2) + 1, loci_number = nLoc(jerra2), loci_length = 1) +
  feat_mutation(par_prior("theta", runif(1, 0, 0.5))) +
  feat_growth(rate = par_prior("gr.rate", runif(1, -4, 1)), time = 0) +
  feat_growth(rate = 0, time = par_prior("time.decline", runif(1, 0, 2))) +
  sumstat_sfs()

# run simulations in parallel
# function to simulate and output
system.time(
  sim_data <- simulate(sim.mod, nsim = 12000, cores = 30)
)
sim_param <- create_abc_param(sim_data, sim.mod)
sim_sumstat <- create_abc_sumstat(sim_data, sim.mod)


# estimate the posterior distribution of theta
posterior <- abc(target.sfs, sim_param, sim_sumstat, 0.05, method = "rejection")
hist(posterior, breaks = 20)

save(posterior, file = "posterior.RData")

# the accepted fits
dim(posterior$unadj.values)

# Euclidian distance of the accepted fits
history(posterior$dist[posterior$region])

# the best fitting
best <- which(posterior$dist[posterior$region] == min(posterior$dist[posterior$region]))
best.par <- posterior$unadj.values[best, ]
best.par


# fitted values
mod <- coal_model(sample_size = nInd(jerra2) + 1, loci_number = nLoc(jerra2), loci_length = 1) +
  feat_mutation(best.par[1]) +
  feat_growth(rate = best.par[2], time = 0) +
  feat_growth(rate = 0, time = best.par[3]) +
  sumstat_sfs()

# simulate
system.time(
  sim_data <- simulate(mod, nsim = 30, cores = 30)
)
sim_param <- create_abc_param(sim_data, mod)
sim_sumstat <- create_abc_sumstat(sim_data, mod)

plot(target.sfs, type = "h")
lines(colMeans(sim_sumstat), col = "red")


