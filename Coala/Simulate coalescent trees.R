
library(coala)
library(ape)
library(abc)

# the coala package uses a 4*N scaling

#-------------------------------------------------------------------------------
# simulate some coalescence trees

# generate coalescence model for 25 diploid individuals
mod <- coal_model(sample_size = 25, loci_number = 1, ploidy = 2) +
  sumstat_trees()

mod

# simulate 4 trees
mod.sim <- simulate(mod, nsim = 4)

# plot trees
par(mfrow = c(2, 2))
for(i in 1:4) {
  tree <- read.tree(text = mod.sim[[i]]$trees[[1]])
  plot(tree, show.tip.label = F, direction = "downwards")
}

#-------------------------------------------------------------------------------

# specify scaled mutation rate theta = 4*N*mu

mod <- coal_model(sample_size = 25, loci_number = 1, ploidy = 2) +
  feat_mutation(5) +
  sumstat_seg_sites() +
  sumstat_sfs()

# simulate
mod.sim <- simulate(mod, nsim = 1)

# dimensions of the SNP matrix
dim(mod.sim$seg_sites[[1]])

# extract the SNP matrix and calculate the site frequency spectra
mat <- get_snps(mod.sim$seg_sites[[1]])
mat

table(colSums(mat))

# get the sfs output from coala
mod.sim$sfs
length(mod.sim$sfs)

#-----------------------------------------------------------------------------
# lots of simulations of sfs
mod <- coal_model(sample_size = 25, loci_number = 1, ploidy = 2) +
  feat_mutation(5) +
  sumstat_sfs()

# simulate
mod.sim <- simulate(mod, nsim = 1000)

# extract sfs for each simulation
a <- lapply(mod.sim, function(x) x$sfs)
# put into a matrix
b <- matrix(unlist(a), nrow = length(a), byrow = T)

# average
plot(colMeans(b), type = "h")

# expected sfs = theta/i = 4*N*mu/i

xx <- 1:50
exp <- 5/(xx)
lines(exp ~ xx, col = "red")

#-----------------------------------------------------------------------------
# ABC example

# simulate a sfs 
n.ind <- 10
n.loci <- 100
mutation.rate <- 5

mod <- coal_model(sample_size = n.ind, loci_number = n.loci) +
  feat_mutation(mutation.rate) +
  sumstat_sfs()

# target sfs
mod.sim <- simulate(mod, nsim = 1, seed = 123)
target.sfs <- mod.sim$sfs
target.sfs

# set up the model to simulate with unknown theta parameter
# specify uniform distribution for prior on theta
sim.mod <- coal_model(sample_size = n.ind, loci_number = n.loci) +
  feat_mutation(par_prior("theta", runif(1, 1, 10))) +
  sumstat_sfs()

sim_data <- simulate(sim.mod, nsim = 2000)

# Get the parameter values from the simulation
sim_param <- create_abc_param(sim_data, sim.mod)
head(sim_param)

# Get the summary statistics from the simulation
sim_sumstat <- create_abc_sumstat(sim_data, sim.mod)
head(sim_sumstat)

# estimate the posterior distribution of theta
posterior <- abc(target.sfs, sim_param, sim_sumstat, 0.05, method = "rejection")
hist(posterior, breaks = 20)



