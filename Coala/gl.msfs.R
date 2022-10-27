#' Creates a site frequency spectrum based on a dartR or genlight object
#' 
#' @param x dartR/genlight object
#' 
#' @return returns a site frequency spectrum. If the genlight object consists of several populations the site frequency spectrum for each population is returned [=a multidimensional site frequency spectrum]. To get a single sfs for this case, you need to redefine the population via: pop(gl)<- rep("A", nInd(gl)) to a single population.
#' 
#' 


gl.msfs<- function(x, minbinsize=1) {
  
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
  #delete minbinsize
  if (minbinsize>1) sfsf <- sfsf[-c(1:(minbinsize-1))]
  
  sfsl[[i]]<- sfsf
  }
  ll <- sapply(sfsl,length)
  if (all(ll==max(ll))) sfsl <- do.call(rbind,sfsl)
  return(sfsl)
}
