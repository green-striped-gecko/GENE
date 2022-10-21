#' Creates a site frequency spectrum based on a dartR or genlight object
#' 
#' @param x dartR/genlight object
#' 
#' @return returns a site frequency spectrum. If the genlight object consists of several populations the site frequency spectrum for each population is returned [=a multidimensional site frequency spectrum]. To get a single sfs for this case, you need to redefine the population via: pop(gl)<- rep("A", nInd(gl)) to a single population.
#' 
#' 


gl.msfs<- function(x) {
  pp <- seppop(x)
  mi <- max(sapply(pp,nInd))
  mm <- sapply(pp, function(x) table((0.5-abs(0.5-gl.alf(x)[,1]))*nInd(x)*2), simplify = FALSE)
  # cut off zero if necessary....
  mm2 <- sapply(mm, function(x) x[names(x)!="0"], simplify = FALSE)
  #max number of individuals (fill with zeros)
  mm3 <- t(sapply(mm2, function(x) c(x, rep(0, max(0,mi-length(x))))))
  colnames(mm3)<- paste0("d",1:mi)
  return(mm3)
}
