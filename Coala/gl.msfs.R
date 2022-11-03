#' Creates a site frequency spectrum based on a dartR or genlight object
#' 
#' @param x dartR/genlight object
#' @param minbinsize remove bins from the left of the sfs. For example to remove singletons (alleles only occurring once among all individuals) set minbinsize to 2. If set to zero, also monomorphic (d0) loci are returned.
#' @param folded if set to TRUE (default) a folded sfs (minor allele frequency sfs) is returned. If set to FALSE then an unfolded (derived allele frequency sfs) is returned. It is assumed that 0 is homozygote for the reference and 2 is homozygote for the derived allele. So you need to make sure your coding is correct. 
#' @param singlepop switch to force to create a one-dimensional sfs, even though the genlight/dartR object contains more than one population
#' 
#' @return returns a site frequency spectrum. If the dartR/genlight object consists of several populations the multidimensional site frequency spectrum for each population is returned [=a multidimensional site frequency spectrum]. Be aware the multidimensional spectrum works only for a limited number of population and individuals [if too high the table command used internally will through an error as the number of populations and individuals (and therefore dimensions) are too large]. To get a single sfs for a genlight/dartR object with multiple populations, you need to set singlepop to TRUE.
#' 


gl.msfs<- function(x, minbinsize=0, folded=TRUE, singlepop=FALSE) {
  
  if (sum(is.na(as.matrix(x)))>0) cat("Your data contains missing data, better to use gl.impute to fill those gaps meaningful!\n")
  #only a single population....
  if (nPop(x)==0) {
    cat("No population definition provided. I proceed, assuming your genlight/dartR object is a single population.\n")
    pop(x)<- rep("A", nInd(x))
    }
  if (nPop(x)==1 | singlepop==TRUE)
  {
    mi <- nInd(x)
    if (!folded) mi=2*mi  #double the number of slots...
    cs <- colSums(as.matrix(x), na.rm=TRUE)
    if (folded) sfs0 <- table(mi-(abs(mi-cs))) else sfs0<- table(cs)
    sfsf <- rep(0,mi+1)
    sfsf[as.numeric(names(sfs0))+1] <- sfs0
    names(sfsf)<- paste0("d",0:mi)
    #delete minbinsize
    if (minbinsize>0) sfsf <- sfsf[-c(1:(minbinsize))]
    sfs <- sfsf
  } else  #multidimensional
  {
    pp <- seppop(x)
    cs <-list()
    for (i in 1:length(pp))
    {
    mi <- nInd(pp[[i]])
    if (!folded) mi=2*mi  #double the number of slots...
  sfs0 <- colSums(as.matrix(pp[[i]]), na.rm=TRUE)
   if (folded) cs[[i]] <- mi-(abs(mi-sfs0)) else cs[[i]] <- sfs0 
    }
  
  #add zeros and make sure they are consistent  
    
  
  msfs0 <- do.call(table, cs)  
  aa <-array(0,dim=table(pop(x))*2+1)
  dimnames(aa) <- sapply(dim(aa), function(x) paste0("d",0:(x-1)), simplify = F)
  dn <- dimnames(msfs0)
  
  do <- lapply(dn, function(x) 1:length(x))
  do1 <- expand.grid(do)
  do2 <- apply(do1, 2, as.numeric)
  dn1 <-expand.grid(dn)
  dn2 <- apply(dn1,2, as.numeric)
  dn3 <- dn2+1

  for (i in 1:nrow(dn3))
  {
    aa[t(dn3[i,])] <- msfs0[t(do2[i,])] 
  }
  #delete minbinsize
  if (minbinsize>0) 
  {
    tt <- paste0("aa[",paste0(rep("-c(1:minbinsize)",length(dim(aa))), collapse = ","),"]")
    aa <- eval(parse(text=tt))
    
  }
  
 
  
  sfs <- aa
  }
  return(sfs)
}
