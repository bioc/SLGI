getSharedDomains <- function(geneNameV,env) {
  if(is.environment(env)){
    stop("Environment are deprecated from genome annotation. Use database annotation package")
  }
  gnV <- as.matrix(geneNameV)
  env <- as.list(env)
  l <- length(gnV)
  if ( l <= 0 ) {
    stop("Zero or negative length name vector")
  } else {
    domains <- env[gnV]
    i <- 1
    sharedDomains <- domains[[i]]
    i <- i+1
    while ((i <= l) && (length(sharedDomains) > 0) && (!any(is.na(sharedDomains)))) {
      sharedDomains <- intersect(domains[[i]],sharedDomains)
      i <- i+1
    }
  }
  return(sharedDomains)
}

sharedBy <- function(domainL) {
  lens <- sapply(domainL,length)
  delete <- rep(FALSE,length(lens))
  delete[which(lens == 0)] <- TRUE
  delete[sapply(domainL,function(x) any(is.na(x)))] <- TRUE

  dl <- domainL[!delete]
  return(reverseSplit(dl))
}

domainDist <- function(domainL) {
  return(table(sharedBy(domainL)))
}
