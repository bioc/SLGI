getSharedDomains <- function(geneNameV,env) {
  gnV <- as.matrix(geneNameV)
  l <- length(gnV)
  if ( l <= 0 ) {
    stop("Zero or negative length name vector")
  } else {
    domains <- mget(gnV,env,ifnotfound=NA)
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
