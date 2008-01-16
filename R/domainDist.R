domainDist <- function(domainL) {
  return(table(sharedBy(domainL)))
}
