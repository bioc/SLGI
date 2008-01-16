getFASTAname <- function(Fobj) {
  return(strsplit(Fobj$desc,"[> ]")[[1]][2])
}
