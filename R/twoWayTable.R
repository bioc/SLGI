twoWayTable <- function(var1,var2idx) {
  total <- length(var1)
  intTtotal <- sum(var1)
  shareTtotal <- length(var2idx)
  bothT <- sum(var1[var2idx])

  cell22 <- bothT
  cell21 <- intTtotal-bothT
  cell12 <- shareTtotal-bothT
  cell11 <- total-cell22-cell21-cell12
  sampleOR <- (cell11*cell22)/(cell12*cell21)

  tbl <- matrix(c(cell11,cell12,cell21,cell22),nrow=2,
                dimnames=list(var1=c("F","T"),var2=c("F","T")))
  return(tbl)
}
