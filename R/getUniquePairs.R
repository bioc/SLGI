#############################################################
##
## Finding tested pairs and interacting pairs
##
#############################################################

getPairs <- function(iMat) {
    fp <- row(iMat)
    sp <- col(iMat)
    
    pairs <- data.frame(query=rownames(iMat)[fp],
                        array=colnames(iMat)[sp],
                        interact=as.logical(iMat))
    row.names(pairs) <- c(1:nrow(pairs))
    return(pairs)
}

expendMat<- function(iMat,respV=character(0)){
    reported <- respV %in% colnames(iMat)
    lresp <- length(respV)
    testM <- cbind(iMat,matrix(0,nrow=dim(iMat)[1],ncol=sum(!reported)))
    colnames(testM) <- c(colnames(iMat),respV[!reported])
    rownames(testM) <- rownames(iMat)

    return(testM)
}

getTestedPairs <- function(iMat, respV){
    ##--only return the pairs that interact, ignor respV
    if(missing(respV))
    stop("respV is missing, with no default.
The names of the genes that were originally tested are required")
    
    testM <- expendMat(iMat, respV)
    
    
    ##--split testM
    idxc <- which(colnames(testM) %in% rownames(testM))
    idxr <- which(rownames(testM) %in% colnames(testM))
    l <- length(idxc)            ##-number of query genes that are also in the reporter gene list
    subtestM1 <- testM[idxr,idxc]
    orn <- order(rownames(subtestM1))
    ocn <- order(colnames(subtestM1))
    subtestM1 <- subtestM1[orn,ocn]
    
    ##--extract pair list from sub matrix 1
    pairL1 <- getPairs(subtestM1)


    ##--remove duplicate pairs
    delete <- rep(FALSE,times=l*l)
    found <- rep(0, times=l*l)
    asym <- 0
    ##--delete one pair in the duplicated pairs, conditioning on the S.L results
    for (r in 1:l) {                            ##-the upper left triangle
        for (c in r:l) {
            if (subtestM1[r,c] == subtestM1[c,r]){   ##-if symmetric, delect either one
                delete[(r-1)*l+c] <- TRUE
                if(subtestM1[r,c] == 1)
                  found[(c-1)*l+r] <- 2
            }
            else {
                asym <- asym+1
                ifelse (subtestM1[r,c] == 0,         ##-if asymmetric, delete the one with non-S.L result
                        delete[(c-1)*l+r] <- TRUE,
                        delete[(r-1)*l+c] <- TRUE)  ##-subtestM1[c,r] == 0
                ifelse (subtestM1[r,c] == 0,
                        found[(r-1)*l+c] <- 1,
                        found[(c-1)*l+r] <- 1)
            }
        }
    }

    pairL1 <- pairL1[!delete,]
    found <- found[!delete]
    
    if (l == dim(testM)[1]) {
        if (l == dim(testM)[2]) {   
            pairL <- pairL1
            pairL[,3] <- found
            pairL[,4] <- c(recip=TRUE)
        } else {                                  ##-( l == dim(testM)[1] and l < dim(testM)[2] )
            subtestM2 <- testM[,-idxc]
            if(is.matrix(subtestM2)==TRUE){          
                pairL2 <- getPairs(subtestM2)
                pairL <- rbind(pairL1, pairL2)
                pairL[,3] <- as.numeric(pairL[,3])
                pairL[1:nrow(pairL1),3] <- found
                pairL[,4] <- FALSE
                pairL[1:nrow(pairL1),4] <- TRUE
            }else{
              xx = colnames(testM)[!colnames(testM)%in%rownames(testM)]
              pairL <- data.frame("query"= xx,
                                  "array" = names(subtestM2),
                                  "interact" = subtestM2,
                                  "recip"= FALSE)
          }
        }
    }else{
        if (l == dim(testM)[2]) {                 ##-( l < dim(testM)[1] and l == dim(testM)[2] )
          subtestM2 <- testM[-idxr,]
          pairL2 <- getPairs(subtestM2)
          pairL <- rbind(pairL1, pairL2)
      }else {                                  ##-( l < dim(testM)[1] and l < dim(testM)[2] )
          subtestM2 <- testM[-idxr,][,idxc]
          subtestM3 <- testM[,-idxc]
          pairL2 <- getPairs(subtestM2)
          pairL3 <- getPairs(subtestM3)
          pairL <- rbind(pairL1, pairL2,pairL3)
          
      }
        pairL[,3] <- as.numeric(pairL[,3])
        pairL[1:nrow(pairL1),3] <- found
        pairL[,4] <- FALSE
        pairL[1:nrow(pairL1),4] <- TRUE  
    }
    colnames(pairL)[4] <- "recip"
    l <- dim(pairL)[1]
    rownames(pairL) <- c(1:l)
    
    return(pairL)
   
}

getUniquePairs <- function(iMat,respV=character(0),only=FALSE) {

  ##--only return the pairs that interact, ignor respV
  if (only) {
    respV=character(0)
  }

  if (length(respV) > 0) {
    ##--Expend the interaction matrix 
      testM <- expendMat(iMat, respV)
  } else {
    testM <- iMat
  }
  ##--split testM
  idxc <- which(colnames(testM) %in% rownames(testM))
  idxr <- which(rownames(testM) %in% colnames(testM))
  l <- length(idxc)            ##-number of query genes that are also in the reporter gene list
  subtestM1 <- testM[idxr,idxc]
  orn <- order(rownames(subtestM1))
  ocn <- order(colnames(subtestM1))
  subtestM1 <- subtestM1[orn,ocn]
  
  ##--extract pair list from sub matrix 1
  pairL1 <- getPairs(subtestM1)

  ##--remove duplicate pairs
  delete <- rep(FALSE,times=l*l)
 
  asym <- 0
  ##--delete one pair in the duplicated pairs, conditioning on the S.L results
  for (r in 1:l) {                            ##-the upper left triangle
    for (c in r:l) {
      if (subtestM1[r,c] == subtestM1[c,r])   ##-if symmetric, delect either one
        delete[(r-1)*l+c] <- TRUE
      else {
        asym <- asym+1
        ifelse ( subtestM1[r,c] == 0,         ##-if asymmetric, delete the one with non-S.L result
                delete[(c-1)*l+r] <- TRUE,
                delete[(r-1)*l+c] <- TRUE)    ##-subtestM1[c,r] == 0
      }
    }
  }
  ##-number of assymmetric elements
  if (l == dim(testM)[1]) {
      if (l == dim(testM)[2]) {                 ##-( l == dim(testM)[1] and l == dim(testM)[2] )
          pairL <- pairL1[!delete,]
      } else {                                  ##-( l == dim(testM)[1] and l < dim(testM)[2] )
          subtestM2 <- testM[,-idxc]
          if(is.matrix(subtestM2)==TRUE){          
              pairL2 <- getPairs(subtestM2)
              pairL <- rbind(pairL1[!delete,],pairL2)
          }else{
              xx = colnames(testM)[!colnames(testM)%in%rownames(testM)]
              pairL <- data.frame("query"= xx,
                                  "array" = names(subtestM2),
                                  "interact" = as.logical(subtestM2))
          }
      }
  } else { 
      if (l == dim(testM)[2]) {                 ##-( l < dim(testM)[1] and l == dim(testM)[2] )
          subtestM2 <- testM[-idxr,]
          pairL2 <- getPairs(subtestM2)
          pairL <- rbind(pairL1[!delete,],pairL2)
      } else {                                  ##-( l < dim(testM)[1] and l < dim(testM)[2] )
          subtestM2 <- testM[-idxr,][,idxc]
          subtestM3 <- testM[,-idxc]
          pairL2 <- getPairs(subtestM2)
          pairL3 <- getPairs(subtestM3)
          pairL <- rbind(pairL1[!delete,],pairL2,pairL3)
      }
  }
 
  l <- dim(pairL)[1]
  rownames(pairL) <- c(1:l)
  if (only) {
    return(pairL[pairL[,3],][,1:2])
  } else {
  return(pairL)
  }
}
