comemberIn <- function(iMat, interactome){

    interactome <- as(interactome, "matrix")
    rnData <- rownames(iMat)
    cnData <- colnames(iMat)
    if(sum(rnData %in% cnData) != length(rnData))
      stop("input 'data'should be a matrix within identical row and column")
    
    rnInteractome <- rownames(interactome)

    IntData <- rnInteractome %in% rnData
    if(sum(IntData)==0) stop("No commun denominator: no genes are comember of any of these complex")
    
    B <- interactome[IntData, ]
    B <- B[, colSums(B)!=0]
    B <- B[rnData, ]
     
    data2complex <- which(iMat > 0, arr.ind=TRUE)
    result <- vector(mode="list", length=nrow(data2complex))

    for(k in 1:nrow(data2complex)){
        result[[k]] <- names(which(B[data2complex[k,1], ] * B[data2complex[k,2], ] > 0))
    }
    names(result) <- paste(rnData[data2complex[, 1]], "-", cnData[data2complex[, 2]], sep="")
    
    ans <- unlist(result)
    return(ans)
    
}
