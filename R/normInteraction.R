normInteraction <- function(data, genename, interactome){

    ##gene names of all tested genes
    rnInteractome <- row.names(interactome)
    genename2Interactome <- rnInteractome%in%genename

    if(sum(genename2Interactome) == 0) stop("no tested genes are listed in the interactome matrix")

    testedInteractome <- interactome[genename2Interactome, ]
    cntestedI <- colnames(testedInteractome)
    cnData <- colnames(data)
    
    noEmptyComplex <- cntestedI%in%cnData
    testedInteractome <- testedInteractome[, noEmptyComplex]
   
    ## Normalization by the number of possible interaction by complexes
    normFactor <- colSums(testedInteractome)
    U <- sweep(data, 2, normFactor, "/")
    normData <- sweep(U, 1, normFactor, "/")
    normData <- round(normData, 2)
    idx <- which(is.na(normData))
    normData[idx] <- 0
    return(normData)
}
