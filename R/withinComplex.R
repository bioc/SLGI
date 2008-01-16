withinComplex <- function(data, interactome){

    interactome <- as(interactome, "matrix")
    data <- createSquareMatrix(data)
    rnData <- rownames(data)
    cnData <- colnames(data)

    rnInteractome <- rownames(interactome)
    cnInteractome <- colnames(interactome)
    
     ##--Down size data matrix
    dataIntr <- rnData %in% rnInteractome
    if(sum(dataIntr)==0) stop("No commun denominator: no genes are listed in the interactome matrix")
    newDatar <- data[dataIntr,] 

    dataIntc <- cnData %in% rnInteractome
    newData <- newDatar[,dataIntc]
    newData <- newData[,rownames(newData)]

    ##--Down size interactome
    IntData <- rnInteractome%in%row.names(newData)
    SLinteractome <- interactome[IntData,]
    SLinteractome <- SLinteractome[,colSums(SLinteractome)!=0]
    #--Order SLinteractome
    SLinteractome <- SLinteractome[rownames(newData),]

    ##--Building SL by SL pairs sharing complexes matrix (co-membership in one Complex)
    BtB <- tcrossprod(SLinteractome,SLinteractome)
    coM<- newData*BtB
    
    return(coM)
    
}
