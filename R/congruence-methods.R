##------------------------------------------------------------##
##  Count number of commun genetic interaction between a pair
## of query or target genes
##
##------------------------------------------------------------##

getSharedInteraction <- function(iMat, mode="query"){

    if(mode=="targets")
      iMat <- t(iMat)
    
    
    shared = iMat %*% t(iMat)
    diag(shared)=0
    ans = expand.grid(shared)

    pairNames = paste(rep(rownames(shared),each=ncol(shared)),
              rep(colnames(shared), times=nrow(shared)),sep="-")

    nn = strsplit(pairNames, "-")

    ##diagonal
    idx1 = sapply(nn, function(x) x[1]==x[2])

    ##duplicated names
    orderNames = lapply(nn, sort)
    idx2 = duplicated(orderNames)

    ##both
    idx = idx1 | idx2
    
    ans = as.matrix(ans[!idx, ])
    rownames(ans) = pairNames[!idx]
    colnames(ans) = "Shared"
    
    ans
    
}

##------------------------------------------------------------##
##  Congruence score based on Ping Ye et al.
##                        Molecular Systems Biology 1:2005.0026
##
##------------------------------------------------------------##

congruence <- function(iMat,
                       sharedInt,
                       mode="query",
                       universe,
                       padjust=FALSE){

    if(mode=="targets")
      iMat = t(iMat)

    ## identify query genes also used as targets
    rname=rownames(iMat)
    cname=colnames(iMat)
    
    rn= rname%in%cname
    cn= cname%in%rname

    iMat = iMat[rn,cn]

    ## Count interactions found for each query genes
    m = rowSums(iMat)

    ## Count pairs of query(targets) genes
    npairs = nrow(sharedInt)
    pairNames = strsplit(rownames(sharedInt), "-")
   
    ## create matrix
    ## |share int. G1-G2 | G1 interaction | G2 interaction
    totalInt = lapply(pairNames, function(x) c(m[x[1]], m[x[2]]))
    totalInt = mapply(cbind,  totalInt)
    testMat = cbind(sharedInt, t(totalInt))

    ##Hypergeometric test
    numW = testMat[,2]
    numB = universe - numW
    numWdrawn = testMat[,1]
    numDrawn = testMat[,3]
    
    pvals <- phyper(numWdrawn - 1, numW, numB,
                    numDrawn, lower.tail=FALSE)      
    
    ##Adjusting p-values
    if(padjust){
        res <- round(-log(pvals/ (universe^2), 10))
    }else{
        res <- round(-log(pvals, 10))
    }
    names(res) <- rownames(sharedInt)
    res
}  
    



