##-------------------------------------------------------
##
## Evaluate protein comembership within cellular organizational unit
##--------------------------------------------------------
byComplex <-  function(bpL, interactome) {
    if(!is(bpL,"list"))
      stop("bpL must be a list")
    if(!is(interactome,"matrix") )
      stop("interactome must a matrix")
    
    rn <- row.names(interactome)
    bNames <-  names(bpL)
    ans <- rep(0, ncol(interactome))
    for(i in 1:ncol(interactome) ) {
        nIn <-  0
        protInC <- rn[interactome[, i]>0]
        ##bait in Complex i
        wB <- bNames[bNames %in% protInC]
        if(length(wB) == 0 )
          warning("Complex", i, "has a problem (no bait)")
        else {
            ##for each bait
            for(j in wB)
              nIn <-  nIn + length(intersect(bpL[[j]], protInC))
        }
        ans[i] <-  nIn
    }
    ans
    names(ans) <- colnames(interactome)
    return(ans)
}
##-------------------------------------------------------
##
## Compute shared interaction
##--------------------------------------------------------
sharedInt <- function(pairL, interactome, threshold=0){

    pairSL <- pairL[pairL[,3],]
    pairSL <- pairSL[as.character(pairSL[,2])%in%rownames(interactome),]
    
    pairSL <- pairSL[as.character(pairSL[,1])%in%rownames(interactome),]
   
    SLlist <- split(as.character(pairSL[,2]),as.character(pairSL[,1]))
    tested <- unique(c(as.character(pairSL[,1]), as.character(pairSL[,2])))
   
    rn <- rownames(interactome)
    bNames <-  length(SLlist)

   
    ans <- vector(mode="list")
  
    interactomeR <- interactome[intersect(rn, tested), ]
   
    for(i in 1:bNames){
        x <- SLlist[i]
        s <- unlist(x)
        ns <- length(s)
        
        temp <- vector(mode="list", length=ns)
        nomen <- vector(mode="list")
        for(j in 1:ns){
           
            xx <- colSums(interactomeR[c(names(x), s[j]),])
            temp[[j]] <- xx[xx >= threshold] 
            nomen[[j]] <- paste(names(x),s[j], sep="--")
        }
        names(temp) <- nomen
        ans <- c(ans, temp)
        
    }
    
    ans
}  

##-------------------------------------------------------
##
## Compute shared interaction: who's where
## an attempt to answer pleiotropy
##--------------------------------------------------------
sharedIntMat <-function(pairL, interactome){

    xx <- cbind(pairL[,1:2],as.logical(pairL[,3]))
    sh <- sharedInt(xx, interactome,threshold=1)

    lapply(strsplit(names(sh),"--"), function(x) { 
        comp=names(pairL[[paste(x[1],x[2],sep="--")]]); 
        interactome[x,comp]})
    
}

##-------------------------------------------------------
##
## Compute within and between interactions
##--------------------------------------------------------

## interactome: complexData, keggPathways
## intM: interaction matrix or incidence matrix of any genetic interactions we have data on
## example: compM <- compCGraph(Atong, arrayGenes, ScISIverified)
getInteraction <- function(iMat, universe, interactome) {

    interactome <- as(interactome, "matrix")
    ## reduce the interactome to the observed genetic interactions
    
    allSLNames <- union(row.names(iMat), colnames(iMat))
    interactomeSub <- interactome[intersect(allSLNames, rownames(interactome)), ]

    cS <- colSums(interactomeSub)
    if(length(cS>0) != ncol(interactomeSub))
    interactomeSub <- interactomeSub[, cS>0]
    
    ## prepare results
    nComplex <- ncol(interactomeSub)
    bwMat <-  matrix(0, nc=nComplex, nr=nComplex)
    row.names(bwMat) <-  colnames(bwMat) <- colnames(interactomeSub)

    ## count all genetic interactions recorded
    cN = colnames(iMat)
    rN = rownames(iMat)
    for(i in 1:nrow(iMat) ) {
        g1 = rN[i]
        g2 = cN[iMat[i,]>0]
        g2 = setdiff(g2, g1)
        if( length(g2) == 0 ) next
        for( x in g2 ) 
            bwMat = bwMat + outer(interactomeSub[g1,],
    interactomeSub[x,]) + outer(interactomeSub[x,], interactomeSub[g1,])
    }
    ##return(list(ansM = ansM, CDs = interactomeSub))

    ## fix count for double count (record) when baits are also prey and they both found each other
    ## query(bait) genes also prey
    aAndq = intersect(rownames(iMat), universe)
    inA = intersect(aAndq, colnames(iMat))

    inAandC = intersect(inA, rownames(interactomeSub))
    
    ## prepare results
    fixupM = matrix(0, nc=nComplex, nr=nComplex)

    ## count genetic interactions counted twice
    for( i in 1:(length(inAandC)-1) ) {
        g1 = inAandC[i]
        for( j in (i+1):length(inAandC) ) {
            g2 = inAandC[j]
            int = min(iMat[g1, g2], iMat[g2, g1])
            if( int > 0 ){
                fixupM = fixupM +  outer(interactomeSub[g1, ], interactomeSub[g2, ]) +
                  outer(interactomeSub[g2, ], interactomeSub[g1, ]) 
            }
        }
    }

    ## fix diagonal as everything was counted twice
    diag(bwMat) <- diag(bwMat)/2 
    diag(fixupM) <- diag(fixupM)/2 
    ansM2 = bwMat - fixupM

    return(list(bwMat = ansM2, CDs = interactomeSub))
}

##----------------------------------------------------------------##
##                                                                ##
## Genetic interactions to complexes as they might not only be SL ##
##                                                                ##
##----------------------------------------------------------------##
gi2Interactome <- function(iMat, interactome, threshold=0) {
     baitS <- rownames(iMat) %in% rownames(interactome)
     preyS <- colnames(iMat) %in% rownames(interactome)

     
     GIandCD <- iMat[baitS, preyS]
     rS <- rowSums(GIandCD)
     cS <- colSums(GIandCD)
     GIandCD[rS > threshold, cS > threshold]
 }


##----------------------------------------------------------------##
##                                                                ##
## Summaryze the interactions, displaying the GO annotation       ##
##                                                                ##
##----------------------------------------------------------------##
iSummary <- function(iMat, n=10, reverse=FALSE){

    iMat[upper.tri(iMat)] = 0
    wH = which(iMat> 0)
    nM = nrow(iMat)
    col = ((wH - 1) %/% nM) + 1
    
    row = ((wH-1) %% nM) + 1
    
    ans=rep(NA, length(wH))
    for(i in 1:length(wH)) ans[i] = iMat[row[i], col[i]]

    
    complexP = vector("list", length=length(ans))
    for(i in 1:length(wH)) complexP[[i]] = c(rownames(iMat)[row[i]],
                                     colnames(iMat)[col[i]])
    ## for annotation
    thecomp = unlist(complexP)
    ##GO
    thego = grep("^GO", thecomp)
    if(length(thego)>0){
        if(require("GO.db"))
        gocomplex = unique(thecomp[thego])
        xx = as.list(GOTERM)
        annot =  xx[gocomplex]
        goTerm = sapply(annot,function(x) if(!is.null(x)){Term(x)}else{NA})
        names(goTerm) = gocomplex
    }
    ##MIPS
    themips = grep("^MIPS", thecomp)
    if(length(themips)>0){
        mipscomplex = unique(thecomp[themips])
        mips = getMipsInfo()
        mipsTerm = mips[mipscomplex]
        mipsTerm = sapply(mipsTerm, function(x) attr(x, "desc"))
    }
    ##KEGG and EBI to work on    
        
    names(complexP) = ans
    largeS = which(ans > n)
    ansL = ans[largeS]

  
      for(i in seq(along=largeS)) {
          cat(sprintf("---------Count: %2d -----------\n", ansL[i]))
          for(j in complexP[[largeS[i]]]) {
              cat(j)
              if(regexpr("^GO", j)>0)
                  annot = goTerm[j]
              if(regexpr("^MIPS", j)>0)
                  annot = mipsTerm[j]
              
              if(!is.null(annot) & !is.na(annot)){
                  cat(" ", annot)         
              }else {
                  cat(" Not found (possibly deprecated)")
              }
              cat("\n") 
          }
      }
    
    if(reverse==FALSE){    
        res <- complexP[largeS]
    } else{
        temp = sapply(complexP[largeS], function(x) paste(x[1], x[2], sep="-"))
        res = as.numeric(names(temp))
        names(res) = temp
    }
    res
}

##----------------------------------------------------------------##
##                                                                ##
## Summaryze the interactions, displaying the GO annotation       ##
##                                                                ##
##----------------------------------------------------------------##
test2Interact <- function(iMat, tMat, interactome){
    
    iMat[upper.tri(iMat)] = 0
    rN <- rownames(iMat)
    cN <- colnames(iMat)
    
    wH = which(iMat> 0)
    nM = nrow(iMat)
    col = ((wH - 1) %/% nM) + 1
    row = ((wH-1) %% nM) + 1
    
    ans=rep(NA, length(wH))
    ansT=rep(NA, length(wH))
    C1=rep(NA, length(wH))
    C2=rep(NA, length(wH))
    
    interactome <- as(interactome, "matrix")
    sizeC <- colSums(interactome)
       
    for(i in 1:length(wH)) {
        ans[i] = iMat[row[i], col[i]]
        C1[i] = rN[row[i]]
        C2[i] = cN[col[i]]
        ansT[i] = tMat[C1[i], C2[i]]
       
    }
   
    complexP <- data.frame("unit1"=C1, "unit2"=C2, "tested"=ansT, "interact"=ans,
                           "sizeC1"= sizeC[C1], "sizeC2"= sizeC[C2])
    rownames(complexP) <- paste(C1, C2, sep="--")
    complexP
    
}


##----------------------------------------------------------------##
##                                                                ##
## get the genes that were not tested  between two complexes      ##
##                                                                ##
##----------------------------------------------------------------##
getNonTested <- function(iMat, pairs, genename, universe, interactome){
    
    combinaisons = vector(mode="list", )
    qname = rownames(iMat)
    tname = union(universe, qname)

    npairs = nrow(pairs)
    
    for(i in 1:npairs){

        pairInt = interactome[, pairs[i, ]]
        pm = rowSums(pairInt)>0

        pairInt =  pairInt[pm, ]
        gnPairInt = rownames(pairInt)
    
        toTest = gnPairInt[!gnPairInt %in% tname]
        
        totest = totest[!totest %in% genename]
        pair = pair[totest, ]
    
        id1 = which(pair[, 1]>0)
        id2 = which(pair[, 2]>0)
    
        combinaisons[[i]]=paste(rep(rownames(pair)[id1],
                      each=length(id2)), rep(rownames(pair)[id2],
                        times=length(id1)), sep="--")
    }
    names(combinaisons) = paste(pairs[, 1], pairs[, 2], sep="--")
    combinaisons
}
