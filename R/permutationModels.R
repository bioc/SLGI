modelSLGI <- function(iMat, universe, interactome, type="intM", perm=50){
   
    ##Reduce the interactome matrix to the one tested in iMat and universe
    nameTest <- union(rownames(iMat), universe)
    nn <- intersect(rownames(interactome), nameTest)
    idx <-  which(colSums(interactome[nn, ])> 0)
    interactomeR <- interactome[, idx]

    ##Observed interactions
    SLa <-  gi2Interactome(iMat, interactomeR)

    obs <-  getInteraction(SLa, universe, interactomeR)$bwMat
    obsi <- as(obs, "vector")
    names(obsi) <- paste(rownames(obs), rep(colnames(obs), each = nrow(obs)), sep="-")

    ##fornat result
    Inti <- vector(mode="list", length=perm)

    ##Permute 1 and 0 in interaction matrix
    if(type == "intM"){
        rowInt <- nrow(SLa)
        colInt <- ncol(SLa)
        rname <- rownames(SLa)
        cname <- colnames(SLa)
        intV <- as(SLa, "vector")
        for(i in 1:perm) {
            u <-  sample(intV)
            intNew <- matrix(u, nrow=rowInt, dimnames=list(rname, cname))
            Inti[[i]] <-  getInteraction(intNew, universe, interactomeR)$bwMat
        }
    }
    ##Permute genes labels in interactome
    if(type == "interactome"){
        geneInteractome <- rownames(interactomeR)
        for(i in 1:perm) {
         rownames(interactomeR) <-  sample(geneInteractome)
         Inti[[i]] <- getInteraction(SLa,  universe, interactomeR)$bwMat
                
     }
    }
    simulated <- mapply(cbind, Inti)
    colnames(simulated) <-  paste("simulation", 1:perm, sep="")
    rownames(simulated) <- names(obsi)

    return(new("siResult", "Observed"=obsi, "Expected"=simulated))
}
