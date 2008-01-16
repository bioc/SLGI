hyperG <- function(data, nbTested, universe){
    ## Here's a diagram based on using GO as the category:
    ##
    ##          SL    notSL
    ##          White   Black
    ## Tested    n11     n12
    ## not       n21     n22
    
    data <- data[data[, 1]!=0,]
     
    ## Num white drawn (n11)
    numWdrawn <- data[, 1]
    numW <- nbTested
    numB <- universe-nbTested
    numDrawn <-  data[, 2]
    
    n21 <- numW - numWdrawn
    n12 <- numDrawn - numWdrawn
    n22 <- numB - n12
    
    odds_ratio <-  (numWdrawn * n22) / (n12 * n21)

    expected <- (numWdrawn + n12) * (numWdrawn + n21)
    expected <- expected / (numWdrawn + n12 + n21 + n22)

    pvals <- vector(mode="numeric",length=nrow(data))
      
    pvals <- phyper(numWdrawn-1, numW, numB, numDrawn, lower.tail=FALSE, log.p = FALSE)
   
    data.frame(data, P=pvals, Odds=odds_ratio, Expected=expected)
}
