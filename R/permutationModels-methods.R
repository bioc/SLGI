## ===========================================================================
## plot
## ---------------------------------------------------------------------------
## Plot function for permutation model results 
## ---------------------------------------------------------------------------
##  plotModelSLGI",
##  ignature = signature(x="siResult", y="missing"),
##  efinition = function(x, main="",...){
##     ans <- compare(x)
##     plot(density(ans), main=main,...)
##  
##  

setMethod("plot",
          signature=signature(x="siResult", y="missing"),
          function(x, n=500, auto.key = TRUE,
                   ylab="Average interaction", xlab="Index",...){

              a <- which(x@Observed >0)
              b <- apply(x@Expected, 1, mean)

              set.seed(123)
              idx = sample(1:length(a), n)

              expected <- b[a][idx]
              observed <- x@Observed[a][idx]

              model <- data.frame("Expected"=expected,"Observed"=observed)
              xyplot(Expected + Observed~ 1:n, data=model, ylab=ylab, xlab=xlab,
                     auto.key = auto.key,...)
              
          }
          )
     
## ===========================================================================
## count proportion
## ---------------------------------------------------------------------------
## Plot function for  permutation model results
## ---------------------------------------------------------------------------
setMethod("compare",
          signature=signature(x="siResult"),
          function(x){
              obsData <- x@Observed[x@Observed!=0]
              simulData <- x@Expected[x@Observed!=0, ]
              
              nIter <- length(obsData)
              nPerm <- ncol(simulData)
              ans <- c()
              i = 0
              count = c()
              for(i in 1:nIter){
                  count = sum(obsData[i] > simulData[i,])
                  if(count > 0){
                      ans[i] = count/nPerm
                  } else {
                      ans[i] = 0
                  }
              }
              ans
          }
          )
          
## ==========================================================================
## show method for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature=signature("siResult"),
          definition=function(object) {
            obs <- length(object@Observed)
            perm <- ncol(object@Expected)
            msg <- paste("siResult object, results of a graph approach on synthetic
             genetic interaction data with ",
                         obs, " observations and ", perm,
                         " permutation applied \n",
                         sep = "")
            cat(msg)
            return(msg)
        })
