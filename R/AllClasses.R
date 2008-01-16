## ===========================================================================
## siResult
## ---------------------------------------------------------------------------
## A container for the results after applying a permutation model on
## synthetic genetic interaction
## ---------------------------------------------------------------------------
setClass("siResult",
         representation(Observed="numeric", Expected="ANY"),       
         prototype=list(Observed=numeric(0), Expected=matrix(0)))
