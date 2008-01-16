topInteraction <- function(data, top=10){

    ## If the significant complexes are the top 10% 
    if(sum(lower.tri(data))!= sum(upper.tri(data))) stop("input 'matrix' must be a square matrice")
    
    lowertriData <- data[lower.tri(data)]
    shared <- round(lowertriData[which(lowertriData>1)], 2)
    top <- quantile(shared, 1-(top/100))
    
    data[upper.tri(data, TRUE)] <- NA
    topshared <- which(data>top, TRUE)
    top <- cbind(row.names(data)[topshared[, 1]], row.names(data)[topshared[, 2]], data[topshared]) 
    
    return(top)
    
}
