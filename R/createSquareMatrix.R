createSquareMatrix <- function(data){

    ## Add mismatching columns to rows
    onlyCol <- colnames(data)[!(colnames(data)%in%rownames(data))]
    addCol2Row <- matrix(0, nrow=length(onlyCol), ncol=ncol(data))
    row.names(addCol2Row) <- onlyCol
    dataR <- rbind(data, addCol2Row)

    ## Add mismatching rows to columns 
    onlyRow <- row.names(data)[!row.names(data)%in%colnames(data)]
    addRow2Col <- matrix(0, nrow=nrow(dataR), ncol=length(onlyRow))
    colnames(addRow2Col) <- onlyRow
    dataRC <- cbind(dataR, addRow2Col)

    ##order row and column
    index<- order(row.names(dataRC))
    dataRC <-dataRC[index, ]
    SquareMatrix <-dataRC[, row.names(dataRC)]    

    return(SquareMatrix)
}
