
  TFdata = read.delim("../extdata/binding_by_gene.tsv", skip=2,
     as.is=TRUE, head=FALSE)

  h1 = scan("../extdata/binding_by_gene.tsv", n=1, what="", sep="\n")
  h1 = strsplit(h1, "\t")[[1]]


  h2 = scan("../extdata/binding_by_gene.tsv", n=1, what="", sep="\n", skip=1)
  h2 = strsplit(h2, "\t")[[1]]

  pvs = which(h1 =="p-value")

  TFN = which(h2 != "")


  TFmat = as.matrix(TFdata[,c(pvs)])

  dimnames(TFmat) = list(TFdata[,1],  h2[TFN])

  save(TFmat, file="../../data/TFmat.rda")




