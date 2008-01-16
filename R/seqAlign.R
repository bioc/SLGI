callMatcher <- function(seqbank,pidx) {
  dput(seqbank[[ pidx[1] ]]$seq,file="seq1",control=NULL)
  dput(seqbank[[ pidx[2] ]]$seq,file="seq2",control=NULL)
  system('sed s/\\"//g seq1 > seq1.new')
  system('sed s/\\"//g seq2 > seq2.new')
  system("matcher seq1.new seq2.new out.matcher")
}

new.alignRes <- function() {
  result.structure <- list(names=c("",""),
                           results=c(Length=0,Identity=0,Similarity=0,Gaps=0,Score=0),
                           sequences=c("")
                           )
  return(result.structure)
}

skiplines <- function(nline,con) {
  k <- 1
  while (k<=nline) {
    l <- readLines(con,n=1,ok=FALSE)
    k <- k+1
  }
  return(k)
}

read6lines <- function(con) {
  l1 <- substring(readLines(con,n=1,ok=FALSE),8)
  l2 <- substring(readLines(con,n=1,ok=FALSE),8)
  l3 <- substring(readLines(con,n=1,ok=FALSE),8)
  l4 <- substring(readLines(con,n=1,ok=FALSE),8)
  l5 <- substring(readLines(con,n=1,ok=FALSE),8)
  l6 <- readLines(con,n=1,ok=FALSE)
  ans <- rbind(l1,l2,l3,l4,l5)
  rownames(ans) <- c()
  return(ans)
}

pasteAs <- function(as1,as2) {
  ans <- c()
  for (i in 1:5) {
    ans[i] <- paste(as1[i],as2[i],sep="")
  }
  return(ans)
}

readMatcherResult <- function(pairNameV) {
  res <- new.alignRes()
  res$names <- pairNameV
  
  atEnd <- FALSE
  con <- file("out.matcher","r")

  #--skip the 1st 16 lines
  k <- try(skiplines(16,con),silent=FALSE)
  if (inherits(k,"try-error")) {
    close(con)
    stop("damaged file,type 1\n")
  }
  
  for (s in names(res$results)) {
    l <- try(readLines(con,n=1,ok=FALSE),silent=FALSE)
    if (inherits(l,"try-error")) {
      close(con)
      stop("damaged file, type 2\n")
    }
    res$results[s] <- strsplit(l,split="[\t\\/]+")[[1]][3]
  }

  #--reading matched sequences
  #--skip next 4 lines
  k <- try(skiplines(4,con),silent=FALSE)
  if (inherits(k,"try-error")) {
    close(con)
    stop("damaged file, type 3\n")
  }
  
  as1 <- try(read6lines(con),silent=FALSE)
  if (inherits(as1,"try-error")) {
    close(con)
    stop("damaged file, type 4\n")
  }  
  while (!atEnd) {
    as2 <- try(read6lines(con),silent=TRUE)
    if (inherits(as2,"try-error")) {
      atEnd <- TRUE
      close(con)
      break
    }
    as1 <- pasteAs(as1,as2)
  }
  res$sequences <- as1
  return(res)
}


seqMatcherAlign <- function(pairNameV,BankIDV,seqBank) {
  l <- length(pairNameV)
  if (l != 2) {
    stop("Not the right length of input pair names.")
  }
  idx <- sapply(pairNameV,function(x) { match(x,BankIDV) } )
  callMatcher(seqBank,idx)
  alignRes <- readMatcherResult(pairNameV)
  return(alignRes)
}

getAlignStats <- function(alignRes) {
#  names <- array(alignRes$names,dim=c(1,length(alignRes$names)))
#  results <- array(alignRes$results,dim=c(1,length(alignRes$results)))
  stats <- data.frame(rbind(c(alignRes$names,alignRes$results)))
#  dim(stats) <- c(1,length(stats))
  return(stats)
}
