COGleader<-function(Data, NumbClust, type="STS", TITLE = "",
                    freq=1/252){
  dataPriceAll<-Data
  AssetNames <- colnames(dataPriceAll@original.data)
  namesIndex <- AssetNames
  dataPrice <- as.zoo(as.matrix(dataPriceAll@original.data))
  colnames(dataPrice) <- c(1:dim(dataPrice)[2])
  N <- dim(dataPrice)[1]
  ndataPrice <- dim(dataPrice)[2]
  D2 <- matrix(0, ndataPrice, ndataPrice)
  DELTA <- deltat(dataPrice)
  if(type=="STS"){
    for(i in 1:(ndataPrice-1)){
      for(j in (i+1):ndataPrice){
        D2[i,j] <- sqrt(sum((diff(dataPrice[,i])/DELTA-diff(dataPrice[,j])/DELTA)^2, na.rm = TRUE))
        D2[j,i] <- D2[i,j]
        cat(j)
      }
    }
    D2n <- D2/max(D2)
    colnames(D2n) <- colnames(dataPrice)
    rownames(D2n) <- colnames(dataPrice)
  }else{
    warning("Alternatives methods will be implemented as soon as possible")
    return(NULL)
  }
  k <- NumbClust
  cl2 <- hclust(as.dist(D2n))

  plot(cl2,
       main=TITLE,
       xlab="", ylim=c(0,1.1),
       ylab = "",sub="",
       cex.main=0.8,
       cex=0.5, cex.axis=0.8)


  l <- rect.hclust(cl2, k=k)

  namesIndex2 <- character()
  for(i in c(1:k)){
    #  group[[i]] <- dataPrice[, l[[i]]]
    if(length(l[[i]]) == 1){
      #   restricData <- merge(restricData, group[[i]])
      namesIndex2 <- c(namesIndex2, namesIndex[l[[i]]])
      cat(i)
    } else {
      maxGroup <- dataPrice[, l[[i]]]
      myData <- setData(maxGroup, delta = freq)
      forl <- llag(myData)
      test <- NULL
      for(j in c(1:dim(forl)[1])){
        test <- c(test, sum(forl[j, ]>0))
        cat(j)
      }
      maxTest <- max(test)
      condTest <- test==maxTest
      #restricData <- merge(restricData,dataPrice[, l[[i]][condTest]])
      namesIndex2 <- c(namesIndex2, namesIndex[l[[i]][condTest]])
    }
  }
  dataPrice1 <- dataPrice
  colnames(dataPrice1) <- namesIndex
  dataPriceRes <- zoo(dataPrice1[, namesIndex2],order.by = time(Data@original.data))
  colnames(dataPriceRes)<-namesIndex2
  finalData <- setData(dataPriceRes)
  return(finalData)
}

