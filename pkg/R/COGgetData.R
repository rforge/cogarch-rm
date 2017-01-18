
COGgetData <- function(Ticker, from = NULL, to = NULL, src ="yahoo", type="Close", delta=NULL, ...){
  numbTicker <- length(Ticker)
  for(i in c(1:numbTicker)){
    if (is.null(from) && is.null(to)){
      getSymbols(Ticker[i], src =src, ...)
    } else {
      if (is.null(from) && !is.null(to)){
        getSymbols(Ticker[i], to=to, src = src, ...)
      }
      if (!is.null(from) && is.null(to)){
        getSymbols(Ticker[i], from=from, src = src, ...)
      }
      if(!is.null(from) && !is.null(to)){
        getSymbols(Ticker[i], from=from, to = to, src = src, ...)
      }
    }
    dummy <- eval(parse(text=Ticker[i]))
    Index <- paste0(Ticker[i],".",type)
    dummy2 <- na.approx(dummy[, Index])
    # pensare a quali metodi di approssimazione l'utente è autorizzato ad usare
    if (i==1){
      dataPrice <- dummy2
    } else {
      dataPrice <- merge(dataPrice, dummy2)
    }
  }

  if (is.null(delta)){
    data <- setData(dataPrice)
  } else {
    data <- setData(dataPrice, delta = delta)
  }
  return(data)
}
