# Function available
estCOGuniv <- function(yuima.cogarch, yuima.data, startvalue,
                       delta=1/252, expModel=TRUE, Trend="const",method){
  dataPriceRes <- yuima.data@original.data
  if(expModel){
    if(Trend=="const"){
      fundsCol <- colnames(dataPriceRes)

      logPrice <- na.approx(log(dataPriceRes))

      dataPriceLogRet <- apply(logPrice, 2, diff)
      meanlogRet <- apply(na.omit(diff(logPrice)), 2, mean)

      my.dim <- dim(diff(logPrice))
      logRet <- diff(logPrice) - matrix(meanlogRet, my.dim[1], my.dim[2], byrow=TRUE)
      detDummy <- rbind(logPrice[1, ], logRet)
      detlogPrice <- apply(detDummy, 2, cumsum)


      my.res <- list()

      # 1. Definition of COGARCH(1,1) in yuima

      Check<-list()
      a<-NULL
      COGARCH11param<-NULL

      for(j in c(1:my.dim[2])){
        my.data <- setData(as.matrix(na.omit(na.approx(detlogPrice[,j]))), delta = delta)
        yuima11 <- setYuima(data = my.data, model = modCog11)
        my.res[[j]] <- tryCatch(qmle(yuima = yuima11, grideq=TRUE,
                                     start = startvalue,
                                     Est.Incr="IncrPar", aggregation = FALSE, method = method ),
                                error = function(yuima11){NULL})
        if(!is.null(my.res[[j]])){
          Check[[j]]<-Diagnostic.Cogarch(my.res[[j]])
          cat("\n stat, pos", c(Check[[j]]$stationary,Check[[j]]$positivity))
          COGARCH11param <- rbind(COGARCH11param, coef(my.res[[j]]))
          a<-rbind(a,as.character(c(Check[[j]]$stationary,Check[[j]]$positivity)))
        }else{
          a<-rbind(a,c("FAILED","FAILED"))
        }
      }
      colnames(a)<-c("stationarity","positivity")
      cat("\n")
      print(a)
      return(list(EstRes = my.res, mu = meanlogRet,
                  Info = a, Diagn = Check ))
    }else{
      return(NULL)
    }

  }else{
    warning("only exp cogarch is implemented")
    return(NULL)
  }

}


SimBoot <- function(S0,my.res, mu, Time, nrip, my.seed=2){
  mod <- my.res@yuima@model
  Incr.L <- my.res@Incr.Lev
  param <- c(coef(my.res), y1 = Diagnostic.Cogarch(my.res)$meanStateVariable)
  set.seed(my.seed)
  delta<-my.res@yuima@sampling@delta
  Ngrid<-Time/delta
  samp <- setSampling(Initial = 0, Terminal = Time, n = Ngrid)
  AllSboot <- matrix(NA, Ngrid, nrip)
  for(i in c(1 : nrip)){
    pos <- as.integer(runif(Ngrid, min = 1, max = length(Incr.L)))
    trajboot <- simulate(mod, true.parameter = param,
                         increment.L = as.matrix(as.numeric(Incr.L[pos])), samp = samp)
    AllSboot[ , i] <- samp@grid[[1]][-1]*mu + as.numeric(trajboot@data@original.data[ ,1])
  }
  logprice = log(S0)+AllSboot
  price=S0*exp(AllSboot)
  rownames(logprice) <- samp@grid[[1]][-1]
  rownames(price) <- samp@grid[[1]][-1]
  return(list(logprice=logprice,price=price))
}

univariateRM<- function(level, S0, my.res, mu,
                        Time, nrip,my.seed = 2){
  LossL <- NULL
  VaRL <- NULL
  ESL <- NULL
  VaRasPerOfInitValueL <- NULL
  ESasPerOfInitValueL <- NULL


  AllSboot <- suppressWarnings(SimBoot(S0, my.res, mu,
                                       tail(Time,1), nrip,my.seed = my.seed)$logprice-log(S0))
  delta<-my.res@yuima@sampling@delta
  for(j in c(1:length(Time))){
    VaR <- NULL
    ES <- NULL
    VaRasPerOfInitValue <- NULL
    ESasPerOfInitValue <- NULL

    LossDum <- S0*(1-exp(AllSboot))
    Loss <- LossDum[Time[j]/delta,]
    VaR <- c(VaR, quantile(Loss, level))
    VaRasPerOfInitValue <- c(VaRasPerOfInitValue, VaR/S0)
    ES <- c(ES,mean(Loss[Loss>=VaR]))
    ESasPerOfInitValue <- c(ESasPerOfInitValue, ES/S0)

    VaRL[[j]] <- VaR
    ESL[[j]] <- ES
    VaRasPerOfInitValueL[[j]] <- VaRasPerOfInitValue
    ESasPerOfInitValueL[[j]] <- ESasPerOfInitValue
  }

  VaR <- unlist(VaRL)
  names(VaR)<-paste0("Time",c(1:length(Time)))
  ES <- unlist(ESL)
  names(ES)<-paste0("Time",c(1:length(Time)))
  perVaR <- unlist(VaRasPerOfInitValueL)
  names(perVaR)<-paste0("Time",c(1:length(Time)))
  perES <- unlist(ESasPerOfInitValueL)
  names(perES)<-paste0("Time",c(1:length(Time)))
  return(list(VaR=VaR,ES=ES,perVaR=perVaR, perES=perES))
}
