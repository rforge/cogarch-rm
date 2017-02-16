portfolioCOG <- function(yuima.cogarch, yuima.data, Inweight, delta, my.start,
                  Maturity, lambda, n.comp, nrip=1000,
                  Interval, V0=100, lev=0.95, Measure = "VaR",
                  method = "Nelder-Mead"){
  modCog11<-yuima.cogarch
  price<-yuima.data@original.data
  logPrice <- na.approx(log(price))
  S0 <- as.numeric(tail(price,n=1L))



  logret <- apply(logPrice, 2, diff)
  mu<-apply(logret,2,mean)
  Maturity -> timestep
  My.dim<-dim(logret)
  logret<-logret- matrix(mu, My.dim[1], My.dim[2], byrow=TRUE)
  dummIca<- as.matrix(apply(logret,2,cumsum))
  ResICA <- fastICA(dummIca,n.comp = n.comp, fun="exp") # We recover the independent
  # signals using the default values for the function fastICA
  # matplot(ResICA$S, type = "l")
  # matplot(ResICA$S%*%ResICA$A, type = "l")
  ## Estimates COGARCH(1,1) on each independen signals
  my.resICA <- list()

  # 1. Definition of COGARCH(1,1) in yuima

  #modCog11 <- setCogarch(p=1, q=1, measure = list(intensity="1", df=list("dnorm(z, 0, 1)")), measure.type="CP")

  # COGARCH11paramICA <- NULL
  #
  # My.dim <- dim(logRet)
  # fitGARCHICA <- list()
  # spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  #                    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
  #
  # # GARCH11param <- NULL
  # GARCH11paramICA <- matrix(NA, nrow= n.comp,ncol=3)
  # for(j in c(1:n.comp)){
  #   fitGARCHICA[[j]] <- ugarchfit(data = diff(ResICA$S[,j]), spec = spec)
  #   #GARCH11param <- rbind(GARCH11param, c(coef(fitGARCH[[j]])))
  #   if(!is.null(coef(fitGARCHICA[[j]])))
  #     GARCH11paramICA[j,] <- c(coef(fitGARCHICA[[j]]))
  # }
  #
  # colnames(GARCH11paramICA) <- c("Omega", "Alpha1", "Beta1")
  #
  # ValuePar<- ParGarToCog(GARCH11paramICA, delta, rownames =as.character(c(1:n.comp)))
  COGARCH11paramICA<-NULL
  for(j in c(1:n.comp)){
    my.data <- setData(as.matrix(ResICA$S[,j]), delta = delta)
    yuima11 <- setYuima(data = my.data, model = modCog11)
    my.resICA[[j]] <- qmle(yuima = yuima11, grideq=TRUE,
                           start = my.start,
                           Est.Incr="IncrPar", aggregation = FALSE, method =method)
    COGARCH11paramICA <- rbind(COGARCH11paramICA, coef(my.resICA[[j]]))
  }

  ## Next step is the construction of a three for a portfolio
  ## As first step, we consider an equally weighted portfolio

  #weight <- matrix(1/My.dim[2],My.dim[2],1)
  weight <- Inweight
  ## simulate sample paths of independent signals.
  comp<-list()
  for(i in c(1:n.comp)){
    comp[[i]]<-SimBoot(S0 = 1, my.resICA[[i]], mu=0, Time = timestep, nrip =nrip)$logprice
  }

  invW <- ResICA$A

  samp <- setSampling(Initial = 0, Terminal = timestep, n = timestep/delta)

  #PortfolioValue <- NULL
  AssetPar <- list()
  for(i in c(1:My.dim[2])){
    AssetPar[[i]] <- matrix(NA,nrow=length(samp@grid[[1]][-1]), ncol=dim(comp[[1]])[2])
  }
  for(i in c(1:dim(comp[[1]])[2])){
    matrSign <- NULL
    for(j in c(1:length(comp)))
      matrSign <- cbind(matrSign,comp[[j]][,i])
    SimulateNoise <-  matrSign  %*% invW
    Assetdum <- list()
    for(j in c(1:My.dim[2])){
      AssetPar[[j]][,i] <- S0[j]*exp(as.numeric(mu[j])*samp@grid[[1]][-1]+SimulateNoise[,j])
    }
  }

  ## Construction Tree

  ## In order to construct a tree, we use an algorithm composed by the following steps.

  # 1. We consider a partition of time composed by three dates: 30, 60, 90
  Asset_0 <- NULL
  Asset <- NULL
  for(i in c(1:length(AssetPar))){
    Asset_0 <- rbind(Asset_0,AssetPar[[i]][1,])
    Asset <- rbind(Asset,AssetPar[[i]][timestep/delta,])
  }
  ## Now we consider two different approaches.
  # The first approach is based on a myopic algorithm.
  # We construct our portfolio on the interval [0,30] and we repeat
  # rolling the strategy.

  # We construct a partition of the state of Nature using the first asset

  PartStat<- seq(min(Asset[1,]), max(Asset[1,]), length.out = Interval+1)
  # In conclusion S_{1,0} < S_{1,1} < ... <S_{i,1}<...<S_{N}

  # Now consider discrete random variables that approximate the behaviour of
  # Assets after 30 days. As reference asset we consider Asset 1.

  # Given the MC trajectories, we aggregate them according the following role
  # 1) In the node i the value of Asset_j  is determined as
  #\[
  #Asset_{j,i} = \inf\left\{Asset_{j} s.t Asset_{1,i}<=Asset_{1}< Asset_{1,i}\right\}
  #\]
  MyAssetPart <- list()
  Prob <- NULL
  cond <-list()
  for(i in c(1: Interval)){
    cond[[i]]<- Asset[1,]>=PartStat[i] & Asset[1,] < PartStat[i+1]
    Prob <- c(Prob, sum(cond[[i]])/dim(Asset)[2])
  }
  for(j in c(1:dim(Asset)[1])){
    dummy <- numeric(length=Interval)
    for(i in c(1:Interval)){
      dummy[i]<- mean(Asset[j,cond[[i]]])
    }
    MyAssetPart[[j]]<-as.numeric(dummy[Prob>0])
  }
  ProbTrue<-Prob[Prob>0]/sum(Prob[Prob>0])

  Word<-Prob[Prob>0]
  for(j in c(1:dim(Asset)[1])){
    Word<-rbind(Word,MyAssetPart[[j]])
  }

  colnames(Word) <- paste0("w",c(1: length(ProbTrue)))
  rownames(Word) <- c("prob", colnames(logret))

  S_T0<- Asset_0[,1]
  names(S_T0)<-colnames(logret)
  ## We can solve the optimization problem
  myExp<-function(param,my.env){
    if(is.matrix(param)){
      if(dim(param)[1]!=1){
        param <- t(param)
      }
    }else{
      param<- t(as.matrix(param))
    }
    param<- my.env$V0*param/my.env$S_T0
    sum((param%*%my.env$Word[2:dim(my.env$Word)[1],])*my.env$Word[1,])
  }

  my.env <- new.env()
  assign("Word", Word, envir=my.env)
  assign("V0", V0, envir=my.env)
  assign("S_T0",S_T0, envir=my.env)
  param <- rep(1/My.dim[2],My.dim[2])

  # Value <- myExp(param,my.env)
  assign("lev",lev,envir=my.env)

  myVaR<-function(param,my.env){
    if(is.matrix(param)){
      if(dim(param)[1]!=1){
        param <- t(param)
      }
    }else{
      param<- t(as.matrix(param))
    }
    param<- my.env$V0*param/my.env$S_T0
    ret <- param%*%my.env$Word[2:dim(my.env$Word)[1],]
    Port<- rbind(ret,my.env$Word[1,])
    myPort <- Port[,order(Port[1,])]
    myPort[2,] <- 1-cumsum(as.numeric(myPort[2,]))
    myPort[1,] <- as.numeric(my.env$V0 - myPort[1,])
    res<-approx(x = myPort[2,], y=myPort[1,],xout = my.env$lev)$y
    return(res)
  }

  myES<-function(param,my.env){
    if(is.matrix(param)){
      if(dim(param)[1]!=1){
        param <- t(param)
      }
    }else{
      param<- t(as.matrix(param))
    }
    param<- my.env$V0*param/my.env$S_T0
    ret <- param%*%my.env$Word[2:dim(my.env$Word)[1],]
    Port<- rbind(ret,my.env$Word[1,])
    myPort <- Port[,order(Port[1,])]
    myPort[2,] <- 1-cumsum(as.numeric(myPort[2,]))
    myPort[1,] <- as.numeric(my.env$V0 - myPort[1,])
    VaR<-approx(x = myPort[2,], y=myPort[1,],xout = my.env$lev)$y
    res <- mean(myPort[1,myPort[1,]>=VaR])
    return(res)
  }
  #VaR95<- myVaR(param,my.env)
  assign("lambda", lambda, envir = my.env)
  assign("Measure",Measure,envir = my.env)
  myObj <- function(param,my.env){
    param[1]<- 1-sum(param[-1])
    Value <- myExp(param,my.env)
    if(my.env$Measure=="VaR"){
      VaR95<- myVaR(param,my.env)
      obj <- -(Value-my.env$lambda*VaR95)
    }else{
      ES95 <-myES(param,my.env)
      obj <- -(Value-my.env$lambda*ES95)
    }
    my.env$realPar1 <- param[1]
    return(obj)
  }

  myOpt<- constrOptim(theta=rep(1/My.dim[2],My.dim[2]), f = myObj, g=NULL,
                      ui = rbind(c(0,rep(-1,My.dim[2]-1)),cbind(0,diag(rep(1,My.dim[2]-1)))),
                      ci=c(-0.98,rep(0.02,My.dim[2]-1)), my.env=my.env)

  myOpt$par[1] <- my.env$realPar1
  condSumWeigth <- sum(myOpt$par[1])
  ExpWeal <-myExp(myOpt$par,my.env)
  ExpVaR <- myVaR(myOpt$par,my.env)
  ExpES <- myES(myOpt$par,my.env)
  res <- list(pesi=myOpt$par,ValoreAttesodiPor =ExpWeal, VaR=ExpVaR, ES = ExpES, Word=Word, Measure=Measure)
  return(res)
}
