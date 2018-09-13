######       ch 7 icin ilk denemeler
### normal simjoint data turetip, long measurement mean and varyansina bakiyrz. jointplot cizidiryorz robustnessi nasil diye
library(joineR)
library(MASS)
library(dplyr)
`simjoint2` <-
  function(n=500,model=c("intslope","int","quad"),sepassoc=FALSE,ntms=5,b1=c(1,1,1,1),b2=c(1,1),gamma=c(1,0.1),sigu,vare=0.01,theta0=-3,theta1=1,censoring=T,censlam=exp(-3),truncation=F,trunctime=max(ntms),gridstep=0.01)
  {
    
    model <- match.arg(model)
    if(model!="intslope"&&model!="int"&&model!="quad")
    {stop(paste("Unknown model", model))}
    ran=2
    if(model=="int"){ran=1}
    if(model=="quad"){ran=3}
    lat=ran
    if(!sepassoc){lat=1}
    if(length(gamma)!=lat){
      cat("Warning: Number of association parameters do not match model choice\n")}
    gamma=rep(gamma,length=ran)
    if(missing(sigu)){sigu<-diag(ran)}
    if(length(sigu)!=ran^2){
      #cat("Warning: Dimension of covariance matrix does not match chosen model\n")
      if(length(sigu)>ran^2){sigu=sigu[1:ran,1:ran]}
      else{sigu=diag(ran)*sigu[1]}}
    if(model=="int"){if(sigu<0){stop("Variance must be positive")}}
    else{if(!isSymmetric(sigu)){
      stop("Covariance matrix is not symmetric")}
      if(any(eigen(sigu)$values<0)||(det(sigu)<=0)){
        stop("Covariance matrix must be positive semi-definite")}}
    
    "getD1"<-function(q,arg)
    {
      D<-matrix(0,q,length(arg))
      for(i in 1:q){D[i,]=arg^(i-1)}
      D
    }
    
    "simdat"<-function(n,model,sepassoc,ntms,b1,b2,gamma,sigu,vare,theta0,theta1,censoring,censlam,truncation,trunctime,gridstep)
    {
      ctsx<-rnorm(n)
      binx<-runif(n,-sqrt(3),sqrt(3))
      X2<-cbind(ctsx,binx)
      id<-1:n
      idl<-rep(id,each=ntms)
      ctsxl<-rep(ctsx,each=ntms)
      binxl<-rep(binx,each=ntms)
      time<-rep(0:(ntms-1),length=n*ntms)
      X1<-cbind(intercept=1,ctsxl,binxl,ltime=time)
      U<-mvrnorm(n,mu=rep(0,ran),Sigma=sigu)
      Ul<-U[rep(1:n,each=ntms),]
      D<-getD1(ran,time)
      DU<-t(D)*Ul
      Y<-X1%*%b1+rowSums(DU)+sqrt(vare)*rnorm(n*ntms)
      u0<-U[,1]
      if(model=="intslope"){u1=U[,2]}
      else{u1=rep(0,n)}
      b2x<-X2%*%b2
      cens<-rep(1,n)
      if(!sepassoc){gamma=rep(gamma[1],ran)}
      if(model!="quad"){
        if(model=="int"){gamma=c(gamma[1],0)}
        uu<-runif(n)
        if(model=="int")
        {
          survtime<--log(uu)/exp(theta0+b2x+gamma[1]*u0)
        }
        else{
          ii<-((theta1+gamma[2]*u1)<0)&(uu<exp(exp(theta0+b2x+gamma[1]*u0)/(theta1+gamma[2]*u1)))
          survtime<-rep(0,n)
          survtime[ii]<-Inf
          survtime[!ii]<-log(1-(theta1+gamma[2]*u1[!ii])*log(uu[!ii])/exp(theta0+b2x[!ii]+gamma[1]*u0[!ii]))/(theta1+gamma[2]*u1[!ii])
        }}
      else{
        tau<-trunctime
        tgrid<-seq(runif(1,0,gridstep),tau,gridstep)
        lam0<-exp(theta0+theta1*tgrid)
        hazt<-gridstep*exp(b2x)%*%lam0
        gD2<-gamma*getD1(ran,tgrid)
        hmat<-exp(U%*%gD2)*hazt
        uu<-matrix(runif(length(hmat)),n,length(tgrid))
        tmat<-matrix(tgrid,n,length(tgrid),byrow=T)
        tmat[hmat<uu]<-tau
        survtime<-apply(tmat,1,min)
        cens[survtime==tau]=0}
      
      if(censoring){censtime=-log(runif(n))/censlam}
      else{censtime=rep(Inf,n)}
      if(model!="quad"){if(truncation){censtime=pmin(censtime,trunctime)}}
      ii<-censtime<survtime
      survtime[ii]<-censtime[ii]
      cens[ii]<-0
      ls<-rep(survtime,each=ntms)
      Y<-Y[ls>time]
      X1<-X1[ls>time,]
      idl<-idl[ls>time]
      time<-time[ls>time]
      cat(100*sum(cens)/n,"% experienced event\n")
      event.rate <- 100*sum(cens)/n
      list(longdat=data.frame(id=idl,Y,time,X1),survdat=data.frame(id,survtime,cens,X2), eventrate=event.rate)
    }
    
    sim<-simdat(n,model,sepassoc,ntms,b1,b2,gamma,sigu,vare,theta0,theta1,censoring,censlam,truncation,trunctime,gridstep)
    list(longitudinal=sim$longdat,survival=sim$survdat, event.rate=sim$eventrate)
  }

jointfit <- function(n=500, gamma=1, model = "intslope", ntms=10, censlam = exp(-3), theta0 = -3, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6), sigu=matrix(c(0.9, 0.1, 0.1, 0.3),2,2), vare=1){
  mydata <- simjoint2(n=n, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 = theta0, theta1 = theta1, b1=b1, b2=b2, sigu = sigu, vare=vare)
  longdat <- mydata$longitudinal
  survdat<- mydata$survival
  longdat2 <- subset(longdat, select=c("id","Y","time"))
  survdat2 <- subset(survdat, select=c("id","survtime","cens"))
  covdat <- subset(survdat, select=c("id", "ctsx", "binx"))
  jointsimdat <-  jointdata(longitudinal=longdat2, baseline=covdat, survival=survdat2, id.col="id", time.col="time")
  gnc <- longdat %>% group_by(id) %>% summarise(ltime = max(ltime))
  gnc$ltime <- gnc$ltime+1
  gnc <- data.frame(gnc)
  eventrate <- mydata$event.rate
  fit<- joint(data=jointsimdat, long.formula=Y~ctsx+binx+time,
              surv.formula=Surv(survtime,cens)~(ctsx+binx), model="intslope")
  jointlongcoef<- fit$coefficients$fixed$longitudinal$b1
  jointsurvcoef1<-  fit$coefficients$fixed$survival[[1]]
  jointsurvcoef2<-  fit$coefficients$fixed$survival[[2]]
  latentassoc1<- fit$coefficients$latent[[1]]
  jointsigma.u1 <- fit$sigma.u[1]
  jointsigma.u2 <- fit$sigma.u[2]
  jointsigma.u4 <- fit$sigma.u[4]
  jointsigma.z <- fit$sigma.z[1]
  meandat=mean(gnc$ltime)
  vardat=var(gnc$ltime)
  coeff <- c(jointlongcoef, jointsurvcoef1, jointsurvcoef2,latentassoc1,jointsigma.u1,jointsigma.u2, jointsigma.u4,jointsigma.z, meandat, vardat, eventrate)
  return(coeff)
}



simrun <- function(n=500,nsim=3, gamma=1, model = "intslope", ntms=10, censlam = exp(-3), theta0 = -3, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6), sigu=matrix(c(0.9, 0.1, 0.1, 0.3),2,2), vare=1){
  result <- matrix(NA, nrow=nsim, ncol=14) #create a matrix to hold outcome
  for(i in 1:nsim) {
    result[i,] <- jointfit(n=n, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 =theta0, theta1 =theta1)
    print(i)
  }
  colnames(result) <- c("long.intercept","long.cont","long.bin","long.time","surv.cont","surv.bin",
                        "association", "sigma.u1","sigma.u2","sigma.u4","sigma.z", "meandat", "vardat", "eventrate")
  list(result=result, b1=b1, b2=b2, gamma=gamma, sigu=sigu, vare=vare )
}


ptm <- proc.time()
A1000.4 <- simrun(n=1000,nsim=1000,  gamma=-1, model = "intslope", ntms=20, censlam = exp(-4.5), theta0 = -7, theta1 = 0.2, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6),sigu=matrix(c(10,1,1,6), 2, 2))
A500.4 <- simrun(n=500,nsim=1000,  gamma=-1, model = "intslope", ntms=20, censlam = exp(-4.5), theta0 = -7, theta1 = 0.2, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6),sigu=matrix(c(10,1,1,6), 2, 2))
A250.4 <- simrun(n=250,nsim=1000,  gamma=-1, model = "intslope", ntms=20, censlam = exp(-4.5), theta0 = -7, theta1 = 0.2, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6),sigu=matrix(c(10,1,1,6), 2, 2))
stoptime<- proc.time() - ptm # Stop the clock
time.ch7 <- stoptime[3]/60 #in minutes.
print(time.ch7)



# true.values <- c(7, -1, 0.5, 0.1, 0.6, -0.6, 0.4, 0.7, 0.9, -0.04, 0.16, 0.01, NA, NA, NA)
# joint.mean <-apply(A250$result, 2, mean)
# joint.sd <- apply(A250$result, 2, sd)
# finalresultintslope<-data.frame(true.values, joint.mean,joint.sd, sep.mean,sep.sd)
# print(finalresultintslope)
# 
# 
# data1 <- jointfit(n=500, gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6))
# 
# jointplot(data1$jointsimdat, Y.col="Y",Cens.col = "cens",mean.profile = TRUE, split = TRUE,ylab="Y",gp1lab="Censored",gp2lab="Event")
# data1$meandat
# data1$vardat
# 
# data2 <- robustness(n=250, gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5)
# 
# jointplot(data2$jointsimdat, Y.col="Y",Cens.col = "cens",mean.profile = TRUE, split = TRUE,ylab="Y",gp1lab="Censored",gp2lab="Event")
# data2$meandat
# data2$vardat
# 
# data3 <- robustness(n=100, gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5)
# 
# jointplot(data3$jointsimdat, Y.col="Y",Cens.col = "cens",mean.profile = TRUE, split = TRUE,ylab="Y",gp1lab="Censored",gp2lab="Event")
# data3$meandat
# data3$vardat
# 
# 
# 
# robustness.deneme <- function(n=n, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 = theta0, theta1 = theta1){
#   mydata <- simjoint(n=n, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 = theta0, theta1 = theta1)
#   longdat <- mydata$longitudinal
#   survdat<- mydata$survival
#   n=dim(survdat)[1]
#   longdat2 <- subset(longdat, select=c("id","Y","time"))
#   survdat2 <- subset(survdat, select=c("id","survtime","cens"))
#   covdat <- subset(survdat, select=c("id", "ctsx", "binx"))
#   jointsimdat <-  jointdata(longitudinal=longdat2, baseline=covdat, survival=survdat2, id.col="id", time.col="time")
#   gnc <- longdat %>% group_by(id) %>% summarise(ltime = max(ltime))
#   gnc$ltime <- gnc$ltime+1
#   gnc <- data.frame(gnc)
#   #  fit<- joint(data=jointsimdat, long.formula=Y~ctsx+binx+time,
#   #            surv.formula=Surv(survtime,cens)~(ctsx+binx), model="intslope")
#   #  fitSE <- jointSE(fit, n.boot = 100)
#   list(meandat=mean(gnc$ltime), vardat=var(gnc$ltime), n=n, jointsimdat=jointsimdat)
# }
# 
# data3 <- robustness.deneme(n=500, gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5)
# 
# jointplot(data3$jointsimdat, Y.col="Y",Cens.col = "cens",mean.profile = TRUE, split = TRUE,ylab="Y",gp1lab="Censored",gp2lab="Event")
# data3$meandat
# data3$vardat
# n <- 250
# ll.est250 <- t(sapply(A250$result, function(m) m$fitSE$Estimate))
# ll.se250 <- t(sapply(A250$result, function(m) m$fitSE$SE))
# n <- 500
# ll.est500 <- t(sapply(A500$result, function(m) m$fitSE$Estimate))
# ll.se500 <- t(sapply(A500$result, function(m) m$fitSE$SE))
# n <- 1000
# ll.est1000 <- t(sapply(A1000$result, function(m) m$fitSE$Estimate))
# ll.se1000 <- t(sapply(A1000$result, function(m) m$fitSE$SE))
# 
# A <- list(b1 = c(1, 1, 1, 1), b2 = c(1, 1), gamma = 1, sigu=diag(diag(2)), vare = 0.01)
# nsim <- 100
# 
# 
# 
# 
# 
# 
# 
# res.mean <- ll.est1000
# res.sd <- ll.se1000
# colnames(res.mean) <- colnames(res.sd)<- c("l.int", "l.cont", "l.bin", "l.time", "s.cont", "s.bin", "gamma", "sigu1", "sigu2", "vare")
# 
# t05=qt(0.975,n-1)
# joint.mean <-apply(res.mean, 2, mean)
# joint.sd <-apply(res.mean, 2, sd)
# coverage.b1 <- colSums((res.mean[,1:4]- t05*res.sd[,1:4]<=A$b1)&(res.mean[,1:4]+t05*res.sd[,1:4]>=A$b1))/nsim
# coverage.b2 <- colSums((res.mean[,5:6]- t05*res.sd[,5:6]<=A$b2)&(res.mean[,5:6]+t05*res.sd[,5:6]>=A$b2))/nsim
# coverage.assoc <- sum((res.mean[,7]- t05*res.sd[,7]<=A$gamma)&(res.mean[,7]+t05*res.sd[,7]>=A$gamma))/nsim
# coverage.sigu<- colSums((res.mean[,8:9]- t05*res.sd[,8:9]<=A$sigu) & (res.mean[,8:9]+t05*res.sd[,8:9]>=A$sigu))/nsim
# coverage.vare <- sum((res.mean[,10]- t05*res.sd[,10]<=A$vare) & (res.mean[,10]+t05*res.sd[,10]>=A$vare))/nsim
# biasb1 <- joint.mean[1:4]-A$b1#bias for long coefficients
# biasb2 <- joint.mean[5:6]-A$b2#bias for surv coefficients
# bias.assoc <- joint.mean[7]-A$gamma#bias for gamma
# bias.sigu <- joint.mean[8:9]-A$sigu #bias for sigu
# bias.vare <- joint.mean[10]-A$vare
# mseb1 <- (joint.sd[1:4])^2+biasb1^2
# mseb2 <- (joint.sd[5:6])^2+biasb2^2
# mse.assoc <- (joint.sd[7])^2+bias.assoc^2
# mse.sigu <- (joint.sd[8:9])^2+bias.sigu^2
# mse.vare <- (joint.sd[10])^2+bias.vare^2
# true.values <- c(A$b1, A$b2, A$gamma, A$sigu, A$vare)
# bias <- c(biasb1, biasb2, bias.assoc, bias.sigu, bias.vare)
# mse <- c(mseb1, mseb2, mse.assoc, mse.sigu, mse.vare)
# coverage <- c(coverage.b1, coverage.b2, coverage.assoc, coverage.sigu, coverage.vare)
# lastres <- data.frame(colnames(res.mean), true.values, bias, mse, coverage)
# 
# 
# 
# 
