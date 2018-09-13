######       ch 7 icin ilk denemeler
### normal simjoint data turetip, long measurement mean and varyansina bakiyrz. jointplot cizidiryorz robustnessi nasil diye
library(joineR)
library(dplyr)


robustness <- function(b1=b1, b2=b2, n=n, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 = theta0, theta1 = theta1, sigu=matrix(c(0.9, 0.1, 0.1, 0.3),2,2), vare=1){
  mydata <- simjoint(b1=b1, b2=b2, n=n, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 = theta0, theta1 = theta1, sigu=sigu, vare=vare)
  longdat <- mydata$longitudinal
  survdat<- mydata$survival
  n=dim(survdat)[1]
  longdat2 <- subset(longdat, select=c("id","Y","time"))
  survdat2 <- subset(survdat, select=c("id","survtime","cens"))
  covdat <- subset(survdat, select=c("id", "ctsx", "binx"))
  jointsimdat <-  jointdata(longitudinal=longdat2, baseline=covdat, survival=survdat2, id.col="id", time.col="time")
  gnc <- longdat %>% group_by(id) %>% summarise(ltime = max(ltime))
  gnc$ltime <- gnc$ltime+1
  gnc <- data.frame(gnc)
  fit<- joint(data=jointsimdat, long.formula=Y~ctsx+binx+time,
              surv.formula=Surv(survtime,cens)~(ctsx+binx), model="intslope")
  fitSE <- jointSE(fit, n.boot = 100)
  list(fitSE=fitSE, meandat=mean(gnc$ltime), vardat=var(gnc$ltime), n=n, jointsimdat=jointsimdat)
}

simrun <- function(n = n, nsim=nsim, model = "intslope", ntms = 10, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6), gamma = 0.4, theta0 = -4, theta1 = 0.5, censlam = exp(-3), sigu=matrix(c(0.9, 0.1, 0.1, 0.3),2,2), vare=1){
   listofdfs <- list()
  for(i in 1:nsim) {
      listofdfs[[i]] <- robustness(b1=b1, b2=b2, n=n, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 =theta0, theta1 =theta1,sigu=sigu, vare=vare)
      print(i)
  } 
  list(result=listofdfs, b1=b1, b2=b2, gamma=gamma)
}
ptm <- proc.time()
A1000 <- simrun(n=1000,nsim=100,  gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6))
A500 <- simrun(n=500,nsim=100,  gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6))
A250 <- simrun(n=250,nsim=100,  gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6), sigu=matrix(c(0.9, 0.1, 0.1, 0.3),2,2), vare=1)
stoptime<- proc.time() - ptm # Stop the clock
time.ch7 <- stoptime[3]/60 #in minutes.
print(time.ch7)

save.image("ch7scen1.RData")

datam <- robustness(n=1000, gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6), sigu=matrix(c(0.9, 0.1, 0.1, 0.3),2,2), vare=1)



# A1000$result[[1]]$fitSE
# 
# data1 <- robustness(n=500, gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5)
#   
jointplot(datam$jointsimdat, Y.col="Y",Cens.col = "cens",mean.profile = TRUE, split = TRUE,ylab="Y",gp1lab="Censored",gp2lab="Event")
datam$meandat
datam$vardat
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
robustness.deneme <- function(n=n,b1=b1, b2=b2, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 = theta0, theta1 = theta1, sigu=sigu, vare=vare){
  mydata <- simjoint(n=n, b1=b1, b2=b2, gamma=gamma, model = model, ntms=ntms, censlam = censlam, theta0 = theta0, theta1 = theta1, sigu=sigu, vare=vare)
  longdat <- mydata$longitudinal
  survdat<- mydata$survival
  n=dim(survdat)[1]
  longdat2 <- subset(longdat, select=c("id","Y","time"))
  survdat2 <- subset(survdat, select=c("id","survtime","cens"))
  covdat <- subset(survdat, select=c("id", "ctsx", "binx"))
  jointsimdat <-  jointdata(longitudinal=longdat2, baseline=covdat, survival=survdat2, id.col="id", time.col="time")
  gnc <- longdat %>% group_by(id) %>% summarise(ltime = max(ltime))
  gnc$ltime <- gnc$ltime+1
  gnc <- data.frame(gnc)
#  fit<- joint(data=jointsimdat, long.formula=Y~ctsx+binx+time,
  #            surv.formula=Surv(survtime,cens)~(ctsx+binx), model="intslope")
#  fitSE <- jointSE(fit, n.boot = 100)
  list(meandat=mean(gnc$ltime), vardat=var(gnc$ltime), n=n, jointsimdat=jointsimdat)
}

data3 <-robustness.deneme(n=1000, gamma=1, model = "intslope", ntms=10, censlam = exp(-2), theta0 = -2, theta1 = 0.5, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6), sigu=matrix(c(0.9, 0.1, 0.1, 0.3),2,2), vare=1)

jointplot(data3$jointsimdat, Y.col="Y",Cens.col = "cens",mean.profile = TRUE, split = TRUE,ylab="Y",gp1lab="Censored",gp2lab="Event")
data3$meandat
data3$vardat
n <- 250
ll.est250 <- t(sapply(A250$result, function(m) m$fitSE$Estimate))
ll.se250 <- t(sapply(A250$result, function(m) m$fitSE$SE))
n <- 500
ll.est500 <- t(sapply(A500$result, function(m) m$fitSE$Estimate))
ll.se500 <- t(sapply(A500$result, function(m) m$fitSE$SE))
n <- 1000
ll.est1000 <- t(sapply(A1000$result, function(m) m$fitSE$Estimate))
ll.se1000 <- t(sapply(A1000$result, function(m) m$fitSE$SE))

A <- list(gamma= 1, b1 = c(7, -1, 0.5, 0.1), b2 = c(0.6, -0.6), sigu=diag(2), vare=0.01)
nsim <- 100


meandata250 <- sapply(A250$result, function(m) m$meandat)
vardata250 <- sapply(A250$result, function(m) m$vardat)

meandata500 <- sapply(A500$result, function(m) m$meandat)
vardata500 <- sapply(A500$result, function(m) m$vardat)

meandata1000 <- sapply(A1000$result, function(m) m$meandat)
vardata1000 <- sapply(A1000$result, function(m) m$vardat)

eventrate250 <- sapply(A250$result, function(m) 100*sum(m$jointsimdat$survival$cens)/250)
eventrate500 <- sapply(A500$result, function(m) 100*sum(m$jointsimdat$survival$cens)/500)
eventrate1000 <- sapply(A1000$result, function(m) 100*sum(m$jointsimdat$survival$cens)/1000)

# longitudinal time sekline bakarken lazim olcak bu kod
AA2 <- A1000$result[[1]]$jointsimdat$longitudinal %>% group_by(id) %>% summarise(time = max(time)) +1
hist(AA2$time)


b1 <- matrix(rep(A$b1,each=nsim),ncol=length(A$b1), byrow=FALSE)
b2 <- matrix(rep(A$b2,each=nsim),ncol=length(A$b2), byrow=FALSE)
gamma <- matrix(rep(A$gamma,each=nsim),ncol=length(A$gamma), byrow=FALSE)
sigu <- matrix(rep(diag(A$sigu),each=nsim),ncol=length(diag(A$sigu)), byrow=FALSE)
vare <- matrix(rep(A$vare,each=nsim),ncol=length(A$vare), byrow=FALSE)



res.mean <- ll.est1000
res.sd <- ll.se1000
colnames(res.mean) <- colnames(res.sd)<- c("l.int", "l.cont", "l.bin", "l.time", "s.cont", "s.bin", "gamma", "sigu1", "sigu2", "vare")

t05=qt(0.975,n-1)
joint.mean <-apply(res.mean, 2, mean)
joint.sd <-apply(res.mean, 2, sd)
coverage.b1 <- colSums(((res.mean[,1:4]- t05*res.sd[,1:4])<=b1) & ((res.mean[,1:4]+t05*res.sd[,1:4])>b1))/nsim
coverage.b2 <- colSums((res.mean[,5:6]- t05*res.sd[,5:6]<=b2)&(res.mean[,5:6]+t05*res.sd[,5:6]>=b2))/nsim
coverage.assoc <- sum((res.mean[,7]- t05*res.sd[,7]<=gamma)&(res.mean[,7]+t05*res.sd[,7]>=gamma))/nsim
coverage.sigu<- colSums((res.mean[,8:9]- t05*res.sd[,8:9]<=sigu) & (res.mean[,8:9]+t05*res.sd[,8:9]>=sigu))/nsim
coverage.vare <- sum((res.mean[,10]- t05*res.sd[,10]<=vare) & (res.mean[,10]+t05*res.sd[,10]>=vare))/nsim
biasb1 <- joint.mean[1:4]-A$b1#bias for long coefficients
biasb2 <- joint.mean[5:6]-A$b2#bias for surv coefficients
bias.assoc <- joint.mean[7]-A$gamma#bias for gamma
bias.sigu <- joint.mean[8:9]-diag(A$sigu) #bias for sigu
bias.vare <- joint.mean[10]-A$vare
mseb1 <- (joint.sd[1:4])^2+biasb1^2
mseb2 <- (joint.sd[5:6])^2+biasb2^2
mse.assoc <- (joint.sd[7])^2+bias.assoc^2
mse.sigu <- (joint.sd[8:9])^2+bias.sigu^2
mse.vare <- (joint.sd[10])^2+bias.vare^2
true.values <- c(A$b1, A$b2, A$gamma, diag(A$sigu), A$vare)
bias <- c(biasb1, biasb2, bias.assoc, bias.sigu, bias.vare)
mse <- c(mseb1, mseb2, mse.assoc, mse.sigu, mse.vare)
coverage <- c(coverage.b1, coverage.b2, coverage.assoc, coverage.sigu, coverage.vare)
lastres <- data.frame(colnames(res.mean), true.values, bias, mse, coverage)



# res.mean <- ll.est250
# res.sd <- ll.se250
# colnames(res.mean) <- colnames(res.sd)<- c("l.int", "l.cont", "l.bin", "l.time", "s.cont", "s.bin", "gamma", "sigu1", "sigu2", "vare")
# 
#   t05=qt(0.975,n-1)
#   joint.mean <-apply(res.mean, 2, mean)
#   joint.sd <-apply(res.mean, 2, sd)
#   coverage.b1 <- colSums((res.mean[,1:4]- t05*res.sd[,1:4]<=A$b1)&(res.mean[,1:4]+t05*res.sd[,1:4]>=A$b1))/nsim
#   coverage.b2 <- colSums((res.mean[,5:6]- t05*res.sd[,5:6]<=A$b2)&(res.mean[,5:6]+t05*res.sd[,5:6]>=A$b2))/nsim
#   coverage.assoc <- sum((res.mean[,7]- t05*res.sd[,7]<=A$gamma)&(res.mean[,7]+t05*res.sd[,7]>=A$gamma))/nsim
#   coverage.sigu<- colSums((res.mean[,8:9]- t05*res.sd[,8:9]<=A$sigu) & (res.mean[,8:9]+t05*res.sd[,8:9]>=A$sigu))/nsim
#   coverage.vare <- sum((res.mean[,10]- t05*res.sd[,10]<=A$vare) & (res.mean[,10]+t05*res.sd[,10]>=A$vare))/nsim
#   biasb1 <- joint.mean[1:4]-A$b1#bias for long coefficients
#   biasb2 <- joint.mean[5:6]-A$b2#bias for surv coefficients
#   bias.assoc <- joint.mean[7]-A$gamma#bias for gamma
#   bias.sigu <- joint.mean[8:9]-A$sigu #bias for sigu
#   bias.vare <- joint.mean[10]-A$vare
#   mseb1 <- (joint.sd[1:4])^2+biasb1^2
#   mseb2 <- (joint.sd[5:6])^2+biasb2^2
#   mse.assoc <- (joint.sd[7])^2+bias.assoc^2
#   mse.sigu <- (joint.sd[8:9])^2+bias.sigu^2
#   mse.vare <- (joint.sd[10])^2+bias.vare^2
#   true.values <- c(A$b1, A$b2, A$gamma, A$sigu, A$vare)
#   bias <- c(biasb1, biasb2, bias.assoc, bias.sigu, bias.vare)
#   mse <- c(mseb1, mseb2, mse.assoc, mse.sigu, mse.vare)
#   coverage <- c(coverage.b1, coverage.b2, coverage.assoc, coverage.sigu, coverage.vare)
#   lastres <- data.frame(colnames(res.mean), true.values, bias, mse, coverage)
