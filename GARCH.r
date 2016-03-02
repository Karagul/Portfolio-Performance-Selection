setwd("D://school/second semester/891/project/project")
library('moments')
library('nnet')
NorthA <- read.table("MSCI North America.txt")
EU <- read.table("MSCI Europe.txt")
Pac <- read.table("MSCI Pacific.txt")
LatinA <- read.table("MSCI Latin America.txt")
Asia <- read.table("MSCI Asia.txt")
p.n <- NorthA$V2
p.e <- EU$V2
p.p <- Pac$V2
p.l <- LatinA$V2
p.a <- Asia$V2
# log returns
T <- length(p.n)
r.n <- 100*log(p.n[2:T]/p.n[1:T-1])
r.e <- 100*log(p.e[2:T]/p.e[1:T-1])
r.p <- 100*log(p.p[2:T]/p.p[1:T-1])
r.l <- 100*log(p.l[2:T]/p.l[1:T-1])
r.a <- 100*log(p.a[2:T]/p.a[1:T-1])
plot(ts(r.a,frequency=260,start=c(1990)),ylab="log returns",main="Asia - Log Returns")



# growth of 1$ - normalized price series
np.n <- p.n/p.n[1]
np.e <- p.e/p.e[1]
np.p <- p.p/p.p[1]
np.l <- p.l/p.l[1]
np.a <- p.a/p.a[1]
max <- max(np.n,np.e,np.p,np.l,np.a)
min <- min(np.n,np.e,np.p,np.l,np.a)
plot(ts(np.a,frequency=260,start=c(1990)),ylim=c(min,max),ylab="normalized price series",main="Normalized Prices - np.a")

hist(r.a,xlab="Returns",breaks=200,main="Asia Returns")

tt <- length(r.a)
# GARCH(1,1)
R.ag <- matrix(0,tt)
hr.ag <- matrix(0,tt)
e1 <- rnorm(tt)
mu.a <- mean(r.a)
hr.ag[1] <- var(r.a)
a0 <- 0.1
# initial value 1
a1 <- 0.1
a2 <- 0.8

# construct GARCH(1,1) returns
for (t in 2:tt){
  hr.ag[t] <- a0+a1*((R.ag[t-1]-mu.a)^2)+a2*hr.ag[t-1]
  R.ag[t] <- mu.a+sqrt(hr.ag[t])*e1[t]
}

# MLE
like <- function(pars){
  mu.a <- pars[1]
  a0 <- pars[2]
  a1 <- pars[3]
  a2 <- pars[4]
  
  llike <- 1000000 # penalty value
  
  if (a1+a2<1){
    hr.ag <- matrix(0,tt)
    hr.ag[1] <- var(r.a)
    for (t in 2:tt){
      hr.ag[t] <- a0+a1*((R.ag[t-1]-mu.a)^2)+a2*hr.ag[t-1]
    }
    llike <- -(tt/2)*log(2*pi)-0.5*sum(log(hr.ag))-sum(0.5*((R.ag-mu.a)^2/hr.ag))
    print(c(pars,llike))
    return (-llike)
  }
  
  return (llike)
}
# parameter 1
pars0 <- c(mu.a,a0,a1,a2)
temp <- optim(pars0,like,method="L-BFGS-B",lower=c(-5,0.0001,0,0),upper=c(5,2,0.99,0.99))
pars.ag1 <- temp$par
pars.ag1


# estimated variance 1
sig2.ag1 <- matrix(0,tt)
sig2.ag1[1] <- var(r.a)
for (t in 2:tt)
  sig2.ag1[t] <- pars.ag1[2]+pars.ag1[3]*(R.ag[t-1]-mu.a)^2+pars.ag1[4]*sig2.ag1[t-1]
sig2.ag1[1]

# initial value 2
a1 <- 0.15
a2 <- 0.8
# parameter 2
pars0 <- c(mu.a,a0,a1,a2)
temp <- optim(pars0,like,method="L-BFGS-B",lower=c(-5,0.0001,0,0),upper=c(5,2,0.99,0.99))
pars.ag2 <- temp$par
pars.ag2


# estimated variance 2
sig2.ag2 <- matrix(0,tt)
sig2.ag2[1] <- var(r.a)
for (t in 2:tt)
  sig2.ag2[t] <- pars.ag2[2]+pars.ag2[3]*(R.ag[t-1]-mu.a)^2+pars.ag2[4]*sig2.ag2[t-1]
sig2.ag2[1]


# initial value 3
a1 <- 0.05
a2 <- 0.9
pars0 <- c(mu.a,a0,a1,a2)
temp <- optim(pars0,like,method="L-BFGS-B",lower=c(-5,0.0001,0,0),upper=c(5,2,0.99,0.99))
pars.ag3<- temp$par
pars.ag3

sig2.ag3 <- matrix(0,tt)
sig2.ag3[1] <- var(r.a)
for (t in 2:tt)
  sig2.ag3[t] <- pars.ag3[2]+pars.ag3[3]*(R.ag[t-1]-mu.a)^2+pars.ag3[4]*sig2.ag3[t-1]
sig2.ag3[1]

# Asia
plot(ts(sig2.ag1,frequency=260,start=c(1990)),ylab=expression(sigma^2),main="Estimated Variance - Asia") 



# North America
tt <- length(r.n)
 # GARCH(1,1)
R.ng <- r.n
hr.ng <- matrix(0,tt)
e1 <- rnorm(tt)
mu.n <- mean(r.n)
hr.ng[1] <- var(r.n)
a0 <- 0.1
a1 <- 0.1
a2 <- 0.8
pars.ng1 <- temp$par
pars.ng1

sig2.ng1 <- matrix(0,tt)
sig2.ng1[1] <- var(r.n)
for (t in 2:tt)
   sig2.ng1[t] <- pars.ng1[2]+pars.ng1[3]*(R.ng[t-1]-mu.n)^2+pars.ng1[4]*sig2.ng1[t-1]
   sig2.ng1[1]

a1 <- 0.15
a2 <- 0.8
pars.ng2 <- temp$par
pars.ng2

sig2.ng2 <- matrix(0,tt)
sig2.ng2[1] <- var(r.n)
for (t in 2:tt)
  sig2.ng2[t] <- pars.ng2[2]+pars.ng2[3]*(R.ng[t-1]-mu.n)^2+pars.ng2[4]*sig2.ng2[t-1]
  sig2.ng2[1]
  
  
a1 <- 0.05
a2 <- 0.9
pars.ng3<- temp$par
pars.ng3

ig2.ng3 <- matrix(0,tt)
sig2.ng3[1] <- var(r.n)
for (t in 2:tt)
  sig2.ng3[t] <- pars.ng3[2]+pars.ng3[3]*(R.ng[t-1]-mu.n)^2+pars.ng3[4]*sig2.ng3[t-1]
  sig2.ng3[1]

# Europe
a0 <- 0.1
a1 <- 0.1
a2 <- 0.8
pars.eg1 <- temp$par
pars.eg1
  

sig2.eg1 <- matrix(0,tt)
sig2.eg1[1] <- var(r.e)
for (t in 2:tt)
  sig2.eg1[t] <- pars.eg1[2]+pars.eg1[3]*(R.eg[t-1]-mu.e)^2+pars.eg1[4]*sig2.eg1[t-1]
  sig2.eg1[1]
  
a1 <- 0.15
a2 <- 0.8
pars.eg2 <- temp$par
pars.eg2

sig2.eg2 <- matrix(0,tt)
sig2.eg2[1] <- var(r.e)
for (t in 2:tt)
   sig2.eg2[t] <- pars.eg2[2]+pars.eg2[3]*(R.eg[t-1]-mu.e)^2+pars.eg2[4]*sig2.eg2[t-1]
   sig2.eg2[1]

   
a1 <- 0.05
a2 <- 0.9
pars.eg3<- temp$par
pars.eg3
   
sig2.eg3 <- matrix(0,tt)
sig2.eg3[1] <- var(r.e)
for (t in 2:tt)
  sig2.eg3[t] <- pars.eg3[2]+pars.eg3[3]*(R.eg[t-1]-mu.e)^2+pars.eg3[4]*sig2.eg3[t-1]
  sig2.eg3[1]

deg.n <- (R.ng-pars.ng1[1])/sqrt(sig2.ng1)
deg.e <- (R.eg -pars.eg2[1])/sqrt(sig2.eg2)
deg.a <- (R.ag-pars.ag1[1])/sqrt(sig2.ag1)
l <- length(deg.n)
deg <- matrix(0,l,3)
deg[,1] <- deg.n
deg[,2] <- deg.e
deg[,3] <- deg.a
cov.deg <- matrix(0,3,3)
for (i in 1:3){for (j in 1:3)cov.deg[i,j] <- cov(deg[,i],deg[,j])
    }
sig.deg <- matrix(0,3)
for (i in 1:3){sig.deg[i] <- sqrt(cov.deg[i,i])
    }
sig.deg
  
cov.deg
rho.deg <- matrix(0,3,3)
for(i in 1:3){for(j in 1:3)rho.deg[i,j] <- cov.deg[i,j]/(sig.deg[i]*sig.deg[j])}
rho.deg
  
# Asia
# GARCH(1,1)
R.acg <- matrix(0,tc)
hr.acg <- matrix(0,tc)
e1 <- rnorm(tc)
mu.ac <- mean(r.ac)
hr.acg[1] <- var(r.ac)

# initial value - full sample estimates
a0 <- pars.ag1[2]
a1 <- pars.ag1[3]
a2 <- pars.ag1[4]
# construct GARCH(1,1) returns
for (t in 2:tc){
  hr.acg[t] <- a0+a1*((R.acg[t-1]-mu.ac)^2)+a2*hr.acg[t-1]
  R.acg[t] <- mu.ac+sqrt(hr.acg[t])*e1[t]
}

# MLE
like <- function(pars){
  mu.ac <- pars[1]
  a0 <- pars[2]
  a1 <- pars[3]
  a2 <- pars[4]
  
  llike <- 1000000 # penalty value
  
  if (a1+a2<1){
    hr.acg <- matrix(0,tc)
    hr.acg[1] <- var(r.ac)
    for (t in 2:tc)
      hr.acg[t] <- a0+a1*((R.acg[t-1]-mu.ac)^2)+a2*hr.acg[t-1]
    llike <- -(tc/2)*log(2*pi)-0.5*sum(log(hr.acg))-sum(0.5*((R.acg-mu.ac)^2/hr.acg))
    print(c(pars,llike))
    return (-llike)
  }
  
  return (llike)
}

pars0 <- c(mu.ac,a0,a1,a2)
temp <- optim(pars0,like,method="L-BFGS-B",lower=c(-5,0.0001,0,0),upper=c(5,2,0.99,0.99))
pars.acg <- temp$par
pars.acg




   