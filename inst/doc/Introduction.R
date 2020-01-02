## -----------------------------------------------------------------------------
library(SC19057)
data("law",package = "bootstrap")
x <- as.matrix(law)
xrange <- seq(from = 500, to = 700, by = 10) 
yrange <- seq(from = 1, to = 4, by = 0.05)
k <- 4 
fit <- knnde(x, k, xrange, yrange)
persp(xrange,yrange,fit,col='orange')


## -----------------------------------------------------------------------------
library(SC19057)
data("law")
x <- law$LSAT
y <- law$GPA
B <- 20
indTest(x,y,B,0.05)


## -----------------------------------------------------------------------------
library(SC19057)
library(Rcpp)
out <- rwMC(2.0,20.0,2000)
plot(out[,1], type="l", main='rwMC(sigma=2)', ylab="x",xlab = "")

## ----echo=FALSE---------------------------------------------------------------
library(knitr)
kable(women)

## -----------------------------------------------------------------------------
myfit<-lm(weight~height,women)
summary(myfit)

## ----echo=FALSE,results='hide'------------------------------------------------
b<-round(coefficients(myfit),4)

## ----women,echo=FALSE---------------------------------------------------------
myfit<-lm(weight~height,women)
plot(women$height,women$weight,xlab = "HEIGHT",ylab = "WEIGHT",main = "Regression of weight on height",col.main="red")
abline(myfit,col=("blue"))


## ----echo=FALSE,results='asis'------------------------------------------------
x<-rnorm(150)
y1<-rt(150,50)
y2<-rt(150,100)
y3<-rt(150,150)

hist(x,main="The Distribution Of X",col.main="red",col="blue",xlim =c(-3,3))
hist(y1,main="The Distribution Of y1",col.main="orange",col="green",xlim =c(-3,3))
hist(y2,main="The Distribution Of y2",col.main="orange",col="green",xlim =c(-3,3))
hist(y3,main="The Distribution Of y3",col.main="orange",col="green",xlim =c(-3,3))


## -----------------------------------------------------------------------------
set.seed(1234)
x<-function(a){
shape<-1/(2*a^2)
t<-rgamma(1e4,1,shape)
x<-sqrt(t)
return(x)
}#generate a random R(σ) sample


hist(x(0.1),prob=TRUE)
hist(x(0.5),prob=TRUE)
hist(x(2),prob=TRUE)
hist(x(5),prob=TRUE)

## -----------------------------------------------------------------------------
set.seed(1234)
n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n,3,1)
x<-function(p1){
p<-p1
u <- runif(n)
k <- as.integer(u < p) #vector of 0’s and 1’s
x <- k * x1 + (1-k) * x2 
}#generte a mixture sample

hist(x(0.75), prob=TRUE);lines(density(x(0.75)),col="red")
hist(x(0.15), prob=TRUE);lines(density(x(0.15)),col="red")
hist(x(0.3), prob=TRUE);lines(density(x(0.3)),col="red")
hist(x(0.45), prob=TRUE);lines(density(x(0.45)),col="red")
hist(x(0.6), prob=TRUE);lines(density(x(0.6)),col="red")
hist(x(0.9), prob=TRUE);lines(density(x(0.9)),col="red")

## -----------------------------------------------------------------------------
set.seed(12345)
wishart <- function(n,sigma){
  d <- nrow(sigma)
  T <- matrix(rep(0,d*d),d,d) 
  for(i in 1:d){
    for(j in 1:i){
        if(i>j) T[i,j] <- rnorm(1)
        if(i==j) T[i,j] <- sqrt(rchisq(1,df = (n-i+1)))
    }
  }
  L <- t(chol(sigma))
  M <- (L%*%T)%*%t(L%*%T)
  M
}

n <- 5
sigma <- matrix(c(1,0.7,0.7,1), nrow = 2, ncol = 2)
d <- 2
M <- wishart(n,sigma)
M

## ----results='hold'-----------------------------------------------------------
set.seed(123)
m<-1e4
t<-runif(m,0,pi/3)
int.m<-mean(sin(t))*pi/3
int.e<-cos(0)-cos(pi/3)
cat("Monte Carlo estimate=",int.m)
cat("\nexact value=",int.e)



## ----results='hold'-----------------------------------------------------------
set.seed(123)
m<-1000
x<-runif(m)
set.seed(123)
x<-runif(m)
y1<-exp(-x)/(1+x^2)
y2<-exp(x-1)/(1+(1-x)^2)
est<-mean((y1+y2)/2)
cat("estimation=",est)
cat("\nvar=",var((y1+y2)/2))
cat("\napproximate reduction in variance=",(var(y1)-var((y1+y2)/2))/var(y1))

## ----results='hold'-----------------------------------------------------------
#to calculate the new estimator with stratified importance sampling
M <- 10000 #number of replicates
k <- 5 #number of strata
r <- M / k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T1 <- numeric(k)
estimates <- matrix(0, N, 2)
fg<-function(u){
  x <- - log(1 - u * (1 - exp(-1)))
  g<-exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
  fg <- g/ (exp(-x) / (1 - exp(-1)))
  return(fg)
}
for (i in 1:N) {
  u <- runif(M) 
  estimates[i, 1] <- mean(fg(u)) #the estimate with importance sampling(f3) only
  for (j in 1:k){
    y <- runif(M/k, (j-1)/k, j/k) #the estimate by stratified importance sampling(f3)
    T1[j] <- mean(fg(y))
    estimates[i, 2] <- mean(T1)
  }
}
apply(estimates, 2, mean)
apply(estimates, 2, var)


## -----------------------------------------------------------------------------
set.seed(123)
m<-1000
c<-numeric(m)
for (i in 1:m) {
  n<-20
  x<-rchisq(n,2)
  alpha<-0.05
  u<-mean(x)+sqrt(var(x)/n)*qt(1-alpha/2,n-1)
  l<-mean(x)+sqrt(var(x)/n)*qt(alpha/2,n-1)
  if(2<=u&&2>=l) c[i]=1
}
cat("coverage probability(t-interval)=",sum(c)/m)

#The example using Chi-square interval
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
x <- rchisq(n, df = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1)
} )
cat("\ncoverage probability(Chi-square interval)=",mean(UCL > 4))


## ----results='hold'-----------------------------------------------------------
m<-1000
n<-1000
p<-c(0.025, 0.05, 0.95, 0.975)
sk0 <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
sk<-numeric(m)
for (i in 1:m) {
  x<-rnorm(n)
  sk[i]<-sk0(x)
}
q<-quantile(sk,probs =p)#quantiles of the estimator

#to compute the s.e of the estimator's quantiles
var.q<-p*(1-p)/(n*dnorm(q,0,sqrt(6/n))^2)
se.q<-sqrt(var.q/m)

sk.n<-rnorm(m,0,sqrt(6/n))
q.n<-quantile(sk.n,probs =p)#the quantiles of the large sample approximation 

print('s.e of esimators:');se.q
print('estimated quantiles:');q
print('quantiles of the large sample approximation:');q.n

## -----------------------------------------------------------------------------
set.seed(123)
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
alpha <- .1
n <- 30
m <- 2500
epsilon <- seq(0,1, .01)
N <- length(epsilon)

cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))#critical value for the skewness test

# against beta distribution
a<-seq(0,5,0.05)
pw<-numeric(N)
for (j in 1:N) { #for each a
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rbeta(n,a[j],5)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pw[j] <- mean(sktests)
}
plot(a, pw, type = "b",xlab = bquote(a), ylim = c(0,1),main='against beta distribution')

#contaminated beta distribution
pw.b<-numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    s <- sample(c(4, 1000), replace = TRUE,
                  size = n, prob = c(1-e, e))
    x <- rbeta(n,s,s)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pw.b[j] <- mean(sktests)
}
plot(epsilon, pw.b, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1),main = 'against contaminated beta distribution')

#against t(v) 
pw.t <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    d.f <- sample(c(1, 5), replace = TRUE,
                  size = n, prob = c(1-e, e))
    x <- rt(n,d.f)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pw.t[j] <- mean(sktests)
}
plot(epsilon, pw.t, type = 'p',pch='$',col='goldenrod',
     xlab = bquote(epsilon), ylim = c(0,1),main = 'aganist t-distribution')


## -----------------------------------------------------------------------------
set.seed(123)
n <- 20
alpha <- .05
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values

#Chisquare
for (j in 1:m) {
  x <- rchisq(n, 1)
  ttest <- t.test(x, alternative = "greater", mu = 1)
  p[j] <- ttest$p.value
}
p.chi <- mean(p < alpha)
se.chi <- sqrt(p.chi * (1 - p.chi) / m)
print(c(p.chi, se.chi))


#Uniform
for (j in 1:m) {
  x <- runif(n, 0,2)
  ttest <- t.test(x, alternative = "greater", mu = 1)
  p[j] <- ttest$p.value
}
p.unif <- mean(p < alpha)
se.unif <- sqrt(p.unif * (1 - p.unif) / m)
print(c(p.unif, se.unif))

#exp
for (j in 1:m) {
  x <- rexp(n, 1)
  ttest <- t.test(x, alternative = "greater", mu = 1)
  p[j] <- ttest$p.value
}
p.exp <- mean(p < alpha)
se.exp <- sqrt(p.exp * (1 - p.exp) / m)
print(c(p.exp, se.exp))


## -----------------------------------------------------------------------------
#an example
n<-10000
x<-rbinom(n,n,p=0.651)
y<-rbinom(n,n,p=0.676)
t.test(x,y,'two.sided',mu=0)

## -----------------------------------------------------------------------------
set.seed(123)
library(bootstrap)
library(lattice)
splom(scor)
cor(scor)
x<-as.matrix(scor)

library(boot) #for boot function

rho<-function(k,l){
r <- function(x, i) {
  #want correlation of columns 1 and 2
  cor(x[i,k], x[i,l])
}
theta <- boot(data = scor, statistic = r, R = 2000)
return(theta)
}
rho(1,2)
rho(3,4)
rho(3,5)
rho(4,5)

## -----------------------------------------------------------------------------
library(boot)
library(moments)
set.seed(1)
mu1<-0
mu2<-(8*gamma(11/2)/gamma(5/2)-3*5*2*5-5^3 )/(2*5)^(3/2)# compute the skewness of chisquare with freedom 5
n<-1e1
m<-1e3
boot.sk <- function(x,i) skewness(x[i])
ci.n<-ci.b<-ci.p<-matrix(NA,m,2)

#sample from N(0,1)
for(i in 1:m){
  x<-rnorm(n)
  de <- boot(data=x,statistic=boot.sk, R = 1000)
  ci <- boot.ci(de,type=c("norm","basic","perc"))
  ci.n[i,]<-ci$norm[2:3]
  ci.b[i,]<-ci$basic[4:5]
  ci.p[i,]<-ci$percent[4:5]
}
cat('norm =',mean(ci.n[,1]<=mu1 & ci.n[,2]>=mu1),
    '\nbasic =',mean(ci.b[,1]<=mu1 & ci.b[,2]>=mu1),
    '\nperc =',mean(ci.p[,1]<=mu1 & ci.p[,2]>=mu1))
cat('\nnorm(right) =',mean(ci.n[,2]<=mu1),
     '\nbasic(right) =',mean(ci.b[,2]<=mu1 ),
     '\nperc(right) =',mean(ci.p[,2]<=mu1))
cat('\nnorm(left) =',mean(ci.n[,1]>=mu1),
     '\nbasic(left) =',mean(ci.b[,1]>=mu1),
     '\nperc(left)=',mean( ci.p[,1]>=mu1))

#sample from chisqare
for(i in 1:m){
  x<-rchisq(n,5)
  de <- boot(data=x,statistic=boot.sk, R = 1000)
  ci <- boot.ci(de,type=c("norm","basic","perc"))
  ci.n[i,]<-ci$norm[2:3]
  ci.b[i,]<-ci$basic[4:5]
  ci.p[i,]<-ci$percent[4:5]
}
cat('\nnorm =',mean(ci.n[,1]<=mu2 & ci.n[,2]>=mu2),
    '\nbasic =',mean(ci.b[,1]<=mu2 & ci.b[,2]>=mu2),
    '\nperc =',mean(ci.p[,1]<=mu2 & ci.p[,2]>=mu2))
cat('\nnorm(right) =',mean(ci.n[,2]<=mu2),
     '\nbasic(right) =',mean(ci.b[,2]<=mu2 ),
     '\nperc(right) =',mean(ci.p[,2]<=mu2 ))
cat('\nnorm(left) =',mean(ci.n[,1]>=mu2),
     '\nbasic(left) =',mean(ci.b[,1]>=mu2),
     '\nperc(left)=',mean( ci.p[,1]>=mu2))

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor)
n<-nrow(scor)
c<-cov(scor)
lambda<-eigen(c)$values
theta.hat<-lambda[1]/sum(lambda)

# to compute jackknife estimates
theta.j<-numeric(n)
for (i in 1:n) {
  c.j<-cov(scor[-i,])
  l.j<-eigen(c.j)$values
  theta.j[i]<-l.j[1]/sum(l.j)
}
bias.j <- (n-1)*(mean(theta.j)-theta.hat)
se.j <- sqrt((n-1)*mean((theta.j-mean(theta.j))^2))
round(c(theta.hat=theta.hat,bias.jack=bias.j, se.jack=se.j),3)

## -----------------------------------------------------------------------------
library(knitr)
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag

# cross validation
e1 <- e2 <- e3 <- e4 <- numeric(n)
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2]*chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2]*chemical[k] +J2$coef[3]*chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2]*chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2)+I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2]*chemical[k] + J4$coef[3]*chemical[k]^2 + J4$coef[4]*chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}


L1 <- lm(magnetic~chemical)
L2 <- lm(magnetic~chemical+I(chemical^2))
L3 <- lm(log(magnetic)~chemical)
L4 <- lm(magnetic~chemical+I(chemical^2)+I(chemical^3))
# to get adjustied R^2 
r1 <- summary(L1)$adj.r.squared
r2 <- summary(L2)$adj.r.squared
r3 <- summary(L3)$adj.r.squared
r4 <- summary(L4)$adj.r.squared

results <- data.frame(Linear=c(mean(e1^2),r1),Quadratic=c(mean(e2^2),r2),
           Exponential=c(mean(e3^2),r3),Cubic=c(mean(e4^2),r4),
           row.names=c(" prediction error","adjustied R-square"))
kable(results)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

set.seed(123)
# Count Five test permutation
count5test_permutation <- function(z) {
n <- length(z)
x <- z[1:(n/2)]
y <- z[-(1:(n/2))]
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y)) 
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0) 
return(as.integer(max(c(outx, outy)) > 5))
}

permutation <- function(z,R) {
  n <- length(z)
  out <- numeric(R)
  for (r in 1: R){
      p <- sample(1:n ,n ,replace = FALSE)
      out[r] <- count5test_permutation(z[p])
  }
  sum(out)/R
}              

n1 <- 20
n2 <- 50
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 1e3

alphahat1 <- mean(replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
count5test(x, y)
}))
alphahat2 <- mean(replicate(m, expr={
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean 
y <- y - mean(y)
z <- c(x,y)
permutation(z,1000) 
})<0.05)

round(c(count5test=alphahat1,count5test_permutation=alphahat2),4)

## -----------------------------------------------------------------------------
#some function used in distance correlation
dCov <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x); m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d); M <- mean(d)
    a <- sweep(d, 1, m); b <- sweep(a, 2, m)
    b + M
  }
  A <- Akl(x); B <- Akl(y)
  sqrt(mean(A * B))
}
ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, -(1:p)] #permute rows of y
  return(nrow(z) * dCov(x, y)^2)
}

library(boot)
library(Ball)
set.seed(123)
n<-20
p.cor1<-p.cor2<-p.ball1<-p.ball2<-numeric(n)

  for(i in 1:n){
    x <- e <- matrix(rnorm(2*i), , 2)
    # for model 1
    y1<-x/4+e
    z1 <- cbind(x,y1)
    boot.obj1 <- boot(data = z1, statistic = ndCov2, R = 999,
                      sim = "permutation", dims = c(2, 2))
    # permutatin: resampling without replacement
    tb1 <- c(boot.obj1$t0, boot.obj1$t)
    p.cor1[i]<- 1-mean(tb1>=tb1[1])
    #ball
    p.ball1 [i]<- 1-bcov.test(z1[,1:2],z1[,3:4],R=999)$p.value
    
    # for model 2
    y2<-x/4*e
    z2 <- cbind(x,y2)
    boot.obj2 <- boot(data = z2, statistic = ndCov2, R = 999,
                      sim = "permutation", dims = c(2, 2))
    # permutatin: resampling without replacement
    tb2 <- c(boot.obj2$t0, boot.obj2$t)
    p.cor2[i]<- 1-mean(tb2>=tb2[1])
    #ball
    p.ball2[i] <- 1-bcov.test(z2[,1:2],z2[,3:4],R=999)$p.value
  }

plot(1:20,p.cor1,'l',xlab = 'n',ylab = 'power',main = 'model 1')
lines(1:20,p.ball1,'l',col='red')
legend(10,0.4,c('distance correlation test','ball covariance test'),col=c('black','red'),
       lty = c(1,1),cex = 0.7 )

plot(1:20,p.cor2,'l',xlab = 'n',ylab = 'power',main = 'model 2')
lines(1:20,p.ball2,'l',col='red')
legend(10,0.4,c('distance correlation test','ball covariance test'),col=c('black','red'),
       lty = c(1,1),cex = 0.7 )


## ----results='hold'-----------------------------------------------------------
#define the density fuction
dl<-function(x) 1/2*exp(-abs(x))
rw.Metropolis <- function(n, sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (dl(y) / dl(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

n <- 4 
N <- 3000
sigma <- c(.01, .1, 1, 10)
x0 <- 5
index <- 0:N


#sigama=0.01
rw1 <- rw.Metropolis(n, sigma[1], x0, N)
y1 <- rw1$x[index]
plot(y1, type="l", main=expression(sigma==0.01), ylab="x",xlab = "")
p1<-1-rw1$k/N

#sigama=0.1
rw2 <- rw.Metropolis(n, sigma[2], x0, N)
y2 <- rw2$x[index]
plot(y2, type="l", main=expression(sigma==0.1), ylab="x",xlab = "")
p2<-1-rw2$k/N

#sigama=1
rw3 <- rw.Metropolis(n, sigma[3], x0, N)
y3 <- rw3$x[index]
plot(y3, type="l", main=expression(sigma==1), ylab="x",xlab = "")
p3<-1-rw3$k/N

#sigma=10
rw4 <- rw.Metropolis(n, sigma[4], x0, N)
y4 <- rw4$x[index]
plot(y4, type="l", main=expression(sigma==10), ylab="x",xlab = "")
p4<-1-rw4$k/N

rate<-data.frame(p1,p2,p3,p4,row.names ='acceptance rate')
names(rate)<-c(0.01,0.1,1,10)
knitr::kable(rate)

## ----results='hold'-----------------------------------------------------------
x<-3 #given a random value
x1<-log(exp(x))
x2<-exp(log(x))
cat("x-log(exp(x))=",x-x1)
cat("\nx-exp(log(x))=",x-x2)
cat("\nlog(exp(x))-exp(log(x))=",x1-x2)
cat("\n")
all.equal(x,x1)
all.equal(x,x2)
all.equal(x1,x2)

## -----------------------------------------------------------------------------
#in the 11.5 case
fc<-function(a,k) sqrt(a^2*k/(1+k-a^2))
f1<-function(u) (1+u^2/(k-1))^(-k/2)
f2<-function(u) (1+u^2/k)^(-(k+1)/2)
f <- function(a){
  2*gamma(k/2)/(sqrt(pi*(k-1))*gamma((k-1)/2))*integrate(f1,0,fc(a,k-1))$value-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))*integrate(f2,0,fc(a,k))$value
}   

K <- c(4:25,100)
n <- length(K)
R <- numeric(n)
for (i in 1:n) {
  k <- K[i]
  R[i] <- uniroot(f,c(0.5,sqrt(k)/2+1))$root
}
print(R)

# in the 11.4case
S<-function(a){
  pt(sqrt((k-1)*a^2/(k-a^2)),k-1)-pt(sqrt((k*a^2)/(k+1-a^2)),k)
}
K1 <- c(4:25,100,500,1000)
n <- length(K1)
inter <- numeric(n)
for (i in 1:n) {
  k <- K1[i]
  inter[i] <- uniroot(S,c(0.5,sqrt(k)/2+1))$root
}
print(inter)

## ----results='hold'-----------------------------------------------------------

library(rootSolve)
N <- 500
na.<-28
nb.<-24
noo<-41
nab<-70
n<-sum(na.,nb.,noo,nab)
#initial estimations for p and q
p.hat<-uniroot(function(p) 2*p-p^2-(na.+nab)/n,c(0,1))$root
q.hat<-uniroot(function(q) 2*q-q^2-(nb.+nab)/n,c(0,1))$root
pa <- c(p.hat,q.hat)

# initial estimates 
tol <- .Machine$double.eps^0.5
pa.old <- pa+1
E <- numeric(N)
for(j in 1:N){
  E[j]<-2*pa[1]*na.*log(pa[1])/(2-pa[1]-2*pa[2])+2*pa[2]*nb.*log(pa[2])/(2-pa[2]-2*pa[1])+2*noo*log(1-pa[1]-pa[2])+na.*(2-2*pa[1]-2*pa[2])*log(2*pa[1]*(1-pa[1]-pa[2]))/(2-pa[1]-2*pa[2])+nb.*(2-2*pa[1]-2*pa[2])*log(2*pa[2]*(1-pa[1]-pa[2]))/(2-pa[2]-2*pa[1])+nab*log(2*pa[1]*pa[2])
  M<-function(x){
    P1<-2*pa[1]*na./((2-pa[1]-2*pa[2])*x[1])-2*noo/(1-x[1]-x[2])+na.*(2-2*pa[1]-2*pa[2])*(1-2*x[1]-x[2])/((2-pa[1]-2*pa[2])*x[1]*(1-x[1]-x[2]))-nb.*(2-2*pa[1]-2*pa[2])/((2-pa[2]-2*pa[1])*(1-x[1]-x[2]))+nab/x[1]
    P2<-2*pa[2]*nb./((2-pa[2]-2*pa[1])*x[2])-2*noo/(1-x[1]-x[2])-na.*(2-2*pa[1]-2*pa[2])/((2-pa[1]-2*pa[2])*(1-x[1]-x[2]))+nb.*(2-2*pa[1]-2*pa[2])*(1-2*x[2]-x[1])/((2-pa[2]-2*pa[1])*x[2]*(1-x[1]-x[2]))+nab/x[2]
    c(P1=P1,P2=P2)
  }
  S<-multiroot(f=M,star=c(.1,.1))
  pa<-S$root
  # update p and q
  if (sum(abs(pa-pa.old)/pa.old)<tol) break
  pa.old<-pa
}
plot(1:N,E,xlab=' ',ylab=' ','l',main='log-maximum likelihood')
cat("p=",pa[1])
cat("\nq=",pa[2])


## -----------------------------------------------------------------------------
data("mtcars")
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

#using fuction "lapply"
m1<-lapply(formulas, lm)
print(m1)

#using loop
for (x in formulas) {
  print(lm(x))
}


## -----------------------------------------------------------------------------
set.seed(123)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows,]
})
#using fuction "lapply"
m2<-lapply(seq_along(bootstraps),function(i) {lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp)})
print(m2)
#using loop
for (x in bootstraps) {
  print(lm(mpg~disp))
}



#without an anonymous function 
m2<-lapply(bootstraps, lm, formula=mpg~disp)
print(m2)
# there are another way which also not using an anonymous function
bootstraps1 <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows,c(1,3)]
})
#using fuction "lapply"(without an anonymous function)
lapply(bootstraps1, lm)
#using loop(without an anonymous function)
for (x in bootstraps1) {
  print(lm(x))
}


## -----------------------------------------------------------------------------

rsq <- function(mod) summary(mod)$r.squared
#R^2 in question 1
lapply(m1, rsq)
#R^2 in question 2
lapply(m2, rsq)

detach(mtcars)

## -----------------------------------------------------------------------------
set.seed(123)
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(seq_along(trials),function(i) {round(trials[[i]]$p.value,4)})

## -----------------------------------------------------------------------------

library(parallel)

num_cores = detectCores()
cluster = makePSOCKcluster(num_cores)

mcsapply = function(cluster,X,FUN,...){
 res = parLapply(cluster,X,FUN,...) 
 simplify2array(res)
}

system.time(mcsapply(cluster, 1:10, function(i) Sys.sleep(1)))
system.time(sapply(1:10, function(i) Sys.sleep(1)))

stopCluster(cluster)


## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
#define cpp function
cppFunction('NumericMatrix rwMC(double sigma, double x0, int N){
  NumericMatrix x(N,2);
  x(0,0) = x0;
  x(0,1) = 1;
  NumericVector u(N);
  u=runif(N);
  for (int i=1;i<N; i++) {
    double y = rnorm(1, x[i-1], sigma)[0];
    double t = exp(-abs(y)) / exp(-abs(x[i-1]));
      if (u[i] <= t) {
        x(i,0) = y; 
        x(i,1) = 1;
      }else {
        x(i,0) = x(i-1,0);
        x(i,1) = 0;
      }
  };
  return x;
}')
#define R function
rwMR <- function(sigma, x0, N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (exp(-abs(y)) / exp(-abs(x[i-1])))) x[i] = y 
    else {
      x[i] = x[i-1]
      k = k+1
    }
  }
  return(list(x = x, k = k))
}
# the plot with sigma=2
out1<-rwMC(2,20,2000)
out2<-rwMR(2,20,2000)
plot(out1[,1], type="l", main='rwMC(sigma=2)', ylab="x",xlab = "")
plot(out2$x,type="l", main='rwMR(sigma=2)', ylab="x",xlab = "")
#compare the generated random numbers
qqplot(out1[,1],out2$x,xlab = 'rwMR',ylab = 'rwMC',main='qqplot of two functions(sigma=2)')
abline(a=0,b=1)
#compare the computation time
ts <- microbenchmark(r=rwMR(2,20,2000),cpp=rwMC(2,20,2000))
summary(ts)

# the plot with sigma=10
out1<-rwMC(10,20,2000)
out2<-rwMR(10,20,2000)
plot(out1[,1], type="l", main='rwMC(sigma=10)', ylab="x",xlab = "")
plot(out2$x,type="l", main='rwMR(sigma=10)', ylab="x",xlab = "")
#compare the generated random numbers
qqplot(out1[,1],out2$x,xlab = 'rwMR',ylab = 'rwMC',main='qqplot of two functions(sigma=10)')
abline(a=0,b=1)
#compare the computation time
ts <- microbenchmark(r=rwMR(10,20,2000),cpp=rwMC(10,20,2000))
summary(ts)

# the plot with sigma=0.5
out1<-rwMC(0.5,20,2000)
out2<-rwMR(0.5,20,2000)
plot(out1[,1], type="l", main='rwMC(sigma=0.5)', ylab="x",xlab = "")
plot(out2$x,type="l", main='rwMR(sigma=0.5)', ylab="x",xlab = "")
#compare the generated random numbers
qqplot(out1[,1],out2$x,xlab = 'rwMR',ylab = 'rwMC',main='qqplot of two functions(sigma=0.5)')
abline(a=0,b=1)
#compare the computation time
ts <- microbenchmark(r=rwMR(0.5,20,2000),cpp=rwMC(0.5,20,2000))
summary(ts)

