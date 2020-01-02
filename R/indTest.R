#' @title A independence test by bootstrap method using R
#' @description This funtion performs independence test using kernel density method and bootstrap method. To use the function, you should make sure you've installed R_package "kedd"
#' @param x,y numeric vectors of data values which need to test
#' @param B the number of bootstrap replications
#' @param alpha the significance level (default 0.05)
#' @return the test result with p_value
#' @import kedd
#' @examples
#' \dontrun{
#' library(kedd)
#' data("faithful")
#' x <- faithful$eruptions
#' y <- faithful$waiting
#' B <- 20
#' indTest(x,y,B,0.01)
#'  }
#' @export
indTest<-function(x,y,B,alpha=0.05){

n <- length(x)
h.x <- (range(x)[2]-range(x)[1])/length(x)
h.y <- (range(y)[2]-range(y)[1])/length(y)

f <- function(xstar,ystar){ # return T_0
  k.x <- matrix(rep(0,n^2),n,n) 
  k.y <- matrix(rep(0,n^2),n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      k.x[i,j] <- kedd::kernel.fun((x[i]-x[j])/h.x)$kx
      k.y[i,j] <- kedd::kernel.fun((y[i]-y[j])/h.y)$kx
      k.y[j,j] <- 0
      k.x[i,i] <- 0 
    }
  }
  f.x <- rowSums(k.x)/(n-1) # leave-one-out KDE for x
  f.y <- rowSums(k.y)/(n-1) # leave-one-out KDE for y
  f.xy <- rowSums(k.x*k.y)/(n-1) # leave-one-out KDE for (x,y)
  Ihat <- mean(f.xy)+mean(f.x)*mean(f.y)-t(f.x)%*%f.y*2/n # estimation for ISE
  s2hat <- 2/(n^2*h.x*h.y)*sum(k.x^2*k.y^2) # estimation for variance \sigma^2
  T_0 <- n*sqrt(h.x*h.y)*Ihat/sqrt(s2hat) # (Fan) T_0 assyptotically follows N(0,1)
  return(T_0)
}

T <- numeric(B)
count <- numeric(B)
for (i in 1:B) {#bootstrap to get T's B replicate
  xstar <- sample(x,n,replace = TRUE)
  ystar <- sample(y,n,replace = TRUE)
  T[i] <- f(xstar,ystar)
  count[i] <- as.integer(abs(T[i])>qnorm(1-alpha/2))
}
p <- mean(count)

if(p<=alpha){
return(print("Not Independent"))}else {
  return(print("Independent"))
}

}
