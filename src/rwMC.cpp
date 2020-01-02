#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler using Rcpp
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp
//' @param sigma the stand deviation of proposal distribution to the random walk Metropolis chains 
//' @param x0 the initial point
//' @param N the number of produced values in chain
//' @return a random sample of \code{the standard Laplace distribution}
//' @useDynLib SC19057
//' @examples
//' \dontrun{
//' out <- rwMC(2.0,20.0,2000)
//' plot(out[,1], type="l", main='rwMC(sigma=2)', ylab="x",xlab = "")
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix rwMC(double sigma, double x0, int N){
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
}
