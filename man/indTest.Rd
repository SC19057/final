% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indTest.R
\name{indTest}
\alias{indTest}
\title{A independence test by bootstrap method using R}
\usage{
indTest(x, y, B, alpha = 0.05)
}
\arguments{
\item{x, y}{numeric vectors of data values which need to test}

\item{B}{the number of bootstrap replications}

\item{alpha}{the significance level (default 0.05)}
}
\value{
the test result with p_value
}
\description{
This funtion performs independence test using kernel density method and bootstrap method. To use the function, you should make sure you've installed R_package "kedd"
}
\examples{
\dontrun{
library(kedd)
data("faithful")
x <- faithful$eruptions
y <- faithful$waiting
B <- 20
indTest(x,y,B,0.01)
 }
}
