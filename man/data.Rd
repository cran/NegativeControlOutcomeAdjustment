
\name{data}
\alias{data}

\title{
  Data for examples
}
\description{
 A simulated toy data set mimicking data from a vaccine trial.
}
\details{ 
  A data frame consisting of 1000 subjects with primary outcome of interest (\code{Y1}) 
   as infected with HPV type 16, secondary outcome (\code{Y2}) as  
    the total number of infections with 20 non-targeted HPV types, 
   vaccine assignment (\code{T}), and two categorical covariates \code{Age}
   and \code{Region} with 13 and 3 levels respectively. 
}
\examples{

 data(data, package="NegativeControlOutcomeAdjustment")

 # Display some of the data
 data[1:5, ]
}
\keyword{data}
\seealso{\code{\link{negativeControlOutcomeAdjustment}}}



