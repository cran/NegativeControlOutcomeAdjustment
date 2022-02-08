
\name{negativeControlOutcomeAdjustment}
\alias{negativeControlOutcomeAdjustment}
\title{NegativeControlOutcomeAdjustment}
\description{
Reducing confounding bias in observational studies of vaccine efficacy using negative control outcomes
}
\usage{
negativeControlOutcomeAdjustment(Y1, Y2, Trt, W=NULL, 
          method=c("Joint-MH", "Joint-NC", "SS-Joint"), minObsPerStrata=20)
}
\arguments{
  \item{Y1}{Binary vector for the outcome of interest.}
  \item{Y2}{Numeric vector of natural numbers for the secondary outcome, which should be unaffected by the treatment.}
  \item{Trt}{Binary vector for the treatment.}
  \item{W}{Categorical confounders to define strata. \code{W} can be a factor, vector, matrix or data frame.
           If \code{NULL}, then only the \code{Joint-NC} method will be used.}
  \item{method}{Any of \code{"Joint-MH"}, \code{"Joint-NC"}, or \code{"SS-Joint"}. See details.
                The default is all three methods.}
  \item{minObsPerStrata}{Minimum number of observations within each stratum for the \code{SS-Joint} method. 
                         Strata with fewer than \code{minObsPerStrata} observations will be
                         removed from the \code{SS-Joint} calculation.
                         If no stratum has at least \code{minObsPerStrata} observations, 
                         then the \code{SS-Joint} method will not be used.
                         The default is 20.}
}
\value{
   A data frame with columns \code{method}, \code{beta_1.hat}, \code{sd.beta_1.hat}, \code{n.strata} and
   \code{error.message}, where \code{beta_1.hat} is the "de-biased" estimate of the 
    treatment effect on \code{Y1}, 
   \code{sd.beta_1.hat} is the sandwich estimate for the standard deviation of \code{beta_1.hat},
   and \code{n.strata} is the number of strata used in the calculation.
 }
\details{
  Observations with missing or non-finite values will be removed from the data. \cr
  For the \code{SS-Joint} method, strata that yield non-finite estimates
  are removed from the calculations, and strata that contain 
  too few observations as defined by \code{minObsPerStrata} above are removed. \cr


  \bold{Joint-MH} 

  \code{Joint-MH} reduces bias due to unmeasured confounders in the estimate of the primary
   outcome relative risk, by using the treatment effect on the secondary outcome. 
  \code{ Joint-MH} first estimates jointly the effects of the treatment on the primary and 
  secondary outcomes from stratification (with Mantel-Haenszel-type weights) on the measured confounders, 
  and then uses the non-zero estimated effect on the secondary outcome to "subtract off" bias due to 
  unmeasured confounders. If the observed confounder used to stratify only has a few large strata, 
  method \code{SS-Joint} could be used instead. 
  If information on observed confounders is not provided, method \code{Joint-NC} is used instead.


  \bold{Joint-NC}

  \code{Joint-NC} reduces bias due to unmeasured confounders in the estimate of the primary outcome
    log-relative risk, by using the treatment effect on the secondary outcome. \code{Joint-NC} first estimates
    jointly the effects of the treatment on the primary and secondary outcomes, by using log-links, 
    and then uses the non-zero estimated effect on the secondary outcome to "subtract off" bias due
    to unmeasured confounders. Information on potential observed confounders is not used with this method.

  
   \bold{SS-Joint}

   \code{SS-Joint} applies method \code{Joint-NC} within each stratum of the observed confounders, 
    in order to reduce bias due to unmeasured confounders in the primary outcome log-relative risk estimate. 
    More precisely, it estimates jointly the effects of the treatment on the primary and secondary outcomes, 
    by using log-links, within each stratum. 
    It then "subtracts off" the non-zero estimated effects on the secondary outcome from the estimated effects 
    on the primary outcome, respectively in each stratum. 
    \code{SS-Joint} finally combines efficiently these stratum-specific estimates by taking a weighted combination, 
    with weights inversely proportional to their respective variances. 
    If no information on observed confounders is provided, method \code{Joint-NC} is used instead. 
    If data is available on a categorical confounder, but with many possibly sparse strata, 
    we recommend using method \code{Joint-MH} instead. 

}

%\seealso{ }

\examples{
  data(data, package="NegativeControlOutcomeAdjustment")

  Y1  <- data[, "Y1"]
  Y2  <- data[, "Y2"]
  Trt <- data[, "T"]

  # With no covariates, only the Joint-NC method is used
  negativeControlOutcomeAdjustment(Y1, Y2, Trt)

  # Age and Region define 39 strata, some of which have fewer than 20 observations. 
  # Other strata that lead to non-finite estimates in the SS-Joint method are also dropped.
  # Warning messages will be issued in these situations. 
  tab <- table(interaction(data$Age, data$Region, drop=TRUE, sep="_"))
  sum(tab < 20)
  negativeControlOutcomeAdjustment(Y1, Y2, Trt, W=data[, c("Age", "Region")])

  # Create two age groups; Age > 18 and Age <= 18 to reduce the number of strata to 6
  Age <- as.numeric(data$Age > 18)
  W   <-  interaction(Age, data$Region, sep="_", drop=TRUE)
  negativeControlOutcomeAdjustment(Y1, Y2, Trt, W=W)

}
 
\references{ 

  Etievant, L., Sampson, J.N., Gail, M.H. (2021). Increasing efficiency and reducing bias when
  assessing HPV vaccination efficacy by using non-targeted HPV strains.  Biometrics, in revision.

}