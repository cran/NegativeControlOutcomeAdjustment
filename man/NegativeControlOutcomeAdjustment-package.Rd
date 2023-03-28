
\name{NegativeControlOutcomeAdjustment-package}
\alias{NegativeControlOutcomeAdjustment}
\docType{package}
\title{
  Estimation of Vaccine Efficacy using Negative Control Outcomes
}
\description{
This package uses observational cohort data from a secondary outcome Y2, 
and from measured strata indicator confounders W, 
to adjust estimates of the log-relative risk of a binary outcome Y1 with a binary exposure T.  
Adjustment for Y2 reduces bias by unmeasured confounders.  
A key assumption is that Y2 is not affected by T ("negative control outcome").  
Another key assumption is that effects of unmeasured confounders on Y1 are proportional 
to their effects on Y2.
This work was motivated by studies of vaccine efficacy against vaccine-targeted strains 
in which measurements of infection against several non-vaccine-targeted strains are also available. 
In such study, Y1 would typically be the indicator of infection with a given vaccine-targeted strain, 
or the indicator of infection with at least one of several vaccine-targeted strains, 
and Y2 would be the total number of infections with the non-vaccine-targeted strains.
}
\details{
Y1 is coded 0 or 1 for absence or presence of the event (e.g. infection with one
 given vaccine-targeted strain, or infection with at least one vaccine-targeted strain).
Y2 is coded with non-negative integers for the total number of events. 
T is coded 0 for unexposed (e.g. received placebo) and 1 for exposed (e.g. received vaccine). 
W is a stratification indicator covering K possible strata. 
It can be entered as a factor, vector, matrix, or data frame. 


The package includes three methods: Joint-NC, Joint-MH, and SS-Joint. 
If there are many strata, the program recommends using the Joint-MH method, 
which works with many, possibly sparse, strata. 
If there are only a few large strata, the SS-Joint method can be used,
 which makes fewer assumptions. 
If W is null (i.e. there are no data on measured strata), the Joint-NC method can be used. 
}
\author{Lola Etievant, Mitchell H. Gail, Joshua N Sampson}
\keyword{package}
\references{
Etievant, L., Sampson, J.N., Gail, M.H. (2022). Increasing efficiency and reducing bias when
 assessing HPV vaccination efficacy by using nontargeted HPV strains.  
 Biometrics. 1-12. https://doi.org/10.1111/biom.13663
}

