\name{plot.PCLasso}
\alias{plot.PCLasso}
\title{
Plot coefficients from a PCLasso object
}
\description{
Produces a plot of the coefficient paths for a fitted \code{PCLasso} object.
}
\usage{
\method{plot}{PCLasso}(x, norm = TRUE, ...)
}
\arguments{
  \item{x}{
Fitted \code{PCLasso} model.
}
  \item{norm}{
If TRUE, plot the norm of each group, rather than the individual coefficients.
}
  \item{\dots}{
Other graphical parameters to \code{plot}.
}
}
\author{
Wei Liu
}
\seealso{
\code{\link{PCLasso}}
}
\examples{
library("survival")

# load data
data(ExpMatrix)
data(survData)
data(PCGroup)

x <- ExpMatrix    
y <- Surv(time=survData[,"time"], event=survData[,"status"])

# fit the PCLasso model
fit1 <- PCLasso(x, y, group = PCGroup, penalty = "grLasso")

# plot the norm of each group
plot(fit1, norm = TRUE)

# plot the individual coefficients
plot(fit1, norm = FALSE)

}

