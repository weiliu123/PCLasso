\name{predict.PCLasso}
\alias{predict.PCLasso}
\title{
Make predictions from a PCLasso model
}
\description{
Similar to other predict methods, this function returns predictions from a 
fitted \code{PCLasso} object.
}
\usage{
\method{predict}{PCLasso}(object, x = NULL, type = c("link", "response", "survival", 
"median", "norm", "coefficients", "vars", "nvars", "vars.unique", 
"nvars.unique", "groups", "ngroups"), lambda, ...)
}

\arguments{
  \item{object}{
Fitted \code{PCLasso} model object.
}
  \item{x}{
Matrix of values at which predictions are to be made. The features (genes) 
contained in \code{x} should be consistent with those contained in \code{x} 
in the \code{PCLasso} function.  Not used for type="coefficients" or for some 
of the type settings in \code{predict}.
}
  \item{type}{
Type of prediction: "link" returns the linear predictors; "response" gives
the risk (i.e., exp(link)); "vars" returns the indices for the nonzero 
coefficients; "vars.unique" returns unique features (genes) with nonzero 
coefficients (If a feature belongs to multiple groups and multiple groups are 
selected, the feature will be repeatedly selected. Compared with "var", 
"var.unique" will filter out repeated features.); "groups" returns the groups 
with at least one nonzero coefficient; "nvars" returns the number of nonzero 
coefficients; "nvars.unique" returens the number of unique features (genes) with
nonzero coefficients; "ngroups" returns the number of groups with at least one 
nonzero coefficient; "norm" returns the L2 norm of the coefficients in each 
group."survival" returns the estimated survival function; "median" estimates 
median survival times. 
}
  \item{lambda}{
Values of the regularization parameter \code{lambda} at which predictions are 
requested. For values of \code{lambda} not in the sequence of fitted models, 
linear interpolation is used.
}
  \item{\dots}{
Arguments to be passed to \code{predict.grpsurv} in the R package \code{grpreg}.
}
}
\details{
See \code{predict.grpsurv} in the R package \code{grpreg} for details.
}
\value{
The object returned depends on \code{type}.
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

set.seed(429006)
train.Idx <- sample(nrow(ExpMatrix), floor(2/3*nrow(ExpMatrix)))
x.train <- ExpMatrix[train.Idx ,]
x.test <- ExpMatrix[-train.Idx ,]
y.train <- survData[train.Idx,]
y.test <- survData[-train.Idx,]

fit1 <- PCLasso(x = x.train, 
                y = Surv(time=y.train[,"time"], event=y.train[,"status"]),
                group = PCGroup)

# predict risk scores of samples in x.test
s <- predict(object = fit1, x = x.test, type="link", lambda=fit1$fit$lambda)

s <- predict(object = fit1, x = x.test, type="link", lambda=fit1$fit$lambda[10])

s <- predict(object = fit1, x = x.test, type="link", lambda=c(0.1, 0.01))

# Nonzero coefficients
sel.groups <- predict(object = fit1, type="groups",
                      lambda = fit1$fit$lambda)
sel.ngroups <- predict(object = fit1, type="ngroups",
                       lambda = fit1$fit$lambda)
sel.vars.unique <- predict(object = fit1, type="vars.unique",
                           lambda = fit1$fit$lambda)
sel.nvars.unique <- predict(object = fit1, type="nvars.unique",
                            lambda = fit1$fit$lambda)
sel.vars <- predict(object = fit1, type="vars",
                    lambda=fit1$fit$lambda)
sel.nvars <- predict(object = fit1, type="nvars",
                     lambda=fit1$fit$lambda)


# For values of lambda not in the sequence of fitted models, 
# linear interpolation is used.
sel.groups <- predict(object = fit1, type="groups",
                      lambda = c(0.1, 0.01))
sel.ngroups <- predict(object = fit1, type="ngroups",
                       lambda = c(0.1, 0.01))
sel.vars.unique <- predict(object = fit1, type="vars.unique",
                           lambda = c(0.1, 0.01))
sel.nvars.unique <- predict(object = fit1, type="nvars.unique",
                            lambda = c(0.1, 0.01))
sel.vars <- predict(object = fit1, type="vars",
                    lambda=c(0.1, 0.01))
sel.nvars <- predict(object = fit1, type="nvars",
                     lambda=c(0.1, 0.01))

}

