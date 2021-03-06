\name{predict.cv.PCLasso}
\alias{predict.cv.PCLasso}
\title{
Make predictions from a cross-validated PCLasso model
}
\description{
Similar to other predict methods, this function returns predictions from a 
fitted "cv.PCLasso" object, using the optimal value chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.PCLasso}(object, x = NULL, type = c("link", "response", 
"survival", "median", "norm", "coefficients", "vars", "nvars", "vars.unique", 
"nvars.unique", "groups", "ngroups"), lambda, ...)
}
\arguments{
  \item{object}{
Fitted \code{cv.PCLasso} model object.
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
requested. For values of  \code{lambda} not in the sequence of fitted models, 
linear interpolation is used.
}
  \item{\dots}{
Arguments to be passed to \code{predict.cv.grpsurv} in the R package 
\code{grpreg}.
}
}
\value{
The object returned depends on \code{type}.
}

\seealso{
\code{\link{cv.PCLasso}}
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

cv.fit1 <- cv.PCLasso(x = x.train, 
                     y = Surv(time=y.train[,"time"], event=y.train[,"status"]),
                     group = PCGroup,
                     nfolds = 5)
                     
# predict risk scores of samples in x.test
s <- predict(object = cv.fit1, x = x.test, type="link", 
             lambda=cv.fit1$cv.fit$lambda.min)

# Nonzero coefficients
sel.groups <- predict(object = cv.fit1, type="groups",
                      lambda = cv.fit1$cv.fit$lambda.min)
sel.ngroups <- predict(object = cv.fit1, type="ngroups",
                       lambda = cv.fit1$cv.fit$lambda.min)
sel.vars.unique <- predict(object = cv.fit1, type="vars.unique",
                     lambda = cv.fit1$cv.fit$lambda.min)
sel.nvars.unique <- predict(object = cv.fit1, type="nvars.unique",
                     lambda = cv.fit1$cv.fit$lambda.min)
sel.vars <- predict(object = cv.fit1, type="vars",
                    lambda=cv.fit1$cv.fit$lambda.min)
sel.nvars <- predict(object = cv.fit1, type="nvars",
                     lambda=cv.fit1$cv.fit$lambda.min)


}


