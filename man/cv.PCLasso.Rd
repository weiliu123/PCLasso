\name{cv.PCLasso}
\alias{cv.PCLasso}
\title{
Cross-validation for \code{PCLasso}
}
\description{
Perform k-fold cross validations for the PCLasso model with grouped 
covariates over a grid of values for the regularization parameter \code{lambda}.
}
\usage{
cv.PCLasso(x, y, group, penalty = c("grLasso", "grMCP", "grSCAD", "gel", 
"cMCP"), nfolds = 5, standardize = TRUE, ...)
}
\arguments{
  \item{x}{
A n x p design matrix of gene expression measurements with n samples and p 
genes, as in \code{PCLasso}.
}
  \item{y}{
The time-to-event outcome, as a two-column matrix or \code{Surv} object, as in 
\code{PCLasso}. The first column should be time on study (follow up time); the 
second column should be a binary variable with 1 indicating that the event has 
occurred and 0 indicating (right) censoring.
}
  \item{group}{
A list of groups as in \code{PCLasso}. The feature (gene) names in \code{group} 
should be consistent with the feature (gene) names in \code{x}. 
}
  \item{penalty}{
The penalty to be applied to the model. For group selection, one of grLasso,
grMCP, or grSCAD. For bi-level selection, one of gel or cMCP. See \code{grpsurv}
in the R package \code{grpreg} for details.
}
  \item{nfolds}{
The number of cross-validation folds. Default is 5.
}
  \item{standardize}{
Logical flag for \code{x} standardization, prior to fitting the model. 
Default is \code{TRUE}.
}
  \item{\dots}{
Arguments to be passed to \code{cv.grpsurv} in the R package \code{grpreg}.
}
}
\details{
The function calls \code{PCLasso} \code{nfolds} times, each time leaving out 
1/\code{nfolds} of the data. The cross-validation error is based on the 
deviance. The numbers for each outcome class are balanced across the folds; 
i.e., the number of outcomes in which y is equal to 1 is the same for each fold,
or possibly off by 1 if the numbers do not divide evenly.   
\code{cv.PCLasso} uses the approach of calculating the full Cox partial 
likelihood using the cross-validated set of linear predictors. See 
\code{cv.grpsurv} in the R package \code{grpreg} for details.
}
\value{
An object with S3 class "cv.PCLasso" containing:
    \item{cv.fit }{An object of class "cv.grpsurv"}
    \item{group.dt }{Groups with  features (genes) not included 
    in \code{x} being filtered out. }
}
\references{
PCLasso: a protein complex-based group lasso-Cox model for accurate prognosis
and risk protein complex discovery. To be published.

Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group
lasso for integrative multi-omics analysis. Journal of computational biology:
a journal of computational molecular cell biology, 22, 73-84.
}
\author{
Wei Liu
}

\seealso{
\code{\link{predict.cv.PCLasso}}
}
\examples{
library("survival")

# load data
data(ExpMatrix)
data(survData)
data(PCGroup)

x = ExpMatrix    
y = Surv(time=survData[,"time"], event=survData[,"status"])

# fit model
cv.fit1 <- cv.PCLasso(x, y, group = PCGroup, penalty = "grLasso", nfolds = 10)

}

