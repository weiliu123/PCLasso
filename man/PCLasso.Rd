\name{PCLasso}
\alias{PCLasso}
\title{
Protein complex-based group lasso-Cox model
}
\description{
Construct a PCLasso model based on a gene expression matrix, survival data, and 
protein complexes.
}
\usage{
PCLasso(x, y, group, penalty = c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"), 
standardize = TRUE, ...)
}
\arguments{
  \item{x}{
A n x p matrix of gene expression measurements with n samples and p genes.
}
  \item{y}{
The time-to-event outcome, as a two-column matrix or \code{Surv} object. The 
first column should be time on study (follow up time); the second column should 
be a binary variable with 1 indicating that the event has occurred and 0 
indicating (right) censoring.
}
  \item{group}{
A list of groups. The feature (gene) names in \code{group} should be consistent
with the feature (gene) names in \code{x}. 
}
  \item{penalty}{
The penalty to be applied to the model. For group selection, one of grLasso,
grMCP, or grSCAD. For bi-level selection, one of gel or cMCP. See \code{grpsurv}
in the R package \code{grpreg} for details.
}
  \item{standardize}{
Logical flag for \code{x} standardization, prior to fitting the model. 
Default is \code{TRUE}.
}
  \item{\dots}{
Arguments to be passed to \code{grpsurv} in the R package \code{grpreg}.
}
}
\details{
The PCLasso model is a prognostic model which selects important predictors at
the protein complex level to achieve accurate prognosis and identify risk
protein complexes. The PCLasso model has three inputs: a gene expression matrix,
survival data, and protein complexes. It estimates the correlation between gene
expression in protein complexes and survival data at the level of protein
complexes.  Similar to the traditional Lasso-Cox model, PCLasso is based on the
Cox PH model and estimates the Cox regression coefficients by maximizing partial
likelihood with regularization penalty. The difference is that PCLasso selects
features at the level of protein complexes rather than individual genes. 
Considering that genes usually function by forming protein complexes, PCLasso 
regards genes belonging to the same protein complex as a group, and constructs 
a l1/l2 penalty based on the sum (i.e., l1 norm) of the l2 norms of the 
regression coefficients of the group members to perform the selection of 
features at the group level. Since a gene may belong to multiple protein 
complexes, that is, there is overlap between protein complexes, the classical
group Lasso-Cox model for non-overlapping groups may lead to false sparse 
solutions. The PCLasso model deals with the overlapping problem of protein 
complexes by constructing a latent group Lasso-Cox model. And by reconstructing 
the gene expression matrix of the protein complexes, the latent group Lasso-Cox
model is transformed into a non-overlapping group Lasso-Cox model in an expanded
space, which can be directly solved using the classical group Lasso method. 
Through the final sparse solution, we can predict the patient's risk score based
on a small set of protein complexes and identify risk protein complexes that are
frequently selected to construct prognostic models.

}
\value{
    An object with S3 class "PCLasso" containing:
    \item{fit }{An object of class "grpsurv"}
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
    \code{\link{predict.PCLasso}}, \code{\link{cv.PCLasso}}
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

}
