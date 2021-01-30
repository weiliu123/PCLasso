# PCLasso
A protein complex-based group lasso-Cox model for accurate prognosis and risk protein complex discovery.

# Description
The PCLasso model is a prognostic model which selects important predictors at the protein complex level to achieve accurate prognosis and identify risk protein complexes. The PCLasso model has three inputs: a gene expression matrix, survival data, and protein complexes. It estimates the correlation between gene expression in protein complexes and survival data at the level of protein complexes.  Similar to the traditional Lasso-Cox model, PCLasso is based on the Cox PH model and estimates the Cox regression coefficients by maximizing partial likelihood with regularization penalty. The difference is that PCLasso selects features at the level of protein complexes rather than individual genes. Considering that genes usually function by forming protein complexes, PCLasso regards genes belonging to the same protein complex as a group, and constructs a l1/l2 penalty based on the sum (i.e., l1 norm) of the l2 norms of the regression coefficients of the group members to perform the selection of features at the group level. Since a gene may belong to multiple protein complexes, that is, there is overlap between protein complexes, the classical group Lasso-Cox model for non-overlapping groups may lead to false sparse solutions. The PCLasso model deals with the overlapping problem of protein complexes by constructing a latent group Lasso-Cox model. And by reconstructing the gene expression matrix of the protein complexes, the latent group Lasso-Cox model is transformed into a non-overlapping group Lasso-Cox model in an expanded space, which can be directly solved using the classical group Lasso method. Through the final sparse solution, we can predict the patient's risk score based on a small set of protein complexes and identify risk protein complexes that are frequently selected to construct prognostic models.

The PCLasso model accepts a gene expression matrix, survival data, and protein 
complexes for the PCLasso model, and makes predictions for new samples and 
identifies risk protein complexes.

\code{PCLasso} constructs a \code{PCLasso} model based on a gene expression 
matrix, survival data, and protein complexes.

\code{predict.PCLasso} makes predictions from a \code{PCLasso} model.

\code{cv.PCLasso} performs k-fold cross validations for the \code{PCLasso} model 
with grouped covariates over a grid of values for the regularization parameter 
\code{lambda}, and returns an optimal value for \code{lambda}.

\code{predict.cv.PCLasso} returns predictions from a fitted \code{cv.PCLasso}
object, using the optimal value chosen for \code{lambda}.

\code{plot.PCLasso} produces a plot of the coefficient paths for a fitted 
\code{PCLasso} object.

\code{plot.cv.PCLasso} plots the cross-validation curve from a \code{cv.PCLasso}
object, along with standard error bars.
}

\references{
PCLasso: a protein complex-based group lasso-Cox model for accurate prognosis
and risk protein complex discovery. To be published.

Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group
lasso for integrative multi-omics analysis. Journal of computational biology:
a journal of computational molecular cell biology, 22, 73-84.
}


