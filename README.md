# PCLasso
A protein complex-based group lasso-Cox model for accurate prognosis and risk protein complex discovery.

# Details
Package: PCLasso

Type: Package

Title: A protein complex-based group lasso-Cox model for accurate prognosis and risk protein complex discovery

Version: 1.0

Date: 2021-01-29

Authors@R: c(person(given = "Wei", family = "Liu", email = "freelw@qq.com", role = c("cre", "aut"),comment = c(ORCID = "0000-0002-5496-3641")))

Depends: R (>= 3.5.0), grpreg

Imports: survival

Description: The PCLasso model is a prognostic model which selects important predictors at the protein complex level to achieve accurate prognosis and identify risk protein complexes. The PCLasso model has three inputs: a gene expression matrix, survival data, and protein complexes. It estimates the correlation between gene expression in protein complexes and survival data at the level of protein complexes.  Similar to the traditional Lasso-Cox model, PCLasso is based on the Cox PH model and estimates the Cox regression coefficients by maximizing partial likelihood with regularization penalty. The difference is that PCLasso selects features at the level of protein complexes rather than individual genes. Considering that genes usually function by forming protein complexes, PCLasso regards genes belonging to the same protein complex as a group, and constructs a l1/l2 penalty based on the sum (i.e., l1 norm) of the l2 norms of the regression coefficients of the group members to perform the selection of features at the group level. Since a gene may belong to multiple protein complexes, that is, there is overlap between protein complexes, the classical group Lasso-Cox model for non-overlapping groups may lead to false sparse solutions. The PCLasso model deals with the overlapping problem of protein complexes by constructing a latent group Lasso-Cox model. And by reconstructing the gene expression matrix of the protein complexes, the latent group Lasso-Cox model is transformed into a non-overlapping group Lasso-Cox model in an expanded space, which can be directly solved using the classical group Lasso method. Through the final sparse solution, we can predict the patient's risk score based on a small set of protein complexes and identify risk protein complexes that are frequently selected to construct prognostic models.

License: Artistic-2.0

# Index of help topics
The PCLasso model accepts a gene expression matrix, survival data, and protein complexes for the PCLasso model, and makes predictions for new samples and identifies risk protein complexes.

PCLasso Construct a PCLasso model based on a gene expression matrix, survival data, and protein complexes.

cv.PCLasso Perform k-fold cross validations for the PCLasso model with grouped covariates over a grid of values for the regularization parameter lambda, and returns an optimal value for lambda.

plot.PCLasso Produce a plot of the coefficient paths for a fitted PCLasso object.

plot.cv.PCLasso Plot the cross-validation curve from a cv.PCLasso object, along with standard error bars.

predict.PCLasso Make predictions from a PCLasso model

predict.cv.PCLasso Make predictions from a cross-validated PCLasso model, using the optimal value chosen for lambda.

ExpMatrix The expression data

PCGroup Protein complexes for PCLasso/cv.PCLasso

survData Survival data

# References
PCLasso: a protein complex-based group lasso-Cox model for accurate prognosis and risk protein complex discovery. To be published.

Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group lasso for integrative multi-omics analysis. Journal of computational biology: a journal of computational molecular cell biology, 22, 73-84.

# Examples
library("PCLasso")

library("survival")

## load data
data(ExpMatrix)

data(survData)

data(PCGroup)

set.seed(429006)

train.Idx <- sample(nrow(ExpMatrix), floor(2/3*nrow(ExpMatrix)))

x.train <- ExpMatrix[train.Idx ,]

x.test <- ExpMatrix[-train.Idx ,]

y.train <- survData[train.Idx,]

y.test <- survData[-train.Idx,]

cv.fit1 <- cv.PCLasso(x = x.train, y = Surv(time=y.train[,"time"], event=y.train[,"status"]), group = PCGroup, nfolds = 5)
                     
## plot the norm of each group
plot(cv.fit1, norm = TRUE)

## plot the individual coefficients
plot(cv.fit1, norm = FALSE)

## plot the cross-validation error (deviance)
plot(cv.fit1, type = "all")


## predict risk scores of samples in x.test
s <- predict(object = cv.fit1, x = x.test, type="link", 
             lambda=cv.fit1$cv.fit$lambda.min)

s <- predict(object = cv.fit1, x = x.test, type="link", 
        lambda=cv.fit1$cv.fit$lambda)

s <- predict(object = cv.fit1, x = x.test, type="link", 
             lambda= c(0.1, 0.01))

## Nonzero coefficients
sel.groups <- predict(object = cv.fit1, type="groups", lambda = cv.fit1$cv.fit$lambda.min)
                      
sel.ngroups <- predict(object = cv.fit1, type="ngroups", lambda = cv.fit1$cv.fit$lambda.min)
                       
sel.vars.unique <- predict(object = cv.fit1, type="vars.unique", lambda = cv.fit1$cv.fit$lambda.min)
                     
sel.nvars.unique <- predict(object = cv.fit1, type="nvars.unique", lambda = cv.fit1$cv.fit$lambda.min)                     


## For values of lambda not in the sequence of fitted models, linear interpolation is used.
sel.groups <- predict(object = cv.fit1, type="groups", lambda = c(0.1,0.05))

sel.ngroups <- predict(object = cv.fit1, type="ngroups", lambda = c(0.1,0.05))

sel.vars.unique <- predict(object = cv.fit1, type="vars.unique", lambda = c(0.1,0.05))

sel.nvars.unique <- predict(object = cv.fit1, type="nvars.unique", lambda = c(0.1,0.05))


