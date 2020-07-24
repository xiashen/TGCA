# TGCA
Total genetic contribution assessment (TGCA) based on genome-wide summary association statistics.


### Statistical Modeling
The TGCA analysis models the Z statistics from genome-wide association studies (GWAS) for each single genetic variant as drawn from a mixture of:

![](http://www.sciweavers.org/upload/Tex2Img_1595584297/eqn.png)

and thereafter makes inference on the TGCA parameter:

![](http://www.sciweavers.org/upload/Tex2Img_1595584489/eqn.png)

assessing the total genetic contribution of the variant on the analyzed set of phenotypes.

### Setup
In order to perform a TGCA analysis in **R**, we recommend installing the **mixtools** package:
```{r}
install.packages('mixtools')
```
for mixture model fitting and the **numDeriv** package:
```{r}
install.packages('numDeriv')
```
for standard error calculation.

### Example

The link to our whole-genome TGCA analysis results is provided in our paper. 
Here using the top SNPs of the two loci _HLA-B_ and _SSH2_, we show how the TGCA inference can be reproduced based on [UK Biobank GWAS results by the Neale's lab](http://www.nealelab.is/uk-biobank).
Example files are available at this GitHub repository.

To simplify the example, we directly provide the first 200 eigenvectors of the estimated phenotypic correlation matrix of 1,376 UK Biobank phenotypes:
```{r}
L <- readRDS('R-0.5_matrix.rds')
```
which can be used to adjust for the phenotypic correlations across the Z statistics before mixture modeling:
```{r}
Z <- as.matrix(readRDS('HLA_SSH2_1376z.rds'))
Zstar <- Z %*% L
```
The `L` matrix can be applied to any variant in the genome for the 1,376 phenotypes. Now we apply the mixture model to the `Zstar` values for the _HLA-B_ locus: 
```{r}
require(mixtools)
m1 <- normalmixEM(Zstar['6:31321915:A:G',], lambda = c(.25, .5, .25), mu = c(-1, 0, 1), sigma = c(2, 1, 2), 
                  mean.constr = c(NA, 0, NA), sd.constr = c(NA, 1, NA), maxit = 9999)
```
The proportions of mixtures, ![](http://www.sciweavers.org/upload/Tex2Img_1595589650/eqn.png), are estimated as:
```{r}
m1$lambda
## [1] 0.19730629 0.49617946 0.30651425
```
and the means and variances of genetic effects as:
```{r}
m1$mu[-2]
## [1] -2.4391990  2.6229547
m1$sigma[-2]
## [1] 0.72591854 1.20545231
```
So that the TGCA parameter is:
```{r}
sum(abs(m1$lambda[-2]*m1$mu[-2]))
## [1] 1.2852423
```

The standard error of these parameters can be derived from the log-likelihood of the mixture distribution:
```{r}
mixture.loglik <- function(param, X = m1$x) {
    lambda <- param[1:3]
    mu <- param[c(4,5)]
    sigma <- param[c(6,7)]
    L = matrix(NA, nrow = length(X), ncol = 3)
    L[,1] = dnorm(X, mean = mu[1], sd = sigma[1])
    L[,2] = dnorm(X, mean = 0, sd = 1)
    L[,3] = dnorm(X, mean = mu[2], sd = sigma[2])
    L[,1] = L[,1]*lambda[1]
    L[,2] = L[,2]*lambda[2]
    L[,3] = L[,3]*lambda[3]
    return(sum(log(rowSums(L))))
}
```
where `param` is the vector containing the seven parameters in the model. The standard errors are derived as:
```{r}
require(numDeriv)
param <- c(m1$lambda, m1$mu[-2], m1$sigma[-2])
H <- hessian(mixture.loglik, param)
se <- sqrt(diag(solve(-H)))
```
Approximating the standard error of the TGCA parameter via Delta method, we have:
```{r}
vXY <- function(mX, mY, vX, vY) mX**2*vY + mY**2*vX + vX*vY
sqrt(vXY(m1$lambda[1], m1$mu[1], se[1]**2, se[4]**2) + vXY(m1$lambda[3], m1$mu[3], se[3]**2, se[5]**2))
## [1] 0.21734431
```
and the Wald test p-value of:
```{r}
pchisq(1.2852423**2/0.21734431**2, 1, lower.tail = FALSE)
## [1] 3.3513147e-09
```

### Reference

Li T, Ning Z, Yang Z, Zhai R, Xu W, Ying K, Wang Y, Chen Y, Shen X (2020). Total genetic contribution assessment across the human genome. _Submitted_.

### Contact

If you have questions, please feel free to email xia (dot) shen (at) ed (dot) ac (dot) uk.

