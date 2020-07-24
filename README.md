# TGCA
Total genetic contribution assessment (TGCA) based on genome-wide summary association statistics.


# Statistical Modeling
The TGCA analysis models the Z statistics from genome-wide association studies (GWAS) for each single genetic variant as drawn from a mixture of:

![](http://www.sciweavers.org/upload/Tex2Img_1595584297/eqn.png)

and thereafter makes inference on the TGCA parameter:

![](http://www.sciweavers.org/upload/Tex2Img_1595584489/eqn.png)

assessing the total genetic contribution of the variant on the analyzed set of phenotypes.

# Setup
In order to perform a TGCA analysis in **R**, we recommend installing the **mixtools** package:
```{r}
install.packages('mixtools')
```
for mixture model fitting and the **numDeriv** package:
```{r}
install.packages('numDeriv')
```
for standard error calculation.

# Example

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
Now we can apply the mixture model to the `Zstar` values for each of the two loci: 

# Phenotypic correlations

# Reference

Li T, Ning Z, Yang Z, Zhai R, Xu W, Ying K, Wang Y, Chen Y, Shen X (2020) _Submitted_.

