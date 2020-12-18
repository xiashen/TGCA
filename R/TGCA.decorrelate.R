#' Decorrelation Z-scores with estimated phenotypic correlations
#' 
#' The function re-weights GWAS Z-scores acorss multiple phenotypes by sample sizes and adjusts their phenotypic correlation before running TGCA.
#' 
#' @param tstat A data frame of GWAS Z-scores across K phenotypes of M SNPs.
#' @param N A vector of the sample size of K phenotypes.
#' @param MAF A vector of minor alle frequency of M SNPs.
#' @param maf.cut A value defining the low-MAF cutoff, where the SNPs with lower MAF values will be used to estimate the phenotypic correlations. Default 5e-4.
#' @param eigen.cut A value defining the proportion of information captured by the decorrelated eigenvectors. Default 0.9.
#'
#' @note You can also use this fuction to re-weight GWAS Z-scores acorss multiple phenotypes by sample size and adjust their phenotypic correlation.
#' 
#' @return A list of two elements, where \code{z.decorrelated} contains a matrix of decorrelated Z-scores corresponding to the M SNPs, and \code{cor.pheno} gives the estimated phenotypic correlation matrix.
#'
#' @author Ting Li, Xia Shen
#'  
#' @references
#' Shen, X., Li, T., Ning, Z. (2020). Improved estimation of phenotypic correlations using summary association statistics. bioRxiv.
#' 
#' @seealso
#' NULL
#'
#' 
#' @examples
#' \dontrun{
#' data(tgca)
#' decor <- TGCA.decorrelate(tstat, MAF, N)
#' image(cor(decor$z.decorrelated))
#' image(decor$cor.pheno)
#' }
#'
#' @export
#' 

TGCA.decorrelate <- function(tstat, N, MAF, maf.cut = 5e-4, eigen.cut = .9) {
    idx <- which(MAF < maf.cut)
    tstat.rp <- tstat[idx,]
    R2 <- cor(tstat.rp, use= 'pairwise.complete.obs')

    svd(R2) -> sv
    qq <- length(which(cumsum(sv$d)/nrow(R2) < eigen.cut))
    L <- sv$u[,1:qq] %*% diag(sqrt(1/sv$d))[1:qq,1:qq]

    w <- sqrt(N)/mean(sqrt(N))

    tstat_weighted <- t(t(tstat)/w)
    tstat_decor <- tstat_weighted %*% L
    return(list(z.decorrelated = tstat_decor, cor.pheno = R2))
}


