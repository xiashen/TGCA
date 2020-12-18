#' Total genetic contribution assessment
#'
#' The function estimates the total genetic contribution of a series of SNPs.
#' 
#' @param zmat A matrix of decorrelated Z-scores corresponding to M SNPs.
#'
#' @note you can also use this fuction to re-weight GWAS z-scores acorss multiple phenotypes by sample size and adjust their phenotypic correlation.
#' 
#' @return A matrix of TGCA results for M SNPs with estimates, standard errors, and p-values of TGCA theta, pi0, pi-, pi+, mu-, mu+, sigma- ,and sigma+.
#'
#' 
#' @author Ting Li, Xia Shen
#' 
#' @references 
#' Li, T., ..., Ning, Z. & Shen, X. Total genetic contribution assessment of the human genome. Submitted (2020).
#'
#' @seealso 
#' NULL
#'
#' 
#' @examples 
#' \dontrun{
#' data(tgca)
#' decor <- TGCA.decorrelate(tstat, MAF, N)
#' idx <- which(MAF < 5e-4)
#' zmat <- decor$z.decorrelated[-idx,]
#' res <- TGCA.scan(zmat)
#' }
#' 
#' @export
#' 

TGCA.scan <- function(zmat) {
require(svMisc)

pimixtool <- matrix(NA, nrow(zmat), 3)
mumixtool <- matrix(NA, nrow(zmat), 2)
sigmamixtool <- matrix(NA, nrow(zmat), 2)
plmixtool <- matrix(NA, nrow(zmat), 1)
loglikelihood <- matrix(NA, nrow(zmat), 1)
restart <- matrix(NA, nrow(zmat), 1)
convergeTime <- matrix(NA, nrow(zmat), 1)

pimixtool.se <-matrix(NA, nrow(zmat), 3)
mumixtool.se <- matrix(NA, nrow(zmat), 2)
sigmamixtool.se <- matrix(NA, nrow(zmat), 2)
plmixtool.se <- matrix(NA, nrow(zmat), 1)

require(mixtools)
for (i in 1:nrow(zmat)) {
 m1 <- try(normalmixEM(zmat[i,], lambda = c(.25, .5, .25), mu = c(-1, 0, 1), sigma = c(2, 1, 2), mean.constr = c(NA, 0, NA), sd.constr = c(NA, 1, NA), maxit = 9999), silent = TRUE)
  
  
  if (!inherits(m1, 'try-error')) {
    restart[i,1] <-m1$restarts
    if (length(m1$all.loglik) < 9999 + 1) {
      pimixtool[i,]<- m1$lambda
      mumixtool[i,]<- m1$mu[-2]
      sigmamixtool[i,]<- m1$sigma[-2]
      plmixtool[i,1] <- sum(abs(m1$lambda[-2]*m1$mu[-2]))
      
      loglikelihood[i,1] <-m1$loglik
      restart[i,1] <-m1$restarts
      convergeTime[i,1] <- length(m1$all.loglik)
      progress(i/nrow(zmat)*100)
      
      param <- c(m1$lambda, m1$mu[-2], m1$sigma[-2])
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
      require(numDeriv)
      H <- hessian(mixture.loglik, param)
      se <- try(diag(sqrt(solve(-H))), silent = TRUE) # 3 lambda + 2 mu + 2 sigma
      if (!inherits(se, 'try-error')) {
        pimixtool.se[i,] <- se[1:3]
        mumixtool.se[i,] <- se[4:5]
        sigmamixtool.se[i,] <- se[6:7]
        vXY <- function(mX, mY, vX, vY) mX**2*vY + mY**2*vX + vX*vY
        plmixtool.se[i,1] <- sqrt(vXY(m1$lambda[1], m1$mu[1], se[1]**2, se[4]**2) + vXY(m1$lambda[3], m1$mu[3], se[3]**2, se[5]**2))
      }
    }
  } 

}

  est <- cbind(pimixtool, mumixtool,sigmamixtool,plmixtool,restart,loglikelihood,convergeTime)
  dimnames(est) <- list(rownames(zmat), c('pi-', 'pi0', 'pi+', 'mu-','mu+','sigma-','sigma+','theta','restart','loglik','convergeTime'))

  se <- cbind(pimixtool.se, mumixtool.se,sigmamixtool.se,plmixtool.se)
  dimnames(se) <- list(rownames(zmat), c('pi-', 'pi0', 'pi+', 'mu-','mu+','sigma-','sigma+','theta'))

  na.omit(se)-> om.se ->final.se
  final.est<-est[rownames(om.se),]
  idx<-which(final.est[,9]==0)

  final.est<-final.est[idx,]
  final.se<-se[rownames(se) %in% rownames(final.est),]

  p_value <- matrix(NA,nrow(final.est),ncol(final.se))
  dimnames(p_value) <- dimnames(final.se)
  p_value <- pchisq((final.est[,c(1:8)]**2)/(final.se**2), 1, lower.tail = FALSE)

  colnames(p_value)->nn
  nn<-nn[c(8,2,1,3,4,5,6,7)]

  colnames(p_value)<-paste0(colnames(p_value),'.P')
  colnames(final.est)<-paste0(colnames(final.est),'.est')
  colnames(final.se)<-paste0(colnames(final.se),'.se')


  df<-cbind(p_value,final.est,final.se)

  df0<-c()
  for(i in 1:length(nn)){
    cat(nn[i],'..')
    idx<-grep(nn[i],colnames(df))
    df1<-df[,idx]
    df0<-cbind(df0,df1)
  }

return(df0)
}

