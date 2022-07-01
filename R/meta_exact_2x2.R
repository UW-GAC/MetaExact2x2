#' @title
#' Exact inference for meta-analysis of 2x2 tables
#'
#' @description
#' A function to compute exact confidence intervals and p-values for
#' meta-analysis of 2x2 tables. The function solves for the cMLE of the log odds
#' ratio and returns exact confidence intervals and p-values using the Blaker
#' method as well as confidence intervals and p-values using the Fisher
#' Information.
#'
#' @details
#' The \code{dat} data.frame must have the following 4 columns with the
#' specified names:
#' \itemize{
#' \item xk: count in the treatment group with the outcome
#' \item tk: total count with the coutcome across both groups (treatment and
#' control)
#' \item mk: count in the treatment group
#' \item nk: count in the control group
#' }
#'
#' The function computes polynomial roots of a hypergeometric distribution for
#' each study. In some instances (e.g. for large studies), the R \code{polyroot}
#' function may have difficulty solving these roots, as indicated by non-zero
#' imaginary components (the \code{root.tol} parameter sets the tolerance for
#' identifying non-zero components.) If this occurs, you can either remove the
#' problematic study from the meta-analysis or calculate the polynomial roots
#' using another software (e.g. Mathematica) and import them using the
#' '\code{lambdas} input to this function.
#'
#' This optional \code{lambdas} input allows for importing roots for a subset of
#' studies; it should be a list object with length equal to the number of rows
#' (studies) in \code{dat}, where the kth element is a vector of lambdas for the
#' kth study in \code{dat}. For every element of \code{lambdas} where a vector
#' of roots is provided, the function will use those roots; for every element of
#' \code{lambdas} that is \code{NULL}, the function will compute the polynomial
#' roots using the \code{polyroot} function.
#'
#' @param dat A data.frame with one row per study and columns as specified in the Details.
#' @param root.tol Tolerance for checking that the imaginary components of the polynomial roots are not too far from 0; default \code{= 1e-12}. See Details.
#' @param lambdas Optional input to provide the lambdas (polynomial roots) for some studies. Should be a list object where the kth element is a vector of lambdas for the kth study in dat and NULL for other studies. The default (\code{NULL}) computes the roots for all studies. See Details.
#' @param lb.interval Optional input to set the interval for \code{uniroot} to solve the Blaker lower bound; must be a length 2 vector if used. The default (\code{NULL}) is to guess a reasonable interval based on the Fisher Information based lower bound.
#' @param ub.interval Optional input to set the interval for \code{uniroot} to solve the Blaker upper bound; must be a length 2 vector if used. The default (\code{NULL}) is to guess a reasonable interval based on the Fisher Information based upper bound.
#'
#' @return Returns a named matrix with 2 rows and 4 columns. The first row contains the values for the exact inference using the Blaker method, and the second row contains the values using the Fisher Information approach. The columns are:
#' \item{est}{The cMLE log odds ratio estimate}
#' \item{ci.lb}{The lower bound of the 95% confidence interval}
#' \item{ci.ub}{The upper bound of the 95% confidence interval}
#' \item{p}{The p-value}
#'
#' @export
meta_exact_2x2 <- function(dat, root.tol = 1e-12, lambdas = NULL, lb.interval = NULL, ub.interval = NULL){

  # check columns
  if(!all(c('xk', 'tk', 'mk', 'nk') %in% colnames(dat))){
    stop('dat must have columns named xk, tk, mk, and nk; please see help page for details')
  }

  # number of studies
  K <- nrow(dat)

  ### Solving for the cMLE log OR ###

  if(!is.null(lambdas)){
    # import provided lambdas
    if(length(lambdas) != K){
      stop('The lambdas list must be the same length as the number of rows in dat')
    }
    lamindiv <- lambdas
  }else{
    lamindiv <- vector("list", length = K)
  }

  for(k in 1:K){
    if(is.null(lamindiv[[k]])){
      root <- polyroot(dhyper(0:(dat$tk[k]), dat$mk[k], dat$nk[k], dat$tk[k]))
      if(any(abs(Im(root)) > root.tol)){
        stop('The polynomial roots of study ', k, ' could not be solved. Either remove that study or import roots using the lambdas input.')
      }
      lamindiv[[k]] <- Re(root)
    }
  }

  lamstudy <- unlist(lamindiv)
  negclik <- function(lpsi){
    if(any(lamstudy==0)){
      -poisbinom::dpoisbinom(sum(dat$xk)-1, 1/(1-lamstudy/exp(lpsi))[-which(lamstudy==0)], log_d=TRUE) }
    else{
      -poisbinom::dpoisbinom(sum(dat$xk), 1/(1-lamstudy/exp(lpsi)), log_d=TRUE)
    }
  }
  negclik <- Vectorize(negclik)
  nlm1 <- nlm(negclik, p=0)

  ### cMLE Fisher Information
  psi0 <- 1 # should this always be 1?
  hess1 <- numDeriv::hessian(negclik, nlm1$est)
  cMLE_Fisher <- c( est = nlm1$est,
                    ci.lb = nlm1$est+qnorm(0.025)/sqrt(hess1),
                    ci.ub = nlm1$est+qnorm(0.975)/sqrt(hess1),
                    p = pchisq((nlm1$est-log(psi0))^2*hess1, df=1, lower.tail=FALSE))

  ### cMLE Blaker
  # get interval for uniroot if not provided
  if(is.null(lb.interval)){
    lb.interval <- c(nlm1$est + qnorm(0.001)/sqrt(hess1), nlm1$est + qnorm(0.2)/sqrt(hess1))
  }
  bin.lb <- uniroot( function(lpsi){ .getp(lpsi, lamstudy, sum(dat$xk)) - 0.05}, lb.interval)
  if(is.null(ub.interval)){
    ub.interval <- c(nlm1$est + qnorm(0.8)/sqrt(hess1), nlm1$est + qnorm(0.999)/sqrt(hess1))
  }
  bin.ub <- uniroot( function(lpsi){ .getp(lpsi, lamstudy, sum(dat$xk)) - 0.05}, ub.interval)
  cMLE_Blaker <- c( est = nlm1$est,
                    ci.lb = bin.lb$root,
                    ci.ub = bin.ub$root,
                    p = .getp(lpsi=0, lamall=lamstudy, sum(dat$xk)))

  rbind(cMLE_Blaker, cMLE_Fisher)
}

# function for getting p-values
.getp <- function(lpsi, lamall, xplus){
  dall <- dbinom( 0:length(lamall), length(lamall), mean(1/(1-lamall/exp(lpsi))) )
  sum(dall[dall<= dall[xplus+1]])
}
