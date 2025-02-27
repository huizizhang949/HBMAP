#' Compute posterior similarity matrix
#'
#' @param mcmc_run_all_output output from \code{HBMAP_mcmc}.
#'
#' @return a list containing:
#' \item{psm.within}{a list of length \eqn{M}. The \eqn{m}-th list contains \eqn{M} posterior similarity matrices
#' recording the co-clustering probabilities between the \eqn{m}-th mouse and other mice.}
#' \item{psm.combined}{a single posterior similarity matrix considering all mice together.}
#' @export
#'
similarity_matrix <- function(mcmc_run_all_output){

  # Input from MCMC
  M <- mcmc_run_all_output$M
  C <- mcmc_run_all_output$C

  C_cumsum = c(0,cumsum(C))

  # Put z in a matrix
  Zmat = matrix(unlist(mcmc_run_all_output$Z_output), length(mcmc_run_all_output$Z_output), sum(C),byrow = TRUE)

  # Posterior similarity between all observations
  psm <- mcclust::comp.psm(Zmat)

  # Posterior similarity between observations across 2 datasets
  psm_within <- list(M)

  for(m1 in c(1:M)){

    psm_within[[m1]] = list(M)

    for(m2 in c(1:M)){

      psm_within[[m1]][[m2]] = psm[(C_cumsum[m1]+1):C_cumsum[m1+1],(C_cumsum[m2]+1):C_cumsum[m2+1]]

    }
  }

  ##-------------------------------- Return both -----------------------------------------

  return(list('psm.within' = psm_within,
              'psm.combined' = psm))

}
