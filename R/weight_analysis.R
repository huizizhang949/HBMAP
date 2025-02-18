#' Title
#'
#' @param N
#' @param Z
#' @param omega_output
#' @param omega_J_M_output
#' @param prior
#' @param mcmc
#' @param post
#' @param mouse_indices
#' @param n_cores
#' @param cluster_type
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
local_weights_analysis <- function(N = 100,
                                   Z,
                                   omega_output,
                                   omega_J_M_output,
                                   prior, mcmc, post,
                                   mouse_indices = NULL,
                                   n_cores=NULL,
                                   cluster_type='FORK',
                                   verbose = FALSE){

  if(!is.list(mcmc)) stop("mcmc must be a list")
  if(!is.list(prior)) stop("prior must be a list")
  if(!is.list(post)) stop("post must be a list")

  if(is.null(n_cores)) {n_cores <- parallel::detectCores()-1}

  J <- max(unlist(Z))

  omega_J_matrix <- do.call(rbind, omega_output)
  mean_omega_J <- colMeans(omega_J_matrix)

  # Proportion of neurons in each cluster across all mice
  if(is.null(mouse_indices)){mouse_indices = c(1:length(Z))}
  Z <- Z[mouse_indices]
  M <- length(mouse_indices)


  w =  mean_omega_J

  # Empirical variance
  w_jm_empirical <- omega_J_M_output

  w_jm_variance_empirical <- lapply(1:length(w_jm_empirical),
                                    function(t) matrix(apply(w_jm_empirical[[t]][,mouse_indices], 1, var),
                                                       nrow = 1))


  w_jm_variance_empirical <- do.call(rbind, w_jm_variance_empirical)

  w_jm_variance_empirical <- colMeans(w_jm_variance_empirical)


  # Simulate a number of allocations under null
  cl <- parallel::makeCluster(n_cores,type=cluster_type)
  # Z_null <- lapply(1:N,
  Z_null <- pbapply::pblapply(1:N, cl=cl,
                   function(n){

                     Z_n <- lapply(1:M,
                                   function(m) extraDistr::rcat(n = length(Z[[m]]),
                                                    prob = w))

                     # Obtain posterior samples of w_jm
                     w_sampling <- HBMAP_mcmc(Y = NULL, mcmc = mcmc, prior = prior,
                                              Z.fix = Z_n, J = J, post = post, verbose = verbose)


                     w_jm_samples <- w_sampling$omega_J_M_output

                     # Compute the variance of each iteration
                     w_jm_variance <- lapply(1:length(w_jm_samples),
                                             function(t) matrix(apply(w_jm_samples[[t]], 1, var),
                                                                nrow = 1))

                     w_jm_variance <- do.call(rbind, w_jm_variance)


                     # Return posterior mean
                     matrix(colMeans(w_jm_variance), nrow = 1)
                     # }
                   })

  parallel::stopCluster(cl)
  # Combine
  Z_null <- do.call(rbind, Z_null)

  # Proportion of times the variance is larger than the null
  proportion_large <- sapply(1:J,
                             function(j) mean(w_jm_variance_empirical[j] > Z_null[,j]))

  return(data.frame('cluster' = 1:J,
                    'probability' = proportion_large))

}
