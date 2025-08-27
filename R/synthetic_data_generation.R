#' Generate data with noise
#'
#' @param mcmc_run_all_output output from \code{HBMAP_mcmc} for the full algorithm.
#' @param Y real data. a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param M number of mice to generate. need to be less than (or equal to) the number of mice in the real data.
#' @param C a vector of numbers of neurons to generate for each mouse.
#' @param regions.name a character vector of region names.
#' @param noise.levels a list of numeric values. For each value, the counts for \eqn{M} mice with a noise level will be generated.
#'
#' @return a list of components:
#' \item{theta}{one set of MCMC samples to generate the data.}
#' \item{synthetic_data}{counts generated from \code{theta} without noise.}
#' \item{Z_synthetic}{allocations generated from \code{theta}.}
#' \item{noisy_synthetic_data}{if \code{noise.levels} is provided, a list of data with added noise.}
#' @export
data_simulation <- function(mcmc_run_all_output,
                            Y,
                            M,
                            C,
                            regions.name = NULL,
                            noise.levels = NULL){



  length_per_chain <- length(mcmc_run_all_output$Z_output)

  J <- mcmc_run_all_output$J
  C_data <- mcmc_run_all_output$C
  R <- mcmc_run_all_output$R
  M_data <- mcmc_run_all_output$M

  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }


  # Sum of counts for each neuron in the real data
  N_CM <- lapply(1:M_data,
                 function(m) colSums(Y[[m]]))
  min_N <- min(unlist(N_CM))

  # theta
  target.index <- sample(length_per_chain,
                         size = 1,
                         replace = FALSE)


  theta <- list(q_j_star = mcmc_run_all_output$q_star_1_J_output[[target.index]],
                gamma_j_star = mcmc_run_all_output$gamma_star_1_J_output[[target.index]],
                omega = mcmc_run_all_output$omega_output[[target.index]],
                alpha = mcmc_run_all_output$alpha_output[target.index],
                Z = mcmc_run_all_output$Z[[target.index]])

  # Distribution to generate the total count for each neuron
  mu = max(mean(unlist(N_CM)-min_N),1)
  disp = mu^2/max(var(unlist(N_CM)-min_N)-mu,1)
  # Cluster-specific parameters
  mu_j_star = sapply(1:J, function(j){
    if(sum(unlist(theta$Z)==j)<1){
      mu_j = mu
    }else{
      mu_j = max(mean(unlist(N_CM)[unlist(theta$Z)==j]-min_N),1)
    }
    mu_j
  })
  disp_j_star = sapply(1:J, function(j){
    if(sum(unlist(theta$Z)==j)<5){
      disp_j = disp
    }else{
      disp_j = mu_j_star[j]^2/max(var(unlist(N_CM)[unlist(theta$Z)==j]-min_N)-mu_j_star[j], 1)
    }
  })

  # Simulate synthetic data
  Y_rep = list()
  Z_rep = list()
  omega_MJ_rep = list()

  for(m in 1:M){
    omega_m = extraDistr::rdirichlet(1, theta$omega*theta$alpha)
    omega_MJ_rep[[m]] = omega_m

    Z_rep_m <- apply(rmultinom(C[m],1,omega_m),2,function(x){which(x==1)})
    Z_rep[[m]] = Z_rep_m

    N_rep_m = rnbinom(C[m], mu = mu_j_star[Z_rep_m], size = disp_j_star[Z_rep_m]) +min_N

    Y_rep_m <- do.call(cbind, lapply(1:C[m],
                                     function(c){


                                       prob.c <- extraDistr::rdirichlet(n = 1,
                                                            theta$q_j_star[Z_rep_m[c],]*theta$gamma_j_star[Z_rep_m[c]])

                                       rmultinom(n = 1, size = N_rep_m[c], prob = prob.c)
                                     }))

    rownames(Y_rep_m) <- regions.name
    Y_rep[[m]] = Y_rep_m
  }

  # Add some noise to the Y_rep
  noisy_data_list = NULL
  if(length(noise.levels)>0){
    noisy_data_list <- lapply(noise.levels,
                              function(noise.level){

                                lapply(1:M,
                                       function(m){

                                         noisy_data_m <- lapply(1:R,
                                                                function(r){

                                                                  noisy_data_m0 <- lapply(1:ncol(Y_rep[[m]]),
                                                                                          function(i){

                                                                                            prob_mi <- rep(0, R)
                                                                                            prob_mi[r] <- (1-noise.level)
                                                                                            prob_mi[-r] <- noise.level/(R-1)

                                                                                            matrix(rmultinom(n = 1, size = Y_rep[[m]][r,i], prob = prob_mi),
                                                                                                   ncol = 1)
                                                                                          })

                                                                  do.call(cbind, noisy_data_m0)
                                                                })

                                         # Add up and round to integer
                                         noisy_data_m <- Reduce('+', noisy_data_m)
                                         rownames(noisy_data_m) <- regions.name

                                         noisy_data_m
                                       })
                              })
  }


  return(list('synthetic_data' = Y_rep,
              'noisy_synthetic_data' = noisy_data_list,
              'Z_synthetic' = Z_rep,
              'theta' = theta))



}
