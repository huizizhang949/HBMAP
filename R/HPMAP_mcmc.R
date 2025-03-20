#' A hierarchical Bayesian mixture for high-throughput axon projections
#'
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' Mandatory to run the full algorithm, or the post-processing step for sampling \eqn{q} and \eqn{\gamma}.
#' @param J integer. Number of clusters to be found. Mandatory to run the full algorithm.
#' @param mcmc a list of MCMC arguments:
#' \itemize{
#' \item \code{number_iter}, number of iterations. Default to 6000.
#' \item \code{burn_in}, number of iterations to discard as burn-in. Default to 3000.
#' \item \code{thinning}, after burn-in, save the sample at every \code{thinning} iterations. Default to 5.
#' \item \code{adaptive_prop}, a small positive constant used in adaptive MCMC algorithm.
#' A larger value will propose a new sample more likely to be different from the previous sample. Default to 0.001.
#' \item \code{auto_save, save_path, save_frequency}, If \code{auto_save} is \code{TRUE} (default is \code{FALSE}), intermediate results will be saved
#' to the path \code{save_path} at every \code{save_frequency} iterations (default is 100).
#' \code{save_path} is mandatory if \code{auto_save} is \code{TRUE}, and should contain the name of the saved results.
#' }
#' @param prior a list giving the prior information:
#' \itemize{
#' \item \code{a_gamma, b_gamma, lb_gamma}, shape, rate parameters and lower bound for \eqn{\gamma_j} in the truncated Gamma prior.
#' Default values are the median of neurons' total projection counts, 2 and 1.
#' \item \code{a, tau, nu}, parameters in the repulsive prior of \eqn{\bm{q}_{1:J}}.
#' The concentration parameters in the Dirichlet prior are given by \code{a}*(average empirical projection strength across mice).
#' \code{tau,nu} are relevant to the repulsive term, where \code{tau} represents the separation between components
#' and \code{nu} denotes how strong the repulsion is. Default values are 2, 0.4, 1/20.
#' \item \code{a_alpha, b_alpha}, shape and rate parameters for the concentration parameter in the mouse-specific Dirichlet process. Both are default to 1.
#' \item \code{a_alpha0, b_alpha0}, shape and rate parameters for the concentration parameter in the baseline Dirichlet process. Both are default to 1.
#' }
#' @param Z.fix a list of allocations (integers). Each is a vector of allocations for individual mouse. Mandatory to run the post-processing step with a fixed clustering.
#' @param post a list of arguments for the post-processing step:
#' \itemize{
#' \item \code{run_q_gamma}, boolean. Whether to sample \eqn{\bm{q}_{1:J}} and \eqn{\bm{\gamma}_{1:J}}. Default to \code{TRUE}.
#' \item \code{run_omega}, boolean. Whether to sample local and global weights. Default to \code{TRUE}.
#' }
#' @param Z.init Initial values for the allocations to run the full algorithm (optional).
#' @param verbose Default to \code{TRUE}, show progress bar.
#'
#' @return The output is a list containing the following components:
#' \item{output_index}{total number of saved MCMC samples.}
#' \item{allocation_probability}{for the full algorithm, the allocation probabilities for each neuron. It is \code{NULL} for the post-processing step.}
#' \item{Z_output}{for the full algorithm, a list of samples for allocations. Each element of the list is a list of allocations for all mice.}
#' \item{omega_J_M_output}{a list of samples for the local weights. Each element of the list is a \eqn{J \times M} matrix.}
#' \item{omega_output}{a list of samples for the global weights. Each element of the list is a vector of length \eqn{J}.}
#' \item{alpha_output, alpha_zero_output}{a vector of samples for two concentration parameters, respectively.}
#' \item{q_star_1_J_output}{a list of samples for projection strengths \eqn{\bm{q}_{1:J}}. Each element of the list is a \eqn{J \times R} matrix.}
#' \item{gamma_star_1_J_output}{a list of samples for dispersion \eqn{\bm{\gamma}_{1:J}}. Each element of the list is a vector of length \eqn{J}.}
#' \item{acceptance_prob}{average acceptance probabilities over iterations for global weights,
#' two concentration parameters in the Dirichlet processes, and cluster-specific projection strengths \eqn{\bm{q}_{1:J}} and dispersion \eqn{\bm{\gamma}_{1:J}}.}
#' \item{M, C, R, J}{number of mice, a vector of numbers of neurons for mice, number of regions, number of clusters.}
#' \item{Z}{for the post-processing step, save the provided fixed clustering. It is \code{NULL} for the full algorithm.}
#' \item{post}{a boolean value indicating whether the post-processing is run or not.}
#' @export
#'
#' @examples
#' mcmc_list = list(number_iter = 5000, thinning = 1, burn_in = 4000, adaptive_prop = 0.0001,
#'       auto_save = FALSE, save_path = NULL, save_frequency = 1000)
#' set.seed(3)
#' mcmc_all_run <- HBMAP_mcmc(Y = Y, J = J, mcmc = mcmc_list, verbose = FALSE)
HBMAP_mcmc <- function(Y = NULL,
                       J = NULL,
                       mcmc = list(),
                       prior = list(),
                       Z.fix = NULL,
                       post = NULL,
                       Z.init = NULL,
                       verbose = TRUE){

  if((is.null(J) || is.null(Y)) & is.null(Z.fix)){
    stop('To sample allocations, please provide Y and J. To run the post-processing step, please provide Z.fix, and Y is only needed to sample q and gamma.')
  }else if(!is.null(Z.fix)){
    run_post <- TRUE
    if(is.null(J)){
      J <- length(unique(unlist(Z.fix)))
    }
    # if(!is.null(J)){
    #   if(J != length(unique(unlist(Z.fix)))) {
    #     stop('J does not match the number of clusters in Z.fix')
    #   }
    # }else{
    #   J <- length(unique(unlist(Z.fix)))
    # }
  }else {
    run_post <- FALSE
  }

  if(!is.list(mcmc)) stop("mcmc must be a list")
  if(!is.list(prior)) stop("prior must be a list")
  if(!is.null(post) & !is.list(post)) stop("post must be a list")

  if(!is.null(mcmc$number_iter) && (!is.numeric(mcmc$number_iter) | (mcmc$number_iter<1))) stop("mcmc$number_iter must be a positive integer")
  if(!is.null(mcmc$thinning) && (!is.numeric(mcmc$thinning) | (mcmc$thinning<1)) && (mcmc$thinning>mcmc$number_iter)) stop("mcmc$thinning must be a positive integer less than number_iter")
  if(!is.null(mcmc$burn_in) && (!is.numeric(mcmc$burn_in) | (mcmc$burn_in<1)) && (mcmc$burn_in>mcmc$number_iter)) stop("mcmc$burn_in must be a positive integer less than number_iter")
  if(!is.null(mcmc$adaptive_prop) && (!is.numeric(mcmc$adaptive_prop) | (mcmc$adaptive_prop<=0))) stop("mcmc$adaptive_prop must be a numeric positive value")
  if(!is.null(mcmc$auto_save) && !is.logical(mcmc$auto_save)) stop("mcmc$auto_save must be a logical value")
  if(!is.null(mcmc$auto_save) && mcmc$auto_save && is.null(mcmc$save_path)) stop("mcmc$save_path must be provided when mcmc$auto_save is TRUE")
  if(!is.null(mcmc$save_frequency) && ((!is.numeric(mcmc$save_frequency) | (mcmc$save_frequency<1)))) stop("mcmc$save_frequency must be a positive integer")

  if(!is.null(prior$a_gamma) && (!is.numeric(prior$a_gamma) | (prior$a_gamma<=0))) stop("prior$a_gamma must be a numerical positive value")
  if(!is.null(prior$b_gamma) && (!is.numeric(prior$b_gamma) | (prior$b_gamma<=0))) stop("prior$b_gamma must be a numerical positive value")
  if(!is.null(prior$lb_gamma) && (!is.numeric(prior$lb_gamma) | (prior$lb_gamma<=0))) stop("prior$lb_gamma must be a numerical positive value")
  if(!is.null(prior$a) && (!is.numeric(prior$a) | (prior$a<=0))) stop("prior$a must be a numerical positive value")
  if(!is.null(prior$tau) && (!is.numeric(prior$tau) | (prior$tau<=0))) stop("prior$tau must be a numerical positive value")
  if(!is.null(prior$nu) && (!is.numeric(prior$nu) | (prior$nu<=0))) stop("prior$nu must be a numerical positive value")
  if(!is.null(prior$a_alpha) && (!is.numeric(prior$a_alpha) | (prior$a_alpha<=0))) stop("prior$a_alpha must be a numerical positive value")
  if(!is.null(prior$b_alpha) && (!is.numeric(prior$b_alpha) | (prior$b_alpha<=0))) stop("prior$b_alpha must be a numerical positive value")
  if(!is.null(prior$a_alpha0) && (!is.numeric(prior$a_alpha0) | (prior$a_alpha0<=0))) stop("prior$a_alpha0 must be a numerical positive value")
  if(!is.null(prior$b_alpha0) && (!is.numeric(prior$b_alpha0) | (prior$b_alpha0<=0))) stop("prior$b_alpha0 must be a numerical positive value")

  if(!is.null(post$run_omega) && !is.logical(post$run_omega)) stop("post$run_omega must be a logical value")
  if(!is.null(post$run_q_gamma) && !is.logical(post$run_q_gamma)) stop("post$run_q_gamma must be a logical value")

  # mcmc settings (default)
  if(is.null(mcmc$number_iter)) {number_iter = 6000} else {number_iter = mcmc$number_iter}
  if(is.null(mcmc$thinning)) {thinning = 5} else {thinning = mcmc$thinning}
  if(is.null(mcmc$burn_in)) {burn_in = 3000} else {burn_in = mcmc$burn_in}
  if(is.null(mcmc$adaptive_prop)) {adaptive_prop = 0.001} else {adaptive_prop = mcmc$adaptive_prop}
  if(is.null(mcmc$auto_save)) {auto_save = FALSE} else {auto_save = mcmc$auto_save}
  if(is.null(mcmc$save_frequency)) {save_frequency = 100} else {save_frequency = mcmc$save_frequency}
  save_path = mcmc$save_path

  # prior settings (default)
  if(is.null(prior$a_gamma)) {a_gamma = median(unlist(lapply(Y,colSums)))} else {a_gamma = prior$a_gamma}
  if(is.null(prior$b_gamma)) {b_gamma = 2} else {b_gamma = prior$b_gamma}
  if(is.null(prior$lb_gamma)) {lb_gamma = 1} else {lb_gamma = prior$lb_gamma}
  if(is.null(prior$a)) {a = 2} else {a = prior$a}
  if(is.null(prior$tau)) {tau = 0.4} else {tau = prior$tau}
  if(is.null(prior$nu)) {nu = 1/20} else {nu = prior$nu}
  if(is.null(prior$a_alpha)) {a_alpha = 1} else {a_alpha = prior$a_alpha}
  if(is.null(prior$b_alpha)) {b_alpha = 1} else {b_alpha = prior$b_alpha}
  if(is.null(prior$a_alpha0)) {a_alpha0 = 1} else {a_alpha0 = prior$a_alpha0}
  if(is.null(prior$b_alpha0)) {b_alpha0 = 1} else {b_alpha0 = prior$b_alpha0}

  # post-processing step (default)
  if(!is.null(post) && is.null(post$run_omega)) {run_omega = TRUE} else {run_omega = post$run_omega}
  if(!is.null(post) && is.null(post$run_q_gamma)) {run_q_gamma = TRUE} else {run_q_gamma = post$run_q_gamma}

  if(run_post && run_q_gamma && is.null(Y)){
    stop('Data Y is needed to sample q and gamma in the post-processing step')
  }
  # --------------- dimensions --------------------

  if(!is.null(Y)){
    M <- length(Y)
    ##-- Need to make sure that all datasets have the same number of regions
    stopifnot(length(unique(sapply(1:M, function(m) nrow(Y[[m]]))))==1)

    ##-- Number of regions
    R <- nrow(Y[[1]])

    ##-- Number of neuron cells - vector of length M
    C <- unlist(lapply(1:M, function(m) ncol(Y[[m]])))

  }else{
    M <- length(Z.fix)

    C <- sapply(Z.fix, length)

    R <- NULL
  }

  # ------------------------ Prepare for outputs -----------------

  Z_output <- NULL
  omega_J_M_output <- NULL
  omega_output <- NULL
  alpha_output <- 0
  alpha_zero_output <- 0
  allocation_prob_output <- NULL
  q_star_1_J_output <- NULL
  gamma_star_1_J_output <- NULL

  #--------------------------------------------------------------------------------

  # Acceptance probability
  acceptance_prob_list <- data.frame(omega = rep(0, number_iter-1),
                                     alpha = rep(0, number_iter-1),
                                     alpha_zero = rep(0, number_iter-1),
                                     q_star = rep(0, number_iter-1),
                                     gamma_star = rep(0, number_iter-1))


  # ----------------------- Step 1: Set starting values ----------------------------
  if(!run_post){
    # full algorithm
    if (is.null(Z.init)){
      # Initialize similar to k-means++
      data.j.prior <- matrix(do.call(cbind, Y), nrow = R)
      pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
      Ycbind = do.call(cbind, Y)
      q_temp = matrix(0,1,R)
      c = sample(1:sum(C), size=1)
      q_temp[1,] = extraDistr::rdirichlet(1,alpha = Ycbind[,c] + 0.1 )
      ps = 1000
      q_prior = extraDistr::rdirichlet(ps, alpha = rep(1/J,J))
      for(j in 2:J){
        q_dist = sapply(1:sum(C), function(c){
          min(sapply(1:dim(q_temp)[1], function(l){
            c(dist(rbind(pdata.j.prior[,c],q_temp[l,]),method = "manhattan"))
          }))
        })
        temp = 10
        if(max(q_dist)<tau){
          q_dist_prior = sapply(1:ps, function(c){
            min(sapply(1:dim(q_temp)[1], function(l){
              c(dist(rbind(q_prior[c,],q_temp[l,]),method = "manhattan"))
            }))
          })
          myprob = exp(temp*q_dist_prior)*(q_dist_prior>tau)
          myprob = myprob/sum(myprob)
          c = sample(1:ps, size =1 , prob = myprob)
          q_temp = rbind(q_temp, q_prior[c,])
        }else{
          myprob = exp(temp*q_dist)*(q_dist>tau)
          myprob = myprob/sum(myprob)
          c =sample(1:sum(C), size=1, prob = myprob)
          q_temp = rbind(q_temp, extraDistr::rdirichlet(1,alpha = Ycbind[,c] + 0.1 ))
        }
        mindist = min(dist(q_temp, method = "manhattan"))
        if(mindist<tau){
          print(paste('min dist between initial q is ', mindist))
        }
      }
      # Allocate z based on closeness to q
      Z_new = lapply(1:M, function(m){
        Z.new.m = sapply(1:C[m], function(c){
          which.min(sapply(1:dim(q_temp)[1], function(l){
            c(dist(rbind(pdata.j.prior[,c],q_temp[l,]),method = "manhattan"))
          }))
        })
      })
    }else{
      Z_new <- Z.init
    }

  }else{
    Z_new <- Z.fix
  }

  if((!run_post) || (run_post & run_omega)){
    # Empirical values based on counts
    omega_hat = function(z=z,J=J){
      omega = rep(1,J)
      tbz = table(z)
      omega[as.numeric(names(tbz))] = omega[as.numeric(names(tbz))] + tbz
      return(omega/sum(omega))
    }

    omega_J_M_new <- matrix(unlist(lapply(Z_new, omega_hat, J=J)), J, M, byrow=FALSE)

    omega_new <- omega_hat(unlist(Z_new), J=J)

    # Empirical values based on DP approximation
    alpha_zero_new  <- mean(unlist(lapply(Z_new,function(x){length(unique(x))}))/log(C))
    alpha_new <-length(unique(unlist(Z_new)))/log(sum(unlist(lapply(Z_new,function(x){length(unique(x))}))))
  }


  # Hyperparameters for q
  if((!run_post) || (run_post & run_q_gamma)){

    # Compute empirical mean of projection strength
    data.j.prior <- matrix(do.call(cbind, Y), nrow = R)
    pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
    data.j.mean <- rowMeans(pdata.j.prior)
    mx <- matrix(data.j.mean/sum(data.j.mean), nrow = 1)
    alpha_h <- a*mx

    # Set q empirically based on Z_new
    q_star_1_J_new <- lapply(1:J,
                             function(j){

                               if(length(which(unlist(Z_new)==j)) != 0){

                                 data.j.prior <- matrix(do.call(cbind, Y)[,which(unlist(Z_new)==j)],
                                                        nrow = R)

                                 pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
                                 data.j.mean <- rowMeans(pdata.j.prior)

                                 mx2 <- matrix(data.j.mean/sum(data.j.mean),
                                               nrow = 1)

                               }else{

                                 mx2 <- matrix(0, nrow = 1, ncol = R)
                                 # mx <- matrix(rdirichlet(1, alpha_h), nrow = 1, ncol = R)
                               }

                               mx2

                             })

    q_star_1_J_new <- do.call(rbind, q_star_1_J_new)

    # mindist = min(dist(q_star_1_J_new, method = "manhattan"))
    # if(mindist<tau){
    #   print(paste('min dist between q is ', mindist))
    # }

    ind = c(1:J)[sapply(1:J, function(j){length(which(unlist(Z_new)==j)) == 0})]
    ind2 = c(1:J)[sapply(1:J, function(j){length(which(unlist(Z_new)==j)) > 0})]
    q_temp = q_star_1_J_new[ind2,]
    # Initialize empty ones that are far from occupied ones
    Ycbind = do.call(cbind, Y)
    for(j in ind){
      # Find the data point further from the current q
      q_dist = sapply(1:sum(C), function(c){
        min(sapply(1:dim(q_temp)[1], function(l){
          c(dist(rbind(pdata.j.prior[,c],q_temp[l,]),method = "manhattan"))
        }))
      })
      c =which.max(q_dist)
      q_star_1_J_new[j,] = extraDistr::rdirichlet(1,alpha = Ycbind[,c] + 0.1 )
      q_temp = rbind(q_temp,q_star_1_J_new[j,])
    }
    q_star_1_J_new <- ifelse(q_star_1_J_new == 0, .Machine$double.eps,  q_star_1_J_new)
    q_star_1_J_new <-  q_star_1_J_new/rowSums(q_star_1_J_new)

    # Set gamma empirically based on Z_new
    gamma_1_J_star_new  <- lapply(1:J,
                                  function(j){

                                    if(length(which(unlist(Z_new)==j)) > 1){

                                      data.j.prior <- matrix(do.call(cbind, Y)[,which(unlist(Z_new)==j)],
                                                             nrow = R)

                                      pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
                                      data.j.var <- apply(pdata.j.prior,1,var)
                                      data.j.mean <- rowMeans(pdata.j.prior)

                                      mx2 <- mean(data.j.mean*(1-data.j.mean)/(data.j.var+ 0.0001)  - 1)
                                      mx2 <- mx2*(mx2>0) + (mx2<=0)
                                      mx2 = max(mx2, lb_gamma +1 )

                                    }else{

                                      mx2 <- cascsim::rtgamma(1,shape = a_gamma, scale = 1/b_gamma, min = lb_gamma)
                                    }

                                    mx2

                                  })

    gamma_1_J_star_new <- do.call(rbind, gamma_1_J_star_new)
  }

  # Total count of acceptance
  omega_count <- 0
  alpha_count <- 0
  alpha_zero_count <- 0
  q_star_count <- 0
  gamma_star_count <- 0
  split_accept = 0
  merge_accept = 0

  #----------------------- Step 2: Prepare for the covariance update -----------------------------

  # omega
  if((!run_post) || (run_post & run_omega)){

    mean_X_component_new <- matrix(0,
                                   nrow = 1,
                                   ncol = J-1)
    tilde_s_component_new <- matrix(0,
                                    nrow = J-1,
                                    ncol = J-1)
    covariance_component_new <- matrix(0,
                                       nrow = J-1,
                                       ncol = J-1)
    scaling_omega = 0.01

    # alpha
    variance_alpha_new <- 0.1

    # alpha_zero
    variance_alpha_zero_new <- 0.1
  }

  # q star
  if((!run_post) || (run_post & run_q_gamma)){

    mean_x_q_new <- lapply(1:J, function(j) matrix(0,
                                                   nrow = 1,
                                                   ncol = R-1))
    tilde_s_q_new <- lapply(1:J, function(j) matrix(0,
                                                    nrow = R-1,
                                                    ncol = R-1))
    covariance_q_new <- lapply(1:J, function(j) matrix(0, nrow = R-1, ncol = R-1))
    q_star_iter_occupied = rep(1,J)
    scaling_q = rep(0.1, J)

    # gamma star
    variance_gamma_new <- rep(0.1,J)
    gamma_star_iter_occupied = rep(1,J)
  }
  # ----------------------- Step 3: Updates -----------------------------

  ##-- Set initial index
  output_index <- 0

  start_time_mcmc <- Sys.time()
  if(verbose) {
    print(paste('Start MCMC:', start_time_mcmc))
    cat('\n')

    pb <- txtProgressBar(min = 1, max = number_iter, style = 3)
  }

  # Iteration starts with number 2
  for(iter in 2:number_iter){

    if(verbose) {setTxtProgressBar(pb, iter)}

    if(iter >= burn_in & iter%%thinning == 0){

      output_index <- output_index + 1
      update <- TRUE
    }else{

      update <- FALSE
    }


    ##----------------------------------- Allocation variables ------------------------------------

    if(!run_post){
      allocation_output <- allocation_variables_dirmult_mcmc_cpp(omega_J_M = omega_J_M_new,
                                                                 q_star_1_J = q_star_1_J_new,
                                                                 gamma_1_J_star = gamma_1_J_star_new,
                                                                 Y = Y)

      Z_new <- allocation_output$Z
      allocation.prob <- allocation_output$allocation.prob
    }


    ##--------------------------- q_star ----------------------------------

    if((!run_post) || (run_post & run_q_gamma)){
      q_output_sim <- q_star_mcmc_cpp(Y = Y,
                                      Z = Z_new,
                                      q_star_1_J = q_star_1_J_new,
                                      gamma_1_J_star = gamma_1_J_star_new,
                                      alpha_h = alpha_h,
                                      tau = tau,
                                      nu = nu,
                                      covariance = covariance_q_new,
                                      mean_x = mean_x_q_new,
                                      tilde_s = tilde_s_q_new,
                                      scaling = scaling_q,
                                      iter_num = q_star_iter_occupied,
                                      adaptive_prop = adaptive_prop)

      covariance_q_new <- q_output_sim$covariance_new
      mean_x_q_new <- q_output_sim$mean_x_new
      tilde_s_q_new <- q_output_sim$tilde_s_new
      scaling_q <- q_output_sim$scaling
      q_star_1_J_new <- q_output_sim$q_star_1_J_new
      q_star_count <- q_output_sim$q_star_count + q_star_count
      acceptance_prob_list$q_star[iter-1] <- q_star_count/((iter-1)*J)
      q_star_iter_occupied = q_output_sim$iter_occupied

      # mindist = min(dist(q_star_1_J_new, method = "manhattan"))
      # if(mindist<tau){
      #   print(paste('min dist between q is ', mindist))
      # }


    ##-------------------------- gamma ---------------------------------

    # gamma_output_sim <- gamma_mcmc(Y = Y,
    #                                Z = Z_new,
    #                                q_star_1_J = q_star_1_J_new,
    #                                gamma_1_J_star = gamma_1_J_star_new,
    #                                a_gamma = a_gamma,
    #                                b_gamma = b_gamma,
    #                                X_mean = X_mean_gamma_new,
    #                                M_2 = M_2_gamma_new,
    #                                variance = variance_gamma_new,
    #                                iter_num = iter,
    #                                adaptive_prop = adaptive_prop)
    #
    # gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
    # variance_gamma_new <- gamma_output_sim$variance_gamma_new
    # M_2_gamma_new <- gamma_output_sim$M_2_gamma_new
    # X_mean_gamma_new <- gamma_output_sim$X_mean_gamma_new
    # gamma_star_count <- gamma_star_count + gamma_output_sim$gamma_count
    # acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)

      gamma_output_sim <- gamma_mcmc_cpp(Y = Y,
                                         Z = Z_new,
                                         q_star_1_J = q_star_1_J_new,
                                         gamma_1_J_star = gamma_1_J_star_new,
                                         a_gamma = a_gamma,
                                         b_gamma = b_gamma,
                                         lb_gamma = lb_gamma,
                                         variance = variance_gamma_new,
                                         iter_num = gamma_star_iter_occupied)

      gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
      variance_gamma_new <- gamma_output_sim$variance_gamma_new
      gamma_star_count <- gamma_star_count + gamma_output_sim$gamma_count
      acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)
      gamma_star_iter_occupied = gamma_output_sim$iter_occupied


    }
    ##------------------------ Dataset specific component probabilities ------------------------------

    if((!run_post) || (run_post & run_omega)){
      dataset_specfic_output <- dataset_specific_mcmc(Z = Z_new,
                                                      omega = omega_new,
                                                      alpha = alpha_new)

      omega_J_M_new <- dataset_specfic_output




    ##------------------------ Component probabilities -----------------------------------------------

      component_output <- component_probabilities_mcmc(omega = omega_new,
                                                       omega_J_M = omega_J_M_new,
                                                       alpha_zero = alpha_zero_new,
                                                       alpha = alpha_new,
                                                       covariance = covariance_component_new,
                                                       mean_x = mean_X_component_new,
                                                       tilde_s = tilde_s_component_new,
                                                       scaling = scaling_omega,
                                                       iter_num = iter,
                                                       adaptive_prop = adaptive_prop)

      omega_new <- component_output$omega_new
      tilde_s_component_new <- component_output$tilde_s_new
      mean_X_component_new <- component_output$mean_x_new
      scaling_omega <- component_output$scaling
      covariance_component_new <- component_output$covariance_new
      omega_count <- omega_count + component_output$accept
      acceptance_prob_list$omega[iter-1] <- omega_count/(iter-1)

    ##-------------------------- Update alpha ----------------------------------

      alpha_output_sim <- alpha_mcmc(omega_J_M = omega_J_M_new,
                                     omega = omega_new,
                                     alpha = alpha_new,
                                     a_alpha = a_alpha,
                                     b_alpha = b_alpha,
                                     variance = variance_alpha_new,
                                     iter_num = iter)

      alpha_new <- alpha_output_sim$alpha_new
      variance_alpha_new <- alpha_output_sim$variance_new
      alpha_count <- alpha_count + alpha_output_sim$accept
      acceptance_prob_list$alpha[iter-1] <- alpha_count/(iter-1)

    ##----------------------- Update alpha_zero --------------------------------

      alpha_zero_output_sim <- alpha_zero_mcmc(omega = omega_new,
                                               alpha_zero = alpha_zero_new,
                                               a_alpha0 = a_alpha0,
                                               b_alpha0 = b_alpha0,
                                               variance = variance_alpha_zero_new,
                                               iter_num = iter)

      alpha_zero_new <- alpha_zero_output_sim$alpha_zero_new
      variance_alpha_zero_new <- alpha_zero_output_sim$variance_new
      alpha_zero_count <- alpha_zero_count + alpha_zero_output_sim$accept
      acceptance_prob_list$alpha_zero[iter-1] <- alpha_zero_count/(iter-1)

    }

    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){

      if(!run_post){
        Z_output[[output_index]] <- Z_new
        allocation_prob_output[[output_index]] <- allocation.prob
      }
      if((!run_post) || (run_post & run_omega)){
        omega_J_M_output[[output_index]] <- omega_J_M_new
        omega_output[[output_index]] <- omega_new
        alpha_output[output_index] <- alpha_new
        alpha_zero_output[output_index] <- alpha_zero_new
      }
      if((!run_post) || (run_post & run_q_gamma)){
        q_star_1_J_output[[output_index]] <- q_star_1_J_new
        gamma_star_1_J_output[[output_index]] <- gamma_1_J_star_new
      }
    }

    if(auto_save == TRUE && iter %% save_frequency == 0){
      my_list <- list('acceptance_prob' = acceptance_prob_list,
                      'M' = M,
                      'C' = C,
                      'R' = R,
                      'J' = J,
                      'Z' = Z.fix,
                      'allocation_probability' = allocation_prob_output,
                      'Z_output' = Z_output,
                      'omega_J_M_output' = omega_J_M_output,
                      'omega_output' = omega_output,
                      'alpha_output' = alpha_output,
                      'alpha_zero_output' = alpha_zero_output,
                      'q_star_1_J_output' = q_star_1_J_output,
                      'gamma_star_1_J_output' = gamma_star_1_J_output,
                      'post' = run_post,
                      'output_index' = output_index)
      save(my_list, file=save_path)

    }


  }

  end_time <- Sys.time()
  diff_time <- difftime(end_time,start_time_mcmc)

  if(verbose) {
    close(pb)

    cat('\n')
    print(paste('End:',end_time))
    cat('\n')
    print(paste('MCMC running time:', round(diff_time, digits = 3),units(diff_time)))

  }


  my_list <- list('acceptance_prob' = acceptance_prob_list,
                  'M' = M,
                  'C' = C,
                  'R' = R,
                  'J' = J,
                  'Z' = Z.fix,
                  'allocation_probability' = allocation_prob_output,
                  'Z_output' = Z_output,
                  'omega_J_M_output' = omega_J_M_output,
                  'omega_output' = omega_output,
                  'alpha_output' = alpha_output,
                  'alpha_zero_output' = alpha_zero_output,
                  'q_star_1_J_output' = q_star_1_J_output,
                  'gamma_star_1_J_output' = gamma_star_1_J_output,
                  'post' = run_post,
                  'output_index' = output_index




  )

  return(my_list)

}



