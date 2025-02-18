##--------------- Simulation of allocations ----------------
allocation_variables_dirmult_mcmc_cpp <- function(omega_J_M,
                                              q_star_1_J,
                                              gamma_1_J_star,
                                              Y){

  # Create an empty list
  M <- length(Y)
  C <- unlist(lapply(Y, ncol))
  J <- nrow(omega_J_M)
  R <- ncol(q_star_1_J)

  ## Initialize Z
  Z <- NULL
  allocation.prob <- NULL

  for(m in 1:M){

    Prob_m <- compute_probability_cpp(Y = Y[[m]], omega_M = omega_J_M[,m], q_star_1_J = q_star_1_J, gamma_1_J_star = gamma_1_J_star)

    ##-- Store Z
    # Z[[m]] <- apply(Prob_m, 1, function(x) extraDistr::rcat(n = 1, prob = x))
    Z[[m]] <- Z_sample_cpp(Prob_m)

    ##-- Store allocation probability
    allocation.prob[[m]] <- Prob_m[cbind(1:C[m], Z[[m]])]

  }

  ##-- Return Z
  return(list('Z' = Z,
              'allocation.prob' = allocation.prob))
}


##-------------------------- q_star ------------------------------

q_replustion = function(q_star,
                        tau,
                        nu){

  # Compute distance between all q_star
  delta.vec = c(stats::dist(q_star, method = "manhattan"))
  # Compute transformation on the log scale
  g_q = - log(1+exp(-(delta.vec -tau)/nu))
  # Compute resulsion term
  h_q = min(g_q)
  aux = -(min(delta.vec) -tau)/nu
  # h_q = -aux - log(exp(-aux) + 1)
  h_q = -aux
  return(h_q)
}


q_star_mcmc_cpp <- function(Y,
                        Z,
                        q_star_1_J,
                        gamma_1_J_star,
                        alpha_h,
                        tau,
                        nu,
                        covariance,
                        mean_x,
                        tilde_s,
                        scaling,
                        iter_num,
                        adaptive_prop,
                        lambda = 0.7){

  J <- nrow(q_star_1_J)
  R <- ncol(q_star_1_J)
  n <- iter_num

  # Define the inputs
  q_star_1_J_old <- q_star_1_J
  covariance_old <- covariance
  mean_X_old <- mean_x
  tilde_s_old <- tilde_s
  X_old <- matrix(0, nrow = J, ncol = R-1) # Dimension = J x R-1
  for(j in 1:J){

    X_old[j,] <- log(q_star_1_J_old[j,1:(R-1)]/q_star_1_J_old[j,R])
  }

  q_star_1_J_new <- q_star_1_J
  tilde_s_new <- tilde_s
  mean_x_new <- mean_x
  covariance_new <- covariance
  q_star_count <- 0

  # For each cluster
  for(j in 1:J){


    if(length(which(unlist(Z)==j)) != 0){

      if(n[j] <= 200){

        # Length of R - 1
        X_j_new <- as.vector(mvtnorm::rmvnorm(n = 1,
                                     mean = X_old[j,],
                                     sigma = diag(scaling[j],nrow = R-1, ncol = R-1)))

      }else{

        # Length of R - 1
        X_j_new <- as.vector(mvtnorm::rmvnorm(n = 1,
                                     mean = X_old[j,],
                                     sigma = scaling[j]*(covariance_old[[j]] + adaptive_prop*diag(1,nrow = R-1, ncol = R-1))))

      }

      # q_star_new
      aux = max(X_j_new)
      q_star_j_new <- c(exp(X_j_new-aux)/(exp(-aux) + sum(exp(X_j_new-aux))),exp(-aux)/(exp(-aux)+sum(exp(X_j_new-aux))))
      if(is.na(q_star_j_new[R])){q_star_j_new[R] = 1 - 10^(-15)}
      # shifted the zero values
      q_star_j_new <- ifelse( q_star_j_new == 0, .Machine$double.eps,  q_star_j_new)
      q_star_j_new <-  q_star_j_new/sum(q_star_j_new)
      q_star_1_J_new_update = q_star_1_J_new
      q_star_1_J_new_update[j,] = q_star_j_new

      # Acceptance probability
      log.accept.prob <- q_star_logprob_cpp(Y_sub = as.matrix((do.call(cbind, Y))[,unlist(Z)==j]),
                                        q_j_star = q_star_j_new,
                                        gamma_j_star = gamma_1_J_star[j],
                                        alpha_h = alpha_h) +
        q_replustion(q_star_1_J_new_update,
                     tau,
                     nu) -

        q_star_logprob_cpp(Y_sub = as.matrix((do.call(cbind, Y))[,unlist(Z)==j]),
                       q_j_star = q_star_1_J_old[j,],
                       gamma_j_star = gamma_1_J_star[j],
                       alpha_h = alpha_h) -

        q_replustion(q_star_1_J_new,
                     tau,
                     nu) +

        sum(log(q_star_j_new)) - sum(log(q_star_1_J_old[j,]))


      accept.prob <- exp(log.accept.prob)

      if(n[j]>100){
        # Update Covariance Structure
        n_cov = n[j]-100
        mean_x_new[[j]] <- mean_X_old[[j]]*(1-1/n_cov) + 1/n_cov*matrix(X_j_new, nrow = 1)
        tilde_s_new[[j]] <- tilde_s_old[[j]] + matrix(X_j_new, ncol = 1)%*%matrix(X_j_new, nrow = 1)
        f_cov = 1
        if(n_cov>1){ f_cov = 1/(n_cov-1)}
        covariance_new[[j]] <- f_cov*tilde_s_new[[j]] - n_cov*f_cov*t(mean_x_new[[j]])%*%mean_x_new[[j]]
      }
      # Update Scaling factor
      scaling[j] = scaling[j] + n[j]^(-lambda)*(min(1, exp(log.accept.prob))- 0.234)
      if(is.na(scaling[j])){
        print('nan issues')}
      if(scaling[j] < 0.0001){scaling[j] <- 0.0001}
      if(scaling[j] > 10){scaling[j] <- 10}
      n[j] = n[j] +1

    }else{

      # If component j is empty

      # Simulate from the prior
      q_star_j_new <- as.vector(extraDistr::rdirichlet(n = 1,
                                                       alpha = alpha_h))
      # shifted the zero values
      q_star_j_new <- ifelse( q_star_j_new == 0, .Machine$double.eps,  q_star_j_new)
      q_star_j_new <-  q_star_j_new/sum( q_star_j_new)
      q_star_1_J_new_update = q_star_1_J_new
      q_star_1_J_new_update[j,] = q_star_j_new

      # X_new
      X_j_new <- log(q_star_j_new[1:(R-1)]/q_star_j_new[R])

      # Acceptance probability
      log.accept.prob <- q_replustion(q_star_1_J_new_update,
                                      tau,
                                      nu) -

        q_replustion(q_star_1_J_new,
                     tau,
                     nu)

      accept.prob <- exp(log.accept.prob)

    }

    if(!is.finite(accept.prob)){
      print('issues with inf/nan')
    }


    # Random Bernoulli
    outcome <- rbinom(n = 1,
                      size = 1,
                      prob = min(1,accept.prob))

    # If outcome is to reject
    if(outcome == 0){

      q_star_j_new <- q_star_1_J_old[j,]
      X_j_new <- X_old[j,]

    }

    # # Update
    # tilde_s_new[[j]] <- tilde_s_old[[j]] + matrix(X_j_new, ncol = 1)%*%matrix(X_j_new, nrow = 1)
    # mean_x_new[[j]] <- mean_X_old[[j]]*(1-1/n) + 1/n*matrix(X_j_new, nrow = 1)
    # covariance_new[[j]] <- 1/(n-1)*tilde_s_new[[j]] - n/(n-1)*t(mean_x_new[[j]])%*%mean_x_new[[j]]
    #
    q_star_count <- q_star_count + outcome
    q_star_1_J_new[j,] <- q_star_j_new
  }



  # Return
  return(list('tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'scaling' = scaling,
              'q_star_count' = q_star_count,
              'q_star_1_J_new' = q_star_1_J_new,
              'iter_occupied' = n))


}



##------------------------------ Gamma ----------------------------------

gamma_mcmc_cpp <- function(Y,
                       Z,
                       q_star_1_J,
                       gamma_1_J_star,
                       a_gamma,
                       b_gamma,
                       lb_gamma,
                       variance,
                       iter_num,
                       lambda =0.7){

  J <- nrow(q_star_1_J)
  n <- iter_num
  R <- ncol(q_star_1_J)

  # previous values
  gamma_1_J_star_old <- gamma_1_J_star
  X_old <- log(gamma_1_J_star_old - lb_gamma)
  variance_old <- variance
  variance_new <- variance

  # Update
  gamma_count <- 0
  gamma_1_J_star_new <- NULL

  # For each component j
  for(j in 1:J){

    if(length(which(unlist(Z)==j)) != 0){

      X_j_new <- rnorm(n = 1,
                       mean = X_old[j],
                       sd = sqrt(variance_old[j]))

      # gamma_j_new
      gamma_j_new <- exp(X_j_new) + lb_gamma

      # Acceptance probability
      log.accept.prob <- gamma_logprob_cpp(Y_sub = as.matrix((do.call(cbind, Y))[,unlist(Z)==j]),
                                       q_j_star = q_star_1_J[j,],
                                       gamma_j_star = gamma_j_new,
                                       a_gamma = a_gamma,
                                       b_gamma = b_gamma,
                                       lb_gamma = lb_gamma)-

        gamma_logprob_cpp(Y_sub = as.matrix((do.call(cbind, Y))[,unlist(Z)==j]),
                      q_j_star = q_star_1_J[j,],
                      gamma_j_star = gamma_1_J_star_old[j],
                      a_gamma = a_gamma,
                      b_gamma = b_gamma,
                      lb_gamma = lb_gamma)-

        log(gamma_1_J_star_old[j]- lb_gamma) + log(gamma_j_new - lb_gamma)

      accept.prob <- exp(log.accept.prob)

      # Update Covariance Structure
      variance_new[j] <- variance_old[j] + n[j]^(-lambda)*(min(1,accept.prob)- 0.44)
      if(variance_new[j] < 0.0001){variance_new[j] <- 0.0001}
      if(variance_new[j] > 1000){variance_new[j] <- 1000 }
      n[j] = n[j] +1

    }else{

      # If component j is empty

      # Simulate from the prior
      gamma_j_new <- cascsim::rtgamma(n = 1,
                            shape = a_gamma,
                            scale = 1/b_gamma,
                            min = lb_gamma)

      # X_new
      X_j_new <- log(gamma_j_new - lb_gamma)

      # Always accept
      accept.prob <- 1

    }

    if(!is.finite(accept.prob)){
      print('issues with inf/nan')
    }

    # Random Bernoulli
    outcome <- rbinom(n = 1,
                      size = 1,
                      prob = min(1,accept.prob))

    # If outcome is to reject
    if(outcome == 0){

      gamma_j_new <- gamma_1_J_star_old[j]
      X_j_new <- X_old[j]
    }

    gamma_count <- gamma_count + outcome
    gamma_1_J_star_new[j] <- gamma_j_new

  }



  # Return
  return(list('gamma_1_J_star_new' = gamma_1_J_star_new,
              'variance_gamma_new' = variance_new,
              'gamma_count' = gamma_count,
              'iter_occupied' = n))

}


##------------ Simulation of Data-Specific component probabilities --------------------
dataset_specific_mcmc <- function(Z,
                                  omega,
                                  alpha){

  J <- length(omega)
  M <- length(Z)

  loop.result <- lapply(1:M, function(m) {

    parameter <- as.vector(table(factor(Z[[m]], levels = c(1:J)))) + alpha*omega

    omega_J_M <- extraDistr::rdirichlet(n = 1,
                                        alpha = parameter)


    # shifted the zero values
    omega_J_M <- ifelse(omega_J_M == 0, .Machine$double.eps, omega_J_M)
    omega_J_M <- omega_J_M/rowSums(omega_J_M)

    return(omega_J_M)
  })

  omega_J_M <- matrix(unlist(loop.result), nrow=J, ncol=M)

  # Returning omega_J_M
  return(omega_J_M)
}



##----------- Simulation of Component probabilities ----------------
component_log_prob <- function(omega,
                               omega_J_M,
                               alpha_zero,
                               alpha){

  J <- nrow(omega_J_M)
  M <- ncol(omega_J_M)

  lprod <- sum((alpha_zero/J-1)*log(omega)) +

    sum(unlist(lapply(1:M, function(m){

      sum(alpha*omega*log(omega_J_M[,m])-lgamma(alpha*omega))
    })))

  # outputs
  return(lprod)
}

component_probabilities_mcmc <- function(omega,
                                         omega_J_M,
                                         alpha_zero,
                                         alpha,
                                         covariance,
                                         mean_x,
                                         tilde_s,
                                         scaling,
                                         iter_num,
                                         adaptive_prop = 0.01,
                                         lambda = 0.7){

  J <- nrow(omega_J_M)
  M <- ncol(omega_J_M)

  # Define the inputs
  omega_old <- omega
  covariance_old <- covariance
  mean_X_d_old <- mean_x
  tilde_s_old <- tilde_s

  X_old <- log(omega_old[1:J-1]/omega_old[J]) # Length = J-1

  # Specify the iteration number
  n <- iter_num

  # Adaptive step
  if(n <= 200){

    X_new <- mvtnorm::rmvnorm(n = 1,
                     mean = X_old,
                     sigma = diag(x=scaling, nrow=J-1, ncol=J-1))

  }else{

    X_new <- mvtnorm::rmvnorm(n = 1,
                     mean = X_old,
                     sigma = scaling*(covariance_old + adaptive_prop*diag(1, nrow = J-1, ncol = J-1)))

  }

  # Compute omega_new (Length = J) from X_new
  omega_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))

  # shifted the zero values
  omega_new <- ifelse(omega_new == 0, .Machine$double.eps, omega_new)
  omega_new <- omega_new/sum(omega_new)

  # Compute acceptance probability
  log_acceptance <- component_log_prob(omega_new,
                                       omega_J_M,
                                       alpha_zero,
                                       alpha) -
    component_log_prob(omega_old,
                       omega_J_M,
                       alpha_zero,
                       alpha) +

    sum(log(omega_new)) - sum(log(omega_old))

  if(!is.finite(log_acceptance)){
    print('issues with inf/nan')
  }

  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob=min(1, exp(log_acceptance)))

  if(outcome == 0){

    X_new <- X_old
    omega_new <- omega_old
  }

  # # Update Covariance Structure
  # tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
  # mean_x_new <- mean_X_d_old*(1-1/n) + 1/n*matrix(X_new, nrow = 1)
  # covariance_new <- 1/(n-1)*tilde_s_new - n/(n-1)*t(mean_x_new)%*%mean_x_new

  # Update Covariance Structure
  tilde_s_new <- tilde_s_old
  mean_x_new <- mean_X_d_old
  covariance_new <- covariance_old

  if(n>100){
    # Update Covariance Structure
    n_cov = n-100
    tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
    mean_x_new <- mean_X_d_old*(1-1/n_cov) + 1/n_cov*matrix(X_new, nrow = 1)
    f_cov = 1
    if(n_cov>1){ f_cov = 1/(n_cov-1)}
    covariance_new <- f_cov*tilde_s_new - n_cov*f_cov*t(mean_x_new)%*%mean_x_new

  }

  # Update Scaling factor
  scaling = scaling + n^(-lambda)*(min(1, exp(log_acceptance))- 0.234)
  if(scaling < 0.0001){scaling <- 0.0001}
  if(scaling > 10){scaling <- 10 }



  # Output
  return(list('omega_new' = omega_new,
              'tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'scaling' = scaling,
              'covariance_new' = covariance_new,
              'accept' = outcome))
}


## ----------------------- Simulation of Concentration parameter alpha -----------------------
alpha_log_prob <- function(omega_J_M,
                           omega,
                           alpha,
                           a_alpha,
                           b_alpha){

  M <- ncol(omega_J_M)

  lprod <- -alpha*b_alpha + M*lgamma(alpha) + (a_alpha-1)*log(alpha) +

    sum(unlist(lapply(1:M, function(m) {
      sum(alpha*omega*log(omega_J_M[,m])-lgamma(alpha*omega))
    })))

  # output
  return(lprod)
}

alpha_mcmc <- function(omega_J_M,
                       omega,
                       alpha,
                       a_alpha,
                       b_alpha,
                       variance,
                       iter_num,
                       lambda = 0.7){

  M <- ncol(omega_J_M)

  # Defining the inputs
  alpha_old <- alpha
  X_old <- log(alpha_old)
  variance_old <- variance

  # Defining the dimensions
  n <- iter_num

  X_new <- rnorm(n = 1,
                 mean = X_old,
                 sd = sqrt(variance_old))

  # alpha_new
  alpha_new <- exp(X_new)

  # Compute log acceptance probability
  log_acceptance <- alpha_log_prob(omega_J_M = omega_J_M,
                                   omega = omega,
                                   alpha = alpha_new,
                                   a_alpha =  a_alpha,
                                   b_alpha =  b_alpha) -

    alpha_log_prob(omega_J_M = omega_J_M,
                   omega = omega,
                   alpha = alpha_old,
                   a_alpha =  a_alpha,
                   b_alpha =  b_alpha) +

    log(alpha_new) - log(alpha_old)

  if(!is.finite(log_acceptance)){
    print('issues with inf/nan')
  }

  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1, exp(log_acceptance)))


  if(outcome == 0){

    X_new <- X_old
    alpha_new <- alpha_old
  }

  # # Update Covariance Structure
  variance_new <- variance_old + n^(-lambda)*(min(1, exp(log_acceptance))- 0.44)
  if(variance_new < 0.0001){variance_new <- 0.0001}
  if(variance_new > 1000){variance_new <- 1000 }

  # Output
  return(list('alpha_new' = alpha_new,
              'variance_new' = variance_new,
              'accept' = outcome))
}

##---------------------- Simulation of concentration parameter alpha_zero -----------------
alpha_zero_log_prob <- function(omega,
                                alpha_zero,
                                a_alpha0,
                                b_alpha0){

  # dimension
  J <- length(omega)

  # log probability
  lprob <- -b_alpha0*alpha_zero + lgamma(alpha_zero) - J*lgamma(alpha_zero/J) + sum(alpha_zero/J*log(omega)) + (a_alpha0-1)*log(alpha_zero)

  # output
  return(lprob)
}

alpha_zero_mcmc <- function(omega,
                            alpha_zero,
                            a_alpha0,
                            b_alpha0,
                            variance,
                            iter_num,
                            lambda = 0.7){

  # dimension
  J <- length(omega)
  n <- iter_num

  # previous values
  alpha_zero_old <- alpha_zero
  X_old <- log(alpha_zero_old)
  variance_old <- variance

  # Simulation
  X_new <- rnorm(n = 1,
                 mean = X_old,
                 sd = sqrt(variance_old))

  # alpha_zero_new
  alpha_zero_new <- exp(X_new)

  # log acceptance probability
  log_acceptance <- alpha_zero_log_prob(omega = omega,
                                        alpha_zero = alpha_zero_new,
                                        a_alpha0 = a_alpha0,
                                        b_alpha0 = b_alpha0) -

    alpha_zero_log_prob(omega = omega,
                        alpha_zero = alpha_zero_old,
                        a_alpha0 = a_alpha0,
                        b_alpha0 = b_alpha0) +

    log(alpha_zero_new) - log(alpha_zero_old)

  if(!is.finite(log_acceptance)){
    print('issues with inf/nan')
  }

  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1, exp(log_acceptance)))



  if(outcome == 0){

    X_new <- X_old
    alpha_zero_new <- alpha_zero_old
  }

  # # Update Covariance Structure
  variance_new <- variance_old + n^(-lambda)*(min(1, exp(log_acceptance))- 0.44)
  if(variance_new < 0.0001){variance_new <- 0.0001}
  if(variance_new > 10){variance_new <- 10 }


  # Output
  return(list('alpha_zero_new' = alpha_zero_new,
              'variance_new' = variance_new,
              'accept' = outcome))
}
