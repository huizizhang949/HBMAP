library(microbenchmark)

# ------ test z full conditional probs --------
# use break point inside the code to reach this step
Prob_m <- compute_probability_cpp(Y = Y[[m]], omega_M = omega_J_M[,m], q_star_1_J = q_star_1_J, gamma_1_J_star = gamma_1_J_star)


Prob1 <- t(sapply(1:C[m],
                       function(c){
                         
                         
                         ##-- log_p_tilde
                         log.p.tilde <- sapply(1:J, function(j){
                           
                           ddirmultinomial(Y[[m]][,c],
                                           n = sum(Y[[m]][,c]),
                                           alpha = q_star_1_J[j,]*gamma_1_J_star[j],
                                           log = TRUE)+
                             log(omega_J_M[j,m])
                         })
                         
                         ##-- logK
                         logK <- -max(log.p.tilde)

                         ##-- probability
                         allo.prob <- exp(log.p.tilde + logK)/sum(exp(log.p.tilde + logK))

                         allo.prob
                       }))

all.equal(c(Prob_m),c(Prob1))

# all.equal(Prob[,1],colSums(Y[[m]]))

a <- sapply(1:C[m],function(cc) Prob_m[cc,Z[[m]]][cc])
identical(a,Prob_m[cbind(1:C[m], Z[[m]])])


# ------- test sampling Z ----------------
set.seed(4)
# draw_cat(c(0.8,0.1,0.1))

k=5
set.seed(6)
table(replicate(1000,draw_cat(rep(1/5,5))))
draw_cat(rdirichlet(1,rep(1/5,5)))

# find_label(c(0,0,0,1,0,0,0,0))

load('huizi/Prob_m.RData')

# compare time
microbenchmark({
  apply(Prob_m, 1, function(x) extraDistr::rcat(n = 1, prob = x))
},
{
  Z_sample_cpp(Prob_m)
}
)


# ------- test ddirmultinomial_cpp --------
microbenchmark(ddirmultinomial(c(1,2,3)+100,6,rep(0.1,3),log = TRUE), 
               ddirmultinomial_cpp(c(1,2,3)+100,6,rep(0.1,3),log = TRUE),check = 'equivalent')

ddirmultinomial(c(1,2,3),6,rep(0.1,3),log = TRUE)
ddirmultinomial_cpp(c(1,2,3),6,rep(0.1,3),log = TRUE)


# ------ test q_star_logprob_cpp -----------
# time comparisons
# use initial values for evaluation of q_logprob
microbenchmark(
  {
    q_star_logprob(Y = data_Hans,
                   Z = clustinit,
                   q_j_star = rep(1/R,R),
                   gamma_j_star = 1,
                   alpha_h = rep(1/R,R),
                   j = 11)
    
  },
  # {
  #   
  #   q_star_logprob_cpp1(Y_bind = do.call(cbind,data_Hans),
  #                      Z_bind = unlist(clustinit),
  #                      q_j_star = rep(1/R,R),
  #                      gamma_j_star = 1,
  #                      alpha_h = rep(1/R,R),
  #                      j = 11)
  #   
  # },
  {
    q_star_logprob_cpp(Y_sub = (do.call(cbind, data_Hans))[,unlist(clustinit)==11],
                       rep(1/R,R),
                       gamma_j_star = 1,
                       alpha_h = rep(1/R,R))
    
  },
  check='equivalent'
)


# time comparisons
# use break point insde the code to reach this function
microbenchmark({
  q_star_logprob(Y = Y,
                 Z = Z,
                 q_j_star = q_star_j_new,
                 gamma_j_star = gamma_1_J_star[j],
                 alpha_h = alpha_h,
                 j = j)
},{
  q_star_logprob_cpp(Y_sub = (do.call(cbind, Y))[,unlist(Z)==j],
                     q_j_star = q_star_j_new,
                     gamma_j_star = gamma_1_J_star[j],
                     alpha_h = alpha_h)
  
  
},
check='equivalent'
)

q_star_logprob(Y = Y,
               Z = Z,
               q_j_star = q_star_j_new,
               gamma_j_star = gamma_1_J_star[j],
               alpha_h = alpha_h,
               j = j)

q_star_logprob_cpp(Y_bind = do.call(cbind,Y),
               Z_bind = unlist(Z),
               q_j_star = q_star_j_new,
               gamma_j_star = gamma_1_J_star[j],
               alpha_h = alpha_h,
               j = j)


q_star_logprob_cpp(Y_sub = (do.call(cbind, Y))[,unlist(Z)==j],
               q_j_star = q_star_j_new,
               gamma_j_star = gamma_1_J_star[j],
               alpha_h = alpha_h)

# ------ test gamma_logprob_cpp -----------
# time comparisons
# use initial values for evaluation of gamma_logprob
microbenchmark(
  {
    gamma_logprob(Y = data_Hans,
                  Z = clustinit,
                  q_j_star = rep(1/R,R),
                  gamma_j_star = 1,
                  a_gamma = a_gamma,
                  b_gamma = b_gamma,
                  j = 11)
    
  },
  # {
  #   
  #   q_star_logprob_cpp1(Y_bind = do.call(cbind,data_Hans),
  #                      Z_bind = unlist(clustinit),
  #                      q_j_star = rep(1/R,R),
  #                      gamma_j_star = 1,
  #                      alpha_h = rep(1/R,R),
  #                      j = 11)
  #   
  # },
  {
    gamma_logprob_cpp(Y_sub = (do.call(cbind, data_Hans))[,unlist(clustinit)==11],
                      rep(1/R,R),
                      gamma_j_star = 1,
                      a_gamma,
                      b_gamma)
    
  },
  check='equivalent'
)
gamma_logprob(Y = Y,
              Z = Z,
              q_j_star = q_star_1_J[j,],
              gamma_j_star = gamma_j_new,
              a_gamma = a_gamma,
              b_gamma = b_gamma,
              j = j)

gamma_logprob_cpp(Y_sub = (do.call(cbind, Y))[,unlist(Z)==j],
              q_j_star = q_star_1_J[j,],
              gamma_j_star = gamma_j_new,
              a_gamma = a_gamma,
              b_gamma = b_gamma)

# ------- compare overall run-time with sara's code (no split and merge, no print)--------
# time comparisons
microbenchmark({
  mcmc_all_hans <- mcmc_run_all_cpp(Y = data_Hans,
                                    J = J,
                                    number_iter = 10,
                                    thinning = 1,
                                    burn_in = 0,
                                    adaptive_prop = 0.0001,
                                    print_Z = FALSE,
                                    # iter_update = 100,
                                    a_gamma = a_gamma,
                                    b_gamma = b_gamma,
                                    a_alpha =  a_alpha,
                                    b_alpha= 1,
                                    a_alpha0 = a_alpha0,
                                    b_alpha0= 1,
                                    Z.init = clustinit)
  
},
{
  mcmc_all_hans1 <- mcmc_run_all(Y = data_Hans,
                                J = J,
                                number_iter = 10,
                                thinning = 1,
                                burn_in = 0,
                                adaptive_prop = 0.0001,
                                print_Z = FALSE,
                                # iter_update = 100,
                                a_gamma = a_gamma,
                                b_gamma = b_gamma,
                                a_alpha =  a_alpha,
                                b_alpha= 1,
                                a_alpha0 = a_alpha0,
                                b_alpha0= 1,
                                Z.init = clustinit)
  
},times = 10
)
mcmc_all_hans <- mcmc_run_all_cpp(Y = data_Hans,
                                  J = J,
                                  number_iter = 100,
                                  thinning = 1,
                                  burn_in = 0,
                                  adaptive_prop = 0.0001,
                                  print_Z = FALSE,
                                  # iter_update = 100,
                                  a_gamma = a_gamma,
                                  b_gamma = b_gamma,
                                  a_alpha =  a_alpha,
                                  b_alpha= 1,
                                  a_alpha0 = a_alpha0,
                                  b_alpha0= 1,
                                  Z.init = clustinit)

mcmc_all_hans <- mcmc_run_all(Y = data_Hans,
                                  J = J,
                                  number_iter = 100,
                                  thinning = 1,
                                  burn_in = 0,
                                  adaptive_prop = 0.0001,
                                  print_Z = FALSE,
                                  # iter_update = 100,
                                  a_gamma = a_gamma,
                                  b_gamma = b_gamma,
                                  a_alpha =  a_alpha,
                                  b_alpha= 1,
                                  a_alpha0 = a_alpha0,
                                  b_alpha0= 1,
                                  Z.init = clustinit)

