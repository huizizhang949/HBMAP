#' Posterior predictive checks based on a single replicate
#'
#' @description
#' This function generates a single replicate of the real data and compares empirical projection strengths for each mouse
#'
#'
#' @param mcmc_run_all_output output from \code{HBMAP_mcmc} for the full algorithm.
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param regions.name a character vector of region names.
#' @param facet_nrow number of rows in the figure. Default to 2.
#'
#' @return a list of ggplots. Each ggplot corresponds to a single mouse.
#' The plot compares the distribution of empirical projection strengths between the real data and the replicate.
#' @export
ppc_single <- function(mcmc_run_all_output,
                       Y,
                       regions.name = NULL,
                       facet_nrow = 2){



  length_per_chain <- length(mcmc_run_all_output$Z_output)

  C <- mcmc_run_all_output$C
  R <- mcmc_run_all_output$R
  M <- mcmc_run_all_output$M

  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }


  # Sum of counts for each neuron cell
  N_CM <- lapply(1:M,
                 function(m) colSums(Y[[m]]))

  # theta
  target.index <- sample(length_per_chain,
                         size = 1,
                         replace = FALSE)


  theta <- list(q_j_star = mcmc_run_all_output$q_star_1_J_output[[target.index]],
                gamma_j_star = mcmc_run_all_output$gamma_star_1_J_output[[target.index]],
                Z = mcmc_run_all_output$Z_output[[target.index]])

  # Simulate replicated data Y

  Y_rep <- lapply(1:M,
                  function(m){

                    do.call(cbind, lapply(1:C[m],
                                          function(c){

                                            prob.c <- extraDistr::rdirichlet(n = 1,
                                                                 theta$q_j_star[theta$Z[[m]][c],]*theta$gamma_j_star[theta$Z[[m]][c]])

                                            rmultinom(n = 1, size = N_CM[[m]][c], prob = prob.c)
                                          }))
                  })




  #------------------------------------------ Some statistics --------------------------------------

  Y_prop <- lapply(1:M,
                   function(m){

                     matrix(as.vector(Y[[m]])/rep(colSums(Y[[m]]), each = R),
                            nrow = R)
                   })

  Y_rep_prop <- lapply(1:M,
                       function(m){

                         matrix(as.vector(Y_rep[[m]])/rep(colSums(Y_rep[[m]]), each = R),
                                nrow = R)
                       })

  Y_prop_df <- data.frame(Region = rep(regions.name, sum(C)),
                          data = 'observed',
                          mouse = rep(paste('mouse', 1:M), c(R*C)),
                          value = as.vector(do.call(cbind,Y_prop)))

  Y_rep_prop_df <- data.frame(Region = rep(regions.name, sum(C)),
                              data = 'replicated',
                              mouse = rep(paste('mouse', 1:M), c(R*C)),
                              value = as.vector(do.call(cbind,Y_rep_prop)))

  # For each of the mouse, draw histograms
  plot.out <- lapply(1:M,
                     function(m){


                       plot.out.m <- rbind(Y_prop_df,
                                           Y_rep_prop_df) %>%
                         dplyr::filter(mouse == paste('mouse', m)) %>%
                         ggplot(mapping = aes(x = value,
                                              fill = data))+
                         geom_histogram(position = 'dodge')+
                         theme_bw()+
                         facet_wrap(~Region, nrow = facet_nrow, scales = 'free')+
                         ggtitle(paste('Mouse', m))+
                         xlab('projection strength')

                       plot.out.m

                     })

  return(plot.out)


}



#' Posterior predictive checks based on multiple replicates
#'
#' @description
#' This function generates several replicates and compares the number of zero counts for each region, as well
#' as the distribution of non-zero counts for each region.
#'
#'
#' @param mcmc_run_all_output output from \code{HBMAP_mcmc} for the full algorithm.
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param N integer. Number of replicates.
#' @param regions.name a character vector of region names.
#'
#' @return a list of components:
#' \item{Y_prop}{Empirical projection strengths from the real data.}
#' \item{Y_rep_prop}{a list of empirical projection strengths for each replicate.}
#' \item{theta}{a list of parameters used to generate the replicates.}
#' \item{non.zero.plot}{a ggplot comparing the boxplots of non-zero counts between replicates and rea -data, for each region.}
#' \item{zero.plot}{a barplot comparing the number of zero counts between replicates and real data, for each region.}
#' @export
ppc_multiple <- function(mcmc_run_all_output,
                        Y,
                        N,
                        regions.name = NULL){



  length_per_chain <- length(mcmc_run_all_output$Z_output)

  C <- mcmc_run_all_output$C
  R <- mcmc_run_all_output$R
  M <- mcmc_run_all_output$M

  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }


  # Sum of counts for each neuron cell
  N_CM <- lapply(1:M,
                 function(m) colSums(Y[[m]]))

  # theta
  target.index <- sample(length_per_chain,
                         size = N,
                         replace = FALSE)


  theta <- lapply(1:N,
                  function(n){

                    list(q_j_star = mcmc_run_all_output$q_star_1_J_output[[target.index[n]]],
                         gamma_j_star = mcmc_run_all_output$gamma_star_1_J_output[[target.index[n]]],
                         Z = mcmc_run_all_output$Z_output[[target.index[n]]])
                  })

  # Simulate replicated data Y
  Y_rep <- lapply(1:N,
                  function(n){

                    Y_rep_n <- lapply(1:M,
                                      function(m){

                                        do.call(cbind, lapply(1:C[m],
                                                              function(c){

                                                                prob.c <- extraDistr::rdirichlet(n = 1,
                                                                                     theta[[n]]$q_j_star[theta[[n]]$Z[[m]][c],]*theta[[n]]$gamma_j_star[theta[[n]]$Z[[m]][c]])

                                                                rmultinom(n = 1, size = N_CM[[m]][c], prob = prob.c)
                                                              }))
                                      })

                    Y_rep_n
                  })



  #------------------------------------------ Some statistics --------------------------------------

  Y_prop <- lapply(1:M,
                   function(m){

                     matrix(as.vector(Y[[m]])/rep(colSums(Y[[m]]), each = R),
                            nrow = R)
                   })

  Y_rep_prop <- lapply(1:N,
                       function(n){

                         lapply(1:M,
                                function(m){

                                  matrix(as.vector(Y_rep[[n]][[m]])/rep(colSums(Y_rep[[n]][[m]]), each = R),
                                         nrow = R)
                                })

                       })

  Y_prop_df <- data.frame(Region = rep(regions.name, sum(C)),
                          data = 'observed',
                          value = as.vector(do.call(cbind,Y_prop)))

  Y_rep_prop_df <- do.call(rbind, lapply(1:N,
                                         function(n){

                                           data.frame(Region = rep(regions.name, sum(C)),
                                                      data = paste('replicated', n),
                                                      value = as.vector(do.call(cbind,Y_rep_prop[[n]])))

                                         }))

  # Convert into factor to keep brain region ordering
  Y_prop_df$Region <- factor(Y_prop_df$Region,
                             levels = regions.name)

  Y_rep_prop_df$Region <- factor(Y_rep_prop_df$Region,
                                 levels = regions.name)

  # Box plot
  non.zero.plot <- rbind(Y_prop_df,
                         Y_rep_prop_df) %>%
    dplyr::filter(value != 0) %>%
    ggplot(mapping = aes(x = Region,
                         y = value,
                         fill = data))+

    geom_boxplot(outlier.shape = NA)+
    theme_bw()+
    ylab('projection strength')

  zero.prop <- rbind(Y_prop_df,
                     Y_rep_prop_df) %>%
    dplyr::filter(value == 0) %>%
    dplyr::group_by(Region, data) %>%
    dplyr::summarise(N = dplyr::n())

  zero.plot <- ggplot(zero.prop, aes(x = Region,
                                     y = N,
                                     fill = data))+
    geom_bar(stat = 'identity',
             position = position_dodge())+
    theme_bw()+
    ylab('Number of zeros')


  # Return replicated Y
  return(list('Y_prop' = Y_prop,
              'Y_rep_prop' = Y_rep_prop,
              'theta' = theta,
              'non.zero.plot' = non.zero.plot,
              'zero.plot' = zero.plot))


}
