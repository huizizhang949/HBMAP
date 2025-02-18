#' Title
#'
#' @param proj_strength
#' @param Z
#'
#' @return
#' @export
#'
#' @examples
reorder_cluster <- function(proj_strength, Z){

  R <- ncol(proj_strength)
  M <- length(Z)
  C <- sapply(Z, length)
  # Relabel clusters based on the strength of projection
  strong.proj <- apply(proj_strength,
                       1,
                       function(x) which(x == max(x)))


  df0 <- lapply(1:R,
                function(r){

                  # Clusters with the strongest projection to region r
                  strong.proj.r <- which(strong.proj==r)


                  if(length(strong.proj.r) != 0){

                    qj_estimate_r <- matrix(proj_strength[strong.proj.r,],
                                            nrow = length(strong.proj.r))

                    output <- strong.proj.r[order(-qj_estimate_r[,r])]
                  }else{

                    output <- NULL

                  }
                }

  )

  # Reorder original label
  label.reorder <- unlist(df0)


  Z <- lapply(1:M,
              function(m){

                sapply(1:C[m], function(c) which(label.reorder == Z[[m]][c]))
              })


  return(list(label.reorder=label.reorder, Z.relabel=Z))

}

#' Title
#'
#' @param post_output
#' @param Z
#' @param regions.name
#' @param epsilon
#' @param cred.int
#'
#' @return
#' @export
#'
#' @examples
mcmc_reorder_cluster <- function(post_output,
                                 Z,
                                 regions.name = NULL,
                                 epsilon = 0.01,
                                 cred.int = 0.9){

  q_star_1_J_output <- post_output$q_star_1_J_output
  gamma_star_1_J_output <- post_output$gamma_star_1_J_output
  omega_J_M_output <- post_output$omega_J_M_output
  omega_output <- post_output$omega_output

  if(!post_output$post || (is.null(q_star_1_J_output) || is.null(gamma_star_1_J_output))){stop('Need to run the post-processing step for sampling q_star and gamma!')}

  # Dimensions
  J <- length(unique(unlist(Z)))
  R <- post_output$R
  C <- sapply(Z, length)
  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }

  ##-- Compute posterior mean
  proj_prob_mean <- Reduce("+", q_star_1_J_output)/length(q_star_1_J_output)
  gamma_mean <- Reduce("+", gamma_star_1_J_output)/length(gamma_star_1_J_output)

  ##-- Credible intervals
  proj_prob_lower <- matrix(0, nrow = J, ncol = R)
  proj_prob_upper <- matrix(0, nrow = J, ncol = R)
  proj_prob_med <- matrix(0, nrow = J, ncol = R)

  for(j in 1:J){
    for(r in 1:R){

      proj_prob_lower[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), (1-cred.int)/2)
      proj_prob_upper[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), 1-(1-cred.int)/2)
      proj_prob_med[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), 0.5)

    }
  }



  #----------------------- Total number of projection regions ---------------------------------

  neuron_projection_df <- NULL
  neuron_projection_df$q_star_1_J_output <- q_star_1_J_output
  neuron_projection_df$gamma_star_1_J_output <- gamma_star_1_J_output
  neuron_projection_df$Z <- Z
  neuron_projection_df$J <- J

  q_tilde_001 <- q_tilde(neuron_projection_df = neuron_projection_df,
                         epsilon = epsilon)


  reorder <- reorder_cluster(proj_strength = proj_prob_mean, Z = Z)

  original.label <- reorder$label.reorder

  Z <- reorder$Z.relabel

  # Change to the new labeling system
  q_star_1_J_output <- lapply(1:length(q_star_1_J_output),
                              function(t) q_star_1_J_output[[t]][original.label,])

  gamma_star_1_J_output <- lapply(1:length(gamma_star_1_J_output),
                                  function(t) gamma_star_1_J_output[[t]][original.label])

  proj_prob_mean <- proj_prob_mean[original.label,]
  proj_prob_lower <- proj_prob_lower[original.label,]
  proj_prob_upper <- proj_prob_upper[original.label,]
  proj_prob_med <- proj_prob_med[original.label,]

  gamma_mean <- gamma_mean[original.label]

  if((!is.null(omega_J_M_output)) && (!is.null(omega_output))){
    omega_J_M_output <- lapply(1:length(omega_J_M_output),
                               function(t) omega_J_M_output[[t]][original.label,])
    omega_output <- lapply(1:length(omega_output),
                           function(t) omega_output[[t]][original.label])
  }

  # Relabel q_tilde
  q_tilde_001 <- q_tilde_001[original.label,]

  colnames(q_tilde_001) <- regions.name
  rownames(q_tilde_001) <- paste('cluster', 1:J)

  cluster.labels <- sapply(1:J,function(j) {
    paste(j, ':', paste(colnames(q_tilde_001)[q_tilde_001[j,] >= 0.5], collapse = ','))
  })
  ##-- Output values
  return(list('J' = J,
              'R' = R,
              'regions.name' = regions.name,
              'q_star_1_J_output' = q_star_1_J_output,
              'proj_prob_mean' = proj_prob_mean,
              'proj_prob_med' = proj_prob_med,
              'proj_prob_lower' = proj_prob_lower,
              'proj_prob_upper' = proj_prob_upper,
              'Z' = Z,
              'gamma_mean' = gamma_mean,
              'gamma_star_1_J_output' = gamma_star_1_J_output,
              # 'cluster.label.summary' = cluster.label.summary,
              'q_tilde_001' = q_tilde_001,
              'cluster.labels' =  cluster.labels,
              'omega_J_M_output' = omega_J_M_output,
              'omega_output' = omega_output)
  )
}




#------ Computation of the probability that the projection strength of neurons -----
#------ in each cluster is greater than a small threshold for each brain region

q_tilde <- function(neuron_projection_df,
                    epsilon = 0.01){

  # Trace of q
  q_trace <- neuron_projection_df$q_star_1_J_output
  gamma_trace <- neuron_projection_df$gamma_star_1_J_output

  # Estimated allocation
  Z <- neuron_projection_df$Z

  # Trace length
  trace.length <- length(q_trace)

  df0 <- lapply(1:nrow(q_trace[[1]]),
                function(j){

                  df0_j <- lapply(1:trace.length,
                                  function(i){

                                    vec0 <- 1-pbeta(q = epsilon,
                                                    shape1 = q_trace[[i]][j,]*gamma_trace[[i]][j],
                                                    shape2 = (1-q_trace[[i]][j,])*gamma_trace[[i]][j],
                                                    lower.tail = TRUE)

                                    matrix(vec0,
                                           nrow = 1)
                                  })

                  df0_j <- do.call(rbind,
                                   df0_j)

                  matrix(round(colMeans(df0_j),3),
                         nrow = 1)

                })

  # Dimension: J x R
  df0 <- do.call(rbind,
                 df0)

  return(df0)
}
