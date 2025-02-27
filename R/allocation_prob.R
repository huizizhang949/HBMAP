#' Plot the allocation probabilities for each motif
#'
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param Z a list of allocations (integers). Each is a vector of allocations for individual mouse.
#' @param cluster.labels a character vector of cluster labels.
#' @param regions.name a character vector of region names.
#' @param allocation_prob a list of allocation probabilities. Output from \code{allocation_probability}.
#' @param cluster.index optional. A vector of cluster labels (integers) to that will be plotted.
#' @param title the title of the plot.
#' @param facet_ncol number of columns in the figure. Default to 5.
#'
#' @return A plot of \eqn{J} panels. Each panel shows the lineplots of empirical projection strengths colored by allocation probabilities.
#'  Clusters labels are shown on the top.
#' @export
projection_probability <- function(Y, Z, cluster.labels = NULL, regions.name = NULL,
                                   allocation_prob,
                                   cluster.index = NULL,
                                   title = '', facet_ncol = 5){


  M <- length(Y)
  R <- dim(Y[[1]])[1]
  J <- length(unique(unlist(Z)))

  if(is.null(cluster.labels)){
    cluster.labels <- 1:J
  }

  if(is.null(cluster.index)){
    cluster.index <- 1:J
  }

  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }

  Y = do.call(cbind, Y)
  df0 <- t(Y)
  df0 = t(apply(df0, 1, function(x){return(x/sum(x))}))

  df <- lapply(cluster.index,
               function(j){

                 data.frame(cluster = cluster.labels[j],
                            pp = as.vector(t(df0[which(unlist(Z) == j),])),
                            region.name = regions.name,
                            cell.index = rep(which(unlist(Z) == j), each = R),
                            # mouse.index = rep(mouse.index[which(unlist(post_output_reorder$Z) == j)], each = R),
                            allocation_prob = rep(unlist(allocation_prob)[which(unlist(Z) == j)], each = R))
               })
  df <- do.call(rbind, df)


  df = df %>%
    dplyr::mutate(cell.index = as.factor(cell.index),
           # mouse.index = as.factor(mouse.index),
           cluster = factor(cluster,
                            levels = unique(df$cluster)),
           region.name = factor(region.name, levels = regions.name))

  # Line graph
  ggplot(df)+
    geom_line(mapping = aes(x = region.name,
                            y = pp,
                            colour = allocation_prob,
                            group = interaction(cell.index, allocation_prob)))+
    facet_wrap(~cluster, ncol = facet_ncol)+
    #scale_color_gradient2(low = "yellow", high = "red") +
    theme_bw()+
    xlab('region')+
    ylab('projection strength') +
    labs(colour = "probability")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}



#' Compute average allocation probabilities across all MCMC samples
#'
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param post_output_reorder output from \code{mcmc_reorder_cluster}.
#'
#' @return a list of numeric vectors. Each vector represents the probabilities of neurons in a mouse belonging to their optimal cluster.
#' @export
allocation_probability = function(Y, post_output_reorder){

  M <- length(Y)
  C <- sapply(1:M, function(m) ncol(Y[[m]]))
  R <- dim(Y[[1]])[1]
  J <- post_output_reorder$J
  Z <- post_output_reorder$Z
  S <- length(post_output_reorder$omega_output)

  allocation_prob = list(M)
  allocation_prob_mat = lapply(1:M, function(m){
    matrix(0,C[m], J)
  })


  for(s in 1:S){
    allocation_prob_mat_s = allocation_prob_given_q_gamma(post_output_reorder$omega_J_M_output[[s]],
                                                          post_output_reorder$q_star_1_J_output[[s]],
                                                          post_output_reorder$gamma_star_1_J_output[[s]],
                                                          Y)
    allocation_prob_mat = lapply(1:M, function(m){
      allocation_prob_mat_s[[m]]/S+allocation_prob_mat[[m]] })
  }


  for(m in 1:M){
    allocation_prob[[m]] = sapply(1:C[m], function(c){
      allocation_prob_mat[[m]][c,Z[[m]][c]]
    })
  }

  return(allocation_prob)

}

##--------------- Allocation probability given one mcmc draw ----------------
allocation_prob_given_q_gamma <- function(omega_J_M,
                                          q_star_1_J,
                                          gamma_1_J_star,
                                          Y){

  M <- length(Y)
  C <- unlist(lapply(Y, ncol))
  J <- nrow(omega_J_M)
  R <- ncol(q_star_1_J)

  # Create an empty list
  allocation.prob <- NULL

  for(m in 1:M){

    Prob_m <- compute_probability_cpp(Y = Y[[m]], omega_M = omega_J_M[,m], q_star_1_J = q_star_1_J, gamma_1_J_star = gamma_1_J_star)

    ##-- Store allocation probability
    allocation.prob[[m]] <- Prob_m

  }

  ##-- Return probability
  return(allocation.prob)
}
