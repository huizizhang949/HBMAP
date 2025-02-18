#' Title
#'
#' @param Y
#' @param post_output_reorder
#' @param allocation_prob
#' @param motif_indices
#' @param ncol
#'
#' @return
#' @export
#'
#' @examples
projection_probability <- function(Y,
                                   post_output_reorder,
                                   allocation_prob,
                                   motif_indices,
                                   ncol = 5){


  M <- length(Y)
  R <- dim(Y[[1]])[1]
  Z <- post_output_reorder$Z

  Y = do.call(cbind, Y)
  df0 <- t(Y)
  df0 = t(apply(df0, 1, function(x){return(x/sum(x))}))

  regions.name <- post_output_reorder$regions.name

  df <- lapply(motif_indices,
               function(j){

                 data.frame(cluster = paste(j, ':', paste(colnames(post_output_reorder$q_tilde_001)[post_output_reorder$q_tilde_001[j,] >= 0.5], collapse = ',')),
                            pp = as.vector(t(df0[which(unlist(post_output_reorder$Z) == j),])),
                            region.name = regions.name,
                            cell.index = rep(which(unlist(post_output_reorder$Z) == j), each = R),
                            # mouse.index = rep(mouse.index[which(unlist(post_output_reorder$Z) == j)], each = R),
                            allocation_prob = rep(unlist(allocation_prob)[which(unlist(post_output_reorder$Z) == j)], each = R))
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
    facet_wrap(~cluster, ncol = ncol)+
    #scale_color_gradient2(low = "yellow", high = "red") +
    theme_bw()+
    xlab('region')+
    ylab('projection strength') +
    labs(colour = "probability")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

# compute average allocation probabilities across all MCMC samples
allocation_probability = function(Y, post_output_reorder){

  M <- length(Y)
  C <- sapply(1:M, function(m) ncol(Y[[m]]))
  R <- dim(Y[[1]])[1]
  J <- post_output_reorder$J
  Z <- post_output_reorder$Z
  S <- length(post_output_reorder$omega_output)

  allocation_prob = list(M)
#' Title
#'
#' @param m
#'
#' @return
#' @export
#'
#' @examples
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

  return(list('allocation_probability' = allocation_prob,
              'allocation_probability_matrix' = allocation_prob_mat))

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
