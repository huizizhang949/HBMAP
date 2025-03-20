#' Identify prominent motifs
#'
#' @description
#' This function finds prominent motifs by computing the posterior probability of the global weights greater than a threshold.
#' Clusters with a high posterior probability are classified as prominent motifs.
#'
#'
#' @param post_output_reorder output from \code{mcmc_reorder_cluster}.
#' @param thresh a value between 0 and 1. The threshold for the global weights.
#' @param prob a value between 0 and 1. Clusters with posterior probability greater than \code{prob} are identified as prominent motifs.
#'
#' @return a list of: 1) a vector cluster labels (integers) identified as prominent motifs.
#' 2) a vector of posterior probabilities that the global weights are greater than \code{thresh}.
#' @export
identify_prominent_motif <- function(post_output_reorder, thresh=0.02, prob = 0.95){

  omega_output <- post_output_reorder$omega_output
  C <- sapply(post_output_reorder$Z, length)

  if(is.null(omega_output)){stop('Need to run the post-processing step for sampling omega!')}

  print(paste('With a threshold of', thresh, 'we identify motifs where we would expect at least', thresh*sum(C), 'neurons in that motif across all mice'))
  omega_mat = matrix(unlist(post_output_reorder$omega_output),
                     nrow = length(post_output_reorder$omega_output),
                     ncol = post_output_reorder$J,
                     byrow = TRUE)
  # Probability of the global weight greater than a threshold
  prob_greater_global = apply(omega_mat>thresh,2, mean)

  # identify prominent motifs with probability greater than a large threshold
  prominent_motifs = which(prob_greater_global>prob)

  return(list(index=prominent_motifs, prob=prob_greater_global))
}
