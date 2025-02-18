# ------- Find prominent motifs -------
#' Title
#'
#' @param post_output_reorder
#' @param thresh
#' @param prob
#'
#' @return
#' @export
#'
#' @examples
identify_prominent_motif <- function(post_output_reorder, thresh=0.02, prob = 0.95){

  omega_J_M_output <- post_output_reorder$omega_J_M_output
  omega_output <- post_output_reorder$omega_output
  C <- sapply(post_output_reorder$Z, length)

  if((is.null(omega_J_M_output) || is.null(omega_output))){stop('Need to run the post-processing step for sampling omega!')}

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
