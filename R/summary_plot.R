#' A summary plot for motifs
#'
#' @import ggplot2
#'
#' @param params a list of arguments:
#' \itemize{
#' \item \code{proj_prob}, posterior mean of projection strengths (a \eqn{J \times R} matrix).
#' \item \code{omega_JM}, posterior mean of local weights (a \eqn{J \times M} matrix).
#' \item \code{omega}, posterior mean of global weights (a vector of length \eqn{J}).
#' }
#' @param eps a list of posterior expected projection strength for all mice. Each component is the output from \code{post_epstrength}.
#' @param prominent_motifs_prob probabilities that the global weights are greater than a threshold. output from \code{identify_prominent_motif}.
#' @param prominent_motifs_thresh a value between 0 and 1. Default to 0.95.
#' Clusters with \code{prominent_motifs_prob} greater than \code{prominent_motifs_thresh} are identified as prominent motifs. see \code{identify_prominent_motif}.
#' @param global_weight_thresh a value between 0 and 1 to threshold the expected global weight.
#' @param data.source.label a character for the title of the legend.
#' @param regions.name a character vector of region names.
#' @param col_bar A vector of color names to represent difference mice.
#' @param col_mat A vector of 2 color names to represent low and high values of projection strengths.
#' @param height height of the barplot on the top, for local and global weights. see \code{anno_barplot}.
#' @param width width of the barplot on the left, for posterior expected projection strength. see \code{anno_barplot}.
#' @param barwidth relative width of the bars, should be smaller than one. see \code{anno_barplot}.
#' @param legend boolean. If \code{TRUE}, a legend will be shown.
#' @param legend_x,legend_y. the \code{x} and \code{y} position of the legends, measured in current viewport. see \code{draw} from \code{ComplexHeatmap} package.
#'
#' @return A summary plot. The middle shows the heatmap of projection strengths for each motif (column), with global and mouse-specific weights attached
#' to each motif visualized above and the mouse-specific expected projection strength to each region on the left.
#' The first block identifies prominent motifs. The second block represents motifs having a lower probability of the global weight being greater
#' than a threshold (\code{thresh} in \code{identify_prominent_motif}) but the expected global weight exceeds \code{global_weight_thresh}.
#' @export
plot_summary <- function(params, eps,
                         prominent_motifs_prob,
                         prominent_motifs_thresh = 0.95,
                         global_weight_thresh = 0.01,
                         data.source.label = 'Mouse', regions.name = NULL,
                         col_bar=NULL, col_mat=c('white','blue'),
                         height = 6, width = 2, barwidth = 0.99,
                         legend = TRUE, legend_x = 0.9, legend_y = 0.95){

  q <- params$proj_prob; omega_JM <- params$omega_JM; omega <- params$omega
  J <- nrow(q); R <- ncol(q); M <- length(eps)
  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }
  if(is.null(col_bar)){
    col_bar <- RColorBrewer::brewer.pal(M, 'Spectral')
  }

  # Matrix of expected projection strength across clusters (regions x clusters)
  df = t(as.matrix(q))
  rownames(df) <- regions.name
  colnames(df) <- 1:J
  # Matrix of expected weight across clusters for the global and local weights for each mouse  (clusters x M+1)
  weightmat = matrix(cbind(omega, omega_JM), nrow = J, ncol = M+1)
  # Matrix of expected projection strength for each mouse (regions x M)
  psmat = as.matrix(do.call(cbind,eps))
  # vector of probability that the global weight is greater than a threshold (length = num of clusters)
  pg = prominent_motifs_prob

  # Colorscale for heatmap (q)
  col_fun = colorRamp2::colorRamp2(c(0, 1), col_mat)


  # Top annotation: expected global and local weights for each mouse
  # Note: if more mice, add another color. Red indicates global weights
  colnames(weightmat) = c("Global", paste(data.source.label, c(1:M)))
  column_ha  = ComplexHeatmap::HeatmapAnnotation(w = ComplexHeatmap::anno_barplot(weightmat,
                                                  beside = TRUE, attach = TRUE,
                                                  barwidth = barwidth,
                                                  gp = grid::gpar(fill = c("red", col_bar)),
                                                  height = grid::unit(height, "cm"))
  )
  # If you only want the global weights in the top annotation
  #column_ha = HeatmapAnnotation(global = anno_barplot(mcmc_omega$mean_omega_J))

  # Left annotation: expected projection strength for each mouse
  row_ha  = ComplexHeatmap::rowAnnotation(p = ComplexHeatmap::anno_barplot(psmat,
                                           beside = TRUE, attach = TRUE,
                                           barwidth = barwidth,
                                           gp = grid::gpar(fill = col_bar),
                                           border = FALSE,
                                           width = grid::unit(width, "cm"))
  )

  # Heat map of projection strength for each cluster with annotation
  print(ComplexHeatmap::Heatmap(df,
          cluster_rows = FALSE,
          row_names_side = "left",
          cluster_columns = FALSE,
          show_column_names = TRUE,
          heatmap_legend_param = list(title = ""),
          col = col_fun,
          column_order = order(weightmat[,1], decreasing = TRUE),
          column_split = as.factor(2*(pg<prominent_motifs_thresh&weightmat[,1]<global_weight_thresh) +
                                     (pg<prominent_motifs_thresh&weightmat[,1]>global_weight_thresh)),
          column_title_gp = grid::gpar(fontsize=0),
          border = TRUE,
          top_annotation = column_ha,
          left_annotation = row_ha
  ))
  if(legend){
    lgd = ComplexHeatmap::Legend(labels = colnames(weightmat), title = "",
                 legend_gp = grid::gpar(fill = c("red",col_bar)),
                 border = "black")
    ComplexHeatmap::draw(lgd,x = grid::unit(legend_x, "npc"), y = grid::unit(legend_y, "npc"),just = c("right", "top"))
  }
}




#' Compute posterior expected projection strength for a new neuron
#'
#' @param m integer. Index of the mouse.
#' @param mcmc_run_all_output output from \code{HBMAP_mcmc} for the full algorithm.
#'
#' @return a numeric vector of posterior expected projection strength to each region, if the new neuron is from mouse \code{m}.
#' @export
post_epstrength = function(m,mcmc_run_all_output){

  R= mcmc_run_all_output$R
  I = length(mcmc_run_all_output$alpha_output)

  #output
  eps_mat = matrix(0,I,R)

  for(i in 1:I){
    wjm = mcmc_run_all_output$omega_J_M_output[[i]][,m]
    qstarjm = mcmc_run_all_output$q_star_1_J_output[[i]]
    gammastarjm = mcmc_run_all_output$gamma_star_1_J_output[[i]]
    eps_mat[i,] = epstrength(qstarjm*gammastarjm, wjm)
  }
  eps = apply(eps_mat,2,mean)
  return(eps)
}

# ----- Compute expected number of counts/N for each region for a DirMult mixture -----
epstrength <- function(alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  eps = matrix(0,R,1)
  for (j in 1:J){
    eps = eps+w[j]*alpha[j,]/sum(alpha[j,])
  }
  return(eps)
}

