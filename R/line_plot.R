#' Plot empirical projection strengths within each cluster
#'
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param Z a list of allocations (integers). Each is a vector of allocations for individual mouse.
#' @param cluster.labels a character vector of cluster labels.
#' @param regions.name a character vector of region names.
#' @param group.index optional. A list of integers denoting the grouping, same dimension as \code{Z}.
#' @param group.name a character for the title of the legend.
#' @param cluster.index optional. A vector of cluster labels (integers) to that will be plotted.
#' @param title the title of the plot.
#' @param facet_ncol number of columns in the figure. Default to 5.
#' @param col optional. A vector of color names for the grouping.
#'
#' @return A plot of \eqn{J} panels. Each panel shows the lineplots of empirical projection strengths within each cluster, with
#' cluster label shown on the top. If grouping is provided, mice are colored by groups.
#' @export
plot_empirical_ps <- function(Y, Z, cluster.labels = NULL, regions.name = NULL,
                              group.index = NULL, group.name = 'group',
                              cluster.index = NULL,
                              title = '', facet_ncol = 5, col = NULL){

  J <- length(unique(unlist(Z)))
  R <- nrow(Y[[1]])
  Y_cbind <- do.call(cbind, Y)
  M <- nrow(Y_cbind)
  C <- sapply(Z, length)
  if(is.null(group.index)){
    group.index <- rep(1,sum(unlist(C)))
  }

  df <- t(Y_cbind)
  # empirical projection strengths: C*M
  df <-  t(apply(df, 1, function(x){return(x/sum(x))}))

  if(is.null(cluster.labels)){
    cluster.labels <- 1:J
  }

  if(is.null(cluster.index)){
    cluster.index <- 1:J
  }

  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }

  bayesian_motif_pp <- lapply(cluster.index,
                              function(j){

                                data.frame(cluster = cluster.labels[j],
                                           pp = as.vector(t(df[which(unlist(Z) == j),])),
                                           region.name = regions.name,
                                           cell.index = rep(which(unlist(Z) == j), each = R),
                                           group.index = rep(group.index[which(unlist(Z) == j)], each = R))
                              })

  bayesian_motif_pp <- do.call(rbind, bayesian_motif_pp)
  gg <- bayesian_motif_pp %>%
    dplyr::mutate(cell.index = as.factor(cell.index),
           group.index = as.factor(group.index),
           cluster = factor(cluster,
                            levels = unique(bayesian_motif_pp$cluster)),
           region.name = factor(region.name, levels=regions.name)) %>%
    ggplot(mapping = aes(x = region.name,
                         y = pp,
                         color = group.index,
                         group = cell.index))+
    geom_line()+
    facet_wrap(~cluster, ncol = facet_ncol)+
    theme_bw()+
    xlab('region')+
    ylab('projection strength')+
    guides(color = guide_legend(title=group.name))+
    ggtitle(title)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  if(length(unique(group.index))==1){
    gg <- gg+theme(legend.position='none')
  }
  if(!is.null(col)){
    gg <- gg+scale_color_manual(values=col)
  }
  print(gg)
}

# ---------- Line plot of estimated projection strength within each cluster ----------
#' Summarize the number of regions each cluster projects to
#'
#' @param post_output_reorder output from \code{mcmc_reorder_cluster}.
#'
#' @return a list of the following components:
#' \item{plot_df}{a data frame with summary statistics for projection strengths, number of regions each cluster projects to and their names.}
#' \item{R, J}{number of regions, number of clusters.}
#' \item{regions.name}{a character vector of region names.}
#' @export
ps_summarize <- function(post_output_reorder){

  J <- post_output_reorder$J
  R <- post_output_reorder$R
  cluster.labels <- post_output_reorder$cluster.labels
  regions.name <- post_output_reorder$regions.name

  # Summary of projection probabilities
  loop.result <- lapply(1:J, function(j){

    data.frame(cluster = paste('cluster', j),
               cluster_label = cluster.labels[j],
               region = post_output_reorder$regions.name,
               projection.mean = post_output_reorder$proj_prob_mean[j,],
               projection.lower = post_output_reorder$proj_prob_lower[j,],
               projection.upper = post_output_reorder$proj_prob_upper[j,],
               projection.med = post_output_reorder$proj_prob_med[j,])
  })

  estimated.projection.df <- do.call(rbind, loop.result)

  # number of regions projected for each cluster
  number_of_projecting_regions <- lapply(1:J,
                                         function(j){

                                           return(data.frame(cluster = paste('cluster', j),
                                                             cluster_label = cluster.labels[j],
                                                             number_of_p = sum(post_output_reorder$q_tilde_001[j,] >= 0.5)))
                                         })
  number_of_projecting_regions <- do.call(rbind, number_of_projecting_regions)


  plot_df <- estimated.projection.df %>%
    dplyr::left_join(number_of_projecting_regions, by=c('cluster','cluster_label'))

  plot_df$number_of_p <- factor(plot_df$number_of_p,
                                levels = 1:R)

  return(list(plot_df = plot_df,
              R = R,
              J = J,
              regions.name = regions.name))
}


#' Plot the posterior mean of projection strengths with uncertainty given by credible intervals
#'
#' @param ps_summary output from \code{ps_summarize}.
#' @param cluster.index optional. A vector of cluster labels (integers) to that will be plotted.
#' @param title the title of the heatmap.
#' @param facet_ncol number of columns in the figure. Default to 5.
#' @param col optional. A vector of color names for the grouping.
#'
#' @return A plot of \eqn{J} panels. Each panel shows the posterior mean of projection strengths and credible intervals for each cluster, with
#' cluster label shown on the top. Clusters are colored by the number of projected regions.
#' @export
plot_estimated_ps <- function(ps_summary,
                              cluster.index = NULL,
                              title = '', facet_ncol = 5, col = NULL){

  plot_df <- ps_summary$plot_df
  J <- ps_summary$J
  R <- ps_summary$R
  regions.name <- ps_summary$regions.name

  if(is.null(cluster.index)){
    cluster.index <- 1:J
  }

  ind = sapply(plot_df$cluster, function(c){c %in% paste('cluster', cluster.index)})
  plot_df <- plot_df[ind,]

  gg <- ggplot(plot_df,mapping = aes(x = factor(region, levels = regions.name),
                               y = projection.med,
                               group = cluster,
                               color = number_of_p)) +
    geom_line()+
    geom_point()+
    geom_errorbar(aes(ymin = projection.lower,
                      ymax = projection.upper),
                  width = 0.1)+
    theme_bw()+
    ylim(c(0,1))+
    ylab('projection strength')+
    xlab('region')+
    labs(color = 'no. of regions', title = title) +
    facet_wrap(~factor(cluster_label, levels = unique(cluster_label)), ncol = facet_ncol)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=10))


    if(!is.null(col)){
      gg+scale_color_manual(values=col)
    }else{
      print(gg)
    }
}

