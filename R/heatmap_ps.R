#' A heatmap of empirical projection strengths grouped by the inferred clustering
#'
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param Z a list of allocations (integers). Each is a vector of allocations for individual mouse.
#' @param regions.name a character vector of region names.
#' @param group.index optional. A list of integers denoting the grouping, same dimension as \code{Z}.
#' @param group.name a character for the title of the legend.
#' @param cluster.index optional. A vector of cluster labels (integers) to that will be plotted.
#' @param col optional. A vector of color names for the grouping.
#' @param title the title of the heatmap.
#'
#' @return a heatmap of empirical projection strengths, with regions on the x-axis and mice on the y-axis.
#' Clusters are separated by horizontal dashed lines. If grouping is provided, mice are colored by groups.
#' @export
heatmap_ps <- function(Y,
                       Z,
                       regions.name = NULL,
                       group.index = NULL,
                       group.name = 'group',
                       cluster.index = NULL,
                       col = NULL,
                       title = ''){


  J <- max(unlist(Z))
  R <- nrow(Y[[1]])
  M <- length(Y)
  C <- sapply(1:M, function(m) ncol(Y[[m]]))
  C_cumsum <- c(0, cumsum(C))

  if(is.null(regions.name)){
    regions.name <- paste('region', 1:R)
  }
  if(is.null(group.index)){
    group.index <- rep(1,sum(unlist(C)))
  }
  if(is.null(cluster.index)){
    cluster.index <- 1:J
  }

  # turn group.index into a list for each dataset
  group.index.list <- lapply(1:M, function(m) group.index[(1+C_cumsum[m]):C_cumsum[m+1]])

  Y.cbind <- do.call(cbind, Y)
  Y.prop.cbind <- matrix(as.vector(Y.cbind)/rep(colSums(Y.cbind), each = R),
                         nrow = R)

  Y.prop <- lapply(1:M, function(m) matrix(Y[[m]]/rep(colSums(Y[[m]]), each = R),
                                           nrow = R))


  # Re-ordered projection strength
  df0 <- lapply(cluster.index,
                function(j){


                  mx0 <- matrix(Y.prop.cbind[,which(unlist(Z)==j)],
                                nrow = R)

                  strong.proj.j <- which(rowSums(mx0) == max(rowSums(mx0)))[1]

                  df0_list <- lapply(1:M, function(m){

                    if(length(which(Z[[m]]==j)) != 0){

                      mx0_m <- matrix(Y.prop[[m]][,which(Z[[m]]==j)],
                                      nrow = R)

                      mx0_m[,order(mx0_m[strong.proj.j,])]
                    }

                  })

                  do.call(cbind, df0_list)

                })

  # group index
  group.index.df <- lapply(cluster.index,
                           function(j){
                             
                             mx0 <- matrix(Y.prop.cbind[,which(unlist(Z)==j)],
                                           nrow = R)
                             
                             strong.proj.j <- which(rowSums(mx0) == max(rowSums(mx0)))[1]

                             df0_list <- lapply(1:M, function(m){


                               if(length(which(Z[[m]]==j)) != 0){
                                 
                                 mx0_m <- matrix(Y1.prop[[m]][,which(Z[[m]]==j)],
                                                 nrow = R1)

                                 g0_m <- group.index.list[[m]][which(Z[[m]]==j)]
                                 # rep(m, length(which(Z[[m]]==j)))
                                 
                                 g0_m[order(mx0_m[strong.proj.j,])]
                               }

                             })

                             unlist(df0_list)

                           })

  group.index.df <- unlist(group.index.df)


  df0 <- do.call(cbind, df0)

  N.j <- sapply(cluster.index,
                function(j) length(which(unlist(Z)==j)))

  # Convert to data frame for plotting
  # subset for cluster.index
  index <- unlist(lapply(1:M,function(m) {
    sapply(Z[[m]],function(x) x %in% cluster.index)
  }))
  Y.prop.df <- data.frame(region = rep(regions.name,
                                       ncol(Y.prop.cbind[,index])),

                          neuron = rep(1:ncol(Y.prop.cbind[,index]),
                                       each = length(regions.name)),

                          projection.strength = as.vector(df0),

                          group = factor(rep(group.index.df, each = R),
                                         levels = 1:length(unique(group.index))))

  gg <- Y.prop.df %>%
    ggplot(mapping = aes(x = factor(region, levels = regions.name),
                         y = neuron,
                         fill = group))+

    geom_tile(mapping = aes(alpha = projection.strength))+
    scale_alpha_identity()+
    theme_bw()+
    xlab('region')+
    labs(fill=group.name)+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size=12),
          plot.background = element_blank() ,
          panel.grid.major = element_blank() ,
          panel.grid.minor = element_blank() ,
          panel.border = element_blank() ,
          panel.background = element_blank()) +
    #draws x and y axis line
    theme(axis.line = element_line(color = 'black'))+
    # rotate region names labels
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylab('neurons')+

    geom_hline(yintercept = cumsum(N.j)[-length(cumsum(N.j))],
               color = 'grey',
               linetype = 'dashed',
               linewidth = 0.5)+

    ggtitle(title)

  if(length(unique(group.index))==1){
    gg <- gg+theme(legend.position='none')
  }
  if(!is.null(col)){
    gg <- gg+scale_fill_manual(values=col)
  }
  print(gg)
}
