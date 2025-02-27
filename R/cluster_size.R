#' Summarize the size of each cluster
#'
#' @description
#' This function shows a barplot of cluster size, with the option to color by groups/mice.
#'
#'
#' @param clustering a list of allocations (integers). Each is a vector of allocations for individual mouse.
#' @param group.index optional. A list of integers denoting the grouping, same dimension as \code{clustering}.
#' @param group.name a character for the title of the legend.
#' @param cluster.index optional. A vector of cluster labels (integers) to that will be plotted.
#' @param title the title of the barplot.
#' @param col optional. A vector of color names for the grouping.
#'
#' @return a barplot with cluster size given by the height of each bar. If grouping is provided, mice are colored by groups.
#' @export
#'
opt.clustering.frequency <- function(clustering,
                                     group.index = NULL,
                                     group.name = 'group',
                                     cluster.index = NULL,
                                     title = '',
                                     col = NULL){

  J <- length(unique(unlist(clustering)))
  M <- length(clustering)
  C <- sapply(clustering, length)
  C_cumsum <- c(0, cumsum(C))
  if(is.null(group.index)){
    group.index <- rep(1,sum(unlist(C)))
  }
  # turn group.index into a list for each dataset
  group.index.list <- lapply(1:M, function(m) group.index[(1+C_cumsum[m]):C_cumsum[m+1]])
  if(is.null(cluster.index)){
    cluster.index <- 1:J
  }

  cluster.names <- 1:J

  loop.result <- lapply(1:M, function(m){

    data.frame(cluster = factor(clustering[[m]],
                                levels = cluster.names),
               group = factor(group.index.list[[m]], levels = sort(unique(group.index)))
               )
  })

  # Frequency table
  z.frequency <- do.call(rbind, loop.result)

  # Bar chart
  gg <- z.frequency %>% dplyr::filter(as.numeric(cluster) %in% cluster.index) %>%
    ggplot(mapping = aes(x = cluster, fill = group))+
    geom_bar(stat="count")+
    theme_bw()+
    ylab('Number of neurons')+
    labs(fill=group.name)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(title)

  if(length(unique(group.index))==1){
    gg <- gg+theme(legend.position='none')
  }
  if(!is.null(col)){
    gg <- gg+scale_fill_manual(values=col)
  }
  print(gg)
}
