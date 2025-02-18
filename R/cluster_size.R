#' Title
#'
#' @param clustering
#' @param group.index
#' @param group.name
#' @param cluster.index
#' @param title
#' @param col
#'
#' @return
#' @export
#'
#' @exampless
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
