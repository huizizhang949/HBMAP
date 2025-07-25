#' Plot posterior similarity matrix
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @param psm.ind a list of posterior similarity matrices between pairwise mice. The first output from \code{similarity_matrix}.
#' @param psm.tot a single posterior similarity for all mice together. The second output from \code{similarity_matrix}.
#' @param group logical. If TRUE, the posterior similarity matrix will be separated by mice.
#' @param index Optional. An integer that indexes the mice. If provided, the posterior similarity matrix will only be plotted for the selected mouse.
#' @param method the method to perform hierarchical clustering to arrange posterior similarity matrix. see \code{hclust}.
#' @param title the title of the first plot.
#'
#' @return If \code{group} is TRUE and \code{index} is not provided, the function returns a heatmap of posterior similarity matrix, with mice separated by black solid lines,
#' the main diagonal blocks correspond to within-mouse similarity. If \code{index} is provided, the function returns a heatmap of posterior similarity matrix for a single mouse.
#' Otherwise a complete heatmap with no separation of mice is returned.
#' @export
#'
plotpsm <- function(psm.ind,
                    psm.tot,
                    group = TRUE,
                    index = NULL,
                    method="complete",
                    title = ''){

  if(any(psm.tot !=t(psm.tot)) | any(psm.tot >1) | any(psm.tot < 0) | sum(diag(psm.tot)) != nrow(psm.tot) ){
    stop("psm.tot must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}


  ##------------------------------ without distinguishing between datasets ------------------------------
  if((!group) & is.null(index)) {
    hc=hclust(as.dist(1-psm.tot), method = method, members = NULL)
    psm_hc=psm.tot
    n=nrow(psm.tot)
    psm_hc[1:n,]=psm_hc[hc$order,]
    psm_hc[,1:n]=psm_hc[,hc$order]


    # Plot 1
    plot.tot <- data.frame(neuron1 = rep(1:n, each = n),
                           neuron2 = rep(1:n, n),
                           posterior.similarity = as.vector(psm_hc)) %>%

      ggplot(mapping = aes(x = neuron1,
                           y = neuron2,
                           fill = posterior.similarity))+
      geom_tile()+
      theme_bw()+
      scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                           values = scales::rescale(c(0, 0.25, 1)),
                           limits = c(0,1))+
      xlab('neurons')+
      ylab('neurons')+
      labs(fill='posterior similarity')+

      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(expand = c(0,0))+
      ggtitle(title)

    return(plot.tot)
  }



  if(group & is.null(index)){
    # Dimension of datasets
    M <- length(psm.ind)
    n <- rep(0,M)
    for(m in 1:M){
      n[m] <- nrow(psm.ind[[m]][[m]])
    }
    n <- c(0,cumsum(n))

    ##----------------------------------- Distinguishing between datasets ----------------------------------
    psm_matrix_output <- matrix(0, nrow=max(n), ncol=max(n))
    hc_diag <- NULL

    ##-- For diagonal
    for(m in 1:M){
      hc_diag[[m]] <- hclust(as.dist(1-psm.ind[[m]][[m]]), method=method, members=NULL)
      psm_hc_diag <- psm.ind[[m]][[m]]
      psm_hc_diag[1:nrow(psm_hc_diag),] <- psm_hc_diag[rev(hc_diag[[m]]$order),]
      psm_hc_diag[,1:nrow(psm_hc_diag)] <- psm_hc_diag[,hc_diag[[m]]$order]
      psm_matrix_output[(n[m]+1):n[m+1], (max(n)-n[m+1]+1):(max(n)-n[m])] <- psm_hc_diag
    }

    ##-- For off-diagonal
    for(m1 in 2:M){
      for(m2 in 1:(m1-1)){
        psm_m1_m2 <- t(psm.ind[[m1]][[m2]])
        psm_m1_m2[1:nrow(psm_m1_m2),] <- psm_m1_m2[rev(hc_diag[[m2]]$order),]
        psm_m1_m2[,1:ncol(psm_m1_m2)] <- psm_m1_m2[,hc_diag[[m1]]$order]
        psm_matrix_output[(n[m2]+1):n[m2+1], (max(n)-n[m1+1]+1):(max(n)-n[m1])] <- psm_m1_m2
      }
    }

    for(m1 in 1:(M-1)){
      for(m2 in (m1+1):M){
        psm_m1_m2 <- psm.ind[[m2]][[m1]]
        psm_m1_m2[1:nrow(psm_m1_m2),] <- psm_m1_m2[rev(hc_diag[[m2]]$order),]
        psm_m1_m2[,1:ncol(psm_m1_m2)] <- psm_m1_m2[,hc_diag[[m1]]$order]
        psm_matrix_output[(n[m2]+1):n[m2+1], (max(n)-n[m1+1]+1):(max(n)-n[m1])] <- psm_m1_m2
      }
    }


    # Plot 2
    plot.ind <- data.frame(neuron1 = rev(rep(1:nrow(psm.tot), each = nrow(psm.tot))),
                           neuron2 = rep(1:nrow(psm.tot), nrow(psm.tot)),
                           posterior.similarity = as.vector(psm_matrix_output)) %>%

      ggplot(mapping = aes(x = neuron1,
                           y = neuron2,
                           fill = posterior.similarity))+
      geom_tile()+
      theme_bw()+
      scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                           values = scales::rescale(c(0, 0.25, 1)),
                           limits = c(0,1))+
      xlab('neurons')+
      ylab('neurons')+
      labs(fill='posterior similarity')+
      geom_vline(xintercept = n[-c(1,M+1)],
                 color = 'black',
                 linewidth = 1)+
      geom_hline(yintercept = n[-c(1,M+1)],
                 color = 'black',
                 linewidth = 1)+
      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(expand = c(0,0))+
      ggtitle(title)

    return(plot.ind)

  }

  if(!is.null(index)){
    m <- index
    hc=hclust(as.dist(1-psm.ind[[m]][[m]]), method = method, members = NULL)
    psm_hc=psm.ind[[m]][[m]]
    n=nrow(psm_hc)
    psm_hc[1:n,]=psm_hc[hc$order,]
    psm_hc[,1:n]=psm_hc[,hc$order]


    # Plot
    plot.sub <- data.frame(neuron1 = rep(1:n, each = n),
                           neuron2 = rep(1:n, n),
                           posterior.similarity = as.vector(psm_hc)) %>%

      ggplot(mapping = aes(x = neuron1,
                           y = neuron2,
                           fill = posterior.similarity))+
      geom_tile()+
      theme_bw()+
      scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                           values = scales::rescale(c(0, 0.25, 1)),
                           limits = c(0,1))+
      xlab('neurons')+
      ylab('neurons')+
      labs(fill='posterior similarity')+

      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(expand = c(0,0))+
      ggtitle(title)

    return(plot.sub)

  }

}

