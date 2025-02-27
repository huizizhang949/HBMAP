#' K-means clustering for MAPseq data
#'
#' @description
#' This function implements k-means clustering for MAPseq data with different data transformations.
#'
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param k number of clusters to be found.
#' @param trans type of data transformation. Default to 'row-sum'. Can also be one of the following:
#' \itemize{
#' \item 'row-sum': normalize each neuron by its total counts.
#' \item 'log row-sum': after row-sum normalization, take log(y+1).
#' \item 'max-row': normalize each neuron by its maximum count.
#' \item 'cosine': normalize each neuron's counts to unit length.
#' }
#' @param restart number of restart in k-means algorithm.
#' @param iter.max the maximum number of iterations allowed in k-means.
#'
#' @return a list of allocations (integers). Each is a vector of allocations for individual mouse.
#' @export
#'
#' @examples
#' set.seed(43)
#' Z <- k_means_exon(Y = Y, k = 20, trans = 'log row-sum', restart = 20, iter.max = 100)
k_means_exon <- function(Y, k, trans = 'row-sum',
                         restart = 50, iter.max = 100){

  # Number of mouse
  M <- length(Y)

  # Number of neurons in each cluster
  C <- sapply(1:M, function(m) ncol(Y[[m]]))

  # Cumulative number of neurons in a vector form
  C_cumsum <- c(0, cumsum(C))

  Y_cbind <- do.call(cbind, Y)

  if(trans == 'row_sum'){
    # Projection probabilities, neurons in the rows
    df <- t(apply(t(Y_cbind), 1, function(x){return(x/sum(x))}))
  }else if (trans == 'log row-sum'){

    df <- t(apply(t(Y_cbind), 1, function(x){return(x/sum(x))}))
    df <- log(df+1)
  }else if (trans == 'max-row'){

    df <- t(apply(t(Y_cbind), 1, function(x){return(x/max(x))}))
  }else if (trans == 'cosine'){

    # df <- t(apply(t(Y_cbind), 1, function(x){return(x/max(x))}))
    df <- t(Y_cbind)
    df <- cosine_normalize(mat = df)


  }


  # Neuron allocation
  allocation <- kmeans(df, centers = k, nstart = restart, iter.max = iter.max)

  # Neuron allocations
  Z <- lapply(1:M,
              function(m) allocation$cluster[(C_cumsum[m]+1):C_cumsum[m+1]])

  return(Z)
}

cosine_normalize <- function(mat) {

  # Normalize each row by its magnitude
  return(mat / sqrt(rowSums(mat^2)))
}


#' Compute average projection strengths within each cluster
#'
#' @param Y a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.
#' @param Z a list of allocations (integers). Each is a vector of allocations for individual mouse.
#'
#' @return a matrix of projection strengths, with rows corresponding to clusters and columns for regions.
#' @export
#'
#' @examples
#' average_ps <- avg_ps(Y = Y, Z = Z)
avg_ps <- function(Y, Z){

  Z_unlist <- unlist(Z)
  Y_cbind <- do.call(cbind, Y)

  M <- length(Y)

  df <- t(Y_cbind)
  # empirical projection strengths: C*M
  df <-  t(apply(df, 1, function(x){return(x/sum(x))}))

  # Calculate the average projection strength of neurons in each cluster
  average_q <- lapply(sort(unique(Z_unlist)),
                      function(j){

                        data_j <- matrix(df[which(Z_unlist == j),], nrow = sum(Z_unlist==j))

                        matrix(colMeans(data_j), nrow = 1)
                      })

  # J * M
  average_q <- do.call(rbind, average_q)

  return(average_q)
}










