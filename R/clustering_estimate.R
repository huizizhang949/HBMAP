#' Update MCMC samples of allocations such that all the cluster labels are consecutive
#'
#' @param mcmc_run_all_output output from \code{HBMAP_mcmc}.
#'
#' @return a matrix of MCMC samples of allocations, with rows for samples and columns for observations.
#' @export
#'
z_trace_updated <- function(mcmc_run_all_output){

  # Dimensions
  C <- mcmc_run_all_output$C
  M <- mcmc_run_all_output$M


  S = length(mcmc_run_all_output$Z_output)
  Zmat = matrix(unlist(mcmc_run_all_output$Z_output), S, sum(C), byrow = TRUE)
  k = apply(Zmat,1,function(x){length(unique(x))})
  z_trace_updated = matrix(0, S, sum(C))

  for(s in S:1){
    #reorder the configuration
    uniq_s=unique(Zmat[s,])
    for( h in 1:k[s]){
      z_trace_updated[s,Zmat[s,]==uniq_s[h]]=h
    }
  }

  return(z_trace_updated)
}


opt.clustering <- function(z_trace,
                           post_similarity,
                           max.k=NULL){

  # M = dim(z_trace)[1]
  # C = dim(z_trace)[2]



  opt.clust <- minVI(psm = post_similarity,
                     cls.draw = z_trace,
                     method = 'all',
                     max.k=max.k)

  opt.clust <- opt.clust$cl[1,]


  return(opt.clust)
}


salso_cluster_estimate <- function(z_trace,
                                   max.k=NULL){



  C = dim(z_trace)[2]

  ##-- Point estimate of z
  if(is.null(max.k)){max.k=C/10}
  z_point_estimate <- salso::salso(x = z_trace,
                                   maxNClusters = max.k)


  return(z_point_estimate)
}




#' Obtain an optimal clustering from MCMC samples
#'
#' @description
#' This function computes an optimal clustering by minimizing posterior expected variation of information.
#'
#'
#' @param z_trace a matrix of MCMC samples of allocations, with rows for samples and columns for observations. Can be output from \code{z_trace_updated}.
#' @param post_similarity a posterior similarity matrix for all mice. Can be output from \code{similarity_matrix}.
#' @param max.k optional. The maximum number of clusters that can be considered by the optimization algorithm
#'
#' @return a vector of estimated allocations (integers).
#' @export
opt.clustering.comb <- function(z_trace,
                                post_similarity,
                                max.k=NULL){



  #------------------------------------------------ Point estimate of z and trace of z from MCMC ------------------------------------------------

  # Result from minvi
  opt.clust.minvi <- opt.clustering(z_trace = z_trace,
                                    post_similarity = post_similarity,
                                    max.k=max.k)

  # Result from dlso
  opt.clust.dlso <- salso_cluster_estimate(z_trace = z_trace,
                                           max.k=max.k)


  #---------------------------------- Difference between each sample and the minvi estimated allocation -----------------------------------------

  vi.min <- apply(z_trace,1,
                   function(t){

                     vi.minvi.t <- clevr::variation_info(true = opt.clust.minvi,
                                                         pred = t, base=2)

                     vi.minvi.t
                   })


  #--------------------------------- Difference between each sample and dlso estimated allocation ------------------------------------------------

  vi.dlso <- apply(z_trace,1,
                    function(t){

                      # Similarity
                      vi.dlso.t <- clevr::variation_info(opt.clust.dlso,
                                                  t, base=2)

                      vi.dlso.t
                    })


  #-------------------------------- Return item ----------------------------------

  # In the case when minvi gives smaller VI than dlso
  if(mean(vi.min) < mean(vi.dlso)){

    # we export the minvi point estimate of z
    return(opt.clust.minvi)

  }else{

    # otherwise export the dlso point estimate of z
    return(opt.clust.dlso)
  }

}
