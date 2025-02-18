# Compute total variation distance in omega_JM between mice and plot a heatmap, with the option of add credible intervals
#' Title
#'
#' @param mcmc_run_all_output
#' @param hpd
#' @param prob
#' @param text_size
#'
#' @return
#' @export
#'
#' @examples
plot_tv_distance <- function(mcmc_run_all_output, hpd=FALSE, prob=0.95, text_size=2){

  M <- mcmc_run_all_output$M

  # function to compute the tv distance comparing all mice to mouse m, x is the matrix for omega_JM
  mytv_dist = function(x,m){
    xdim = dim(x)
    y = matrix(x[,m],xdim[1],xdim[2])
    return(0.5*colSums(abs(x-y)))
  }

  # Compute the posterior mean of tv distance, row m corresponds to comparisons between mouse m and other mice
  tv_mean = matrix(0,M,M)
  for (m in 1:M){
    tv_dist_m = lapply(mcmc_run_all_output$omega_J_M_output,mytv_dist, m )
    tv_dist_m = data.frame(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
    tv_mean[m,] = colMeans(tv_dist_m)
  }
  tv_mean = data.frame(tv_mean, row.names = as.factor(1:M))
  colnames(tv_mean) = as.factor(1:M)

  tv_mean = data.frame(tv_mean, "Mouse.1" = as.factor(1:M))
  colnames(tv_mean)[1:4] <- as.factor(1:M)
  tv_mean <-  tidyr::pivot_longer(tv_mean,
                           cols = !Mouse.1,
                           names_to = "Mouse.2",
                           values_to = "TV"
  )

  gg1 <- ggplot(tv_mean, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
    geom_tile() +
    labs(x = "Mouse",
         y = "Mouse") +
    theme_bw() +
    scale_fill_gradient2(high = "red", mid = "orange", low="white", midpoint = 0.5) +
    geom_text(aes(label = round(tv_mean$TV,3)), color = "black", size = text_size) +
    coord_fixed()


  if(hpd){
    # Compute HPD intervals
    tv_lower = matrix(0,M,M)
    tv_upper = matrix(0,M,M)
    for (m in 1:M){
      tv_dist_m = lapply(mcmc_run_all_output$omega_J_M_output,mytv_dist, m )
      tv_dist_m = coda::as.mcmc(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
      tv_hpd =  coda::HPDinterval((tv_dist_m),prob=prob)
      tv_lower[m,] = tv_hpd[,1]
      tv_upper[m,] = tv_hpd[,2]
    }

    tv_lower = data.frame(tv_lower, row.names = as.factor(1:M))
    colnames(tv_lower) = as.factor(1:M)

    tv_lower = data.frame(tv_lower, "Mouse.1" = as.factor(1:M))
    colnames(tv_lower)[1:M] = as.factor(1:M)
    tv_lower <-  tidyr::pivot_longer(tv_lower,
                              cols = !Mouse.1,
                              names_to = "Mouse.2",
                              values_to = "TV"
    )

    tv_upper = data.frame(tv_upper, row.names = as.factor(1:M))
    names(tv_upper) = as.factor(1:M)

    tv_upper = data.frame(tv_upper, "Mouse.1" = as.factor(1:M))
    names(tv_upper)[1:M] = as.factor(1:M)
    tv_upper <-  tidyr::pivot_longer(tv_upper,
                              cols = !Mouse.1,
                              names_to = "Mouse.2",
                              values_to = "TV"
    )

    labs = paste0(round(tv_mean$TV,3),'\n [',round(tv_lower$TV,2),',', round(tv_upper$TV,2),']')

    gg1 <- ggplot(tv_mean, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
      geom_tile() +
      labs(x = "Mouse",
           y = "Mouse") +
      theme_bw() +
      scale_fill_gradient2(high = "red", mid = "orange", low="white", midpoint = 0.5) +
      geom_text(aes(label = labs), color = "black", size = text_size) +
      coord_fixed()

  }else{
    tv_lower <- NULL; tv_upper <- NULL
  }

  print(gg1)
}
