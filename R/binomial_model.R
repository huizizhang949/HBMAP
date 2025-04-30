#' Binomial Motif Analysis of MAPseq data
#'
#' @param data a list of matrices. Each is a region-by-neuron projection count matrix for an individual brain.
#' @param bin_thresh the binarization threshold. Barcode counts below threshold become 0. Barcode counts above threshold become 1.
#' @param pval_cutoff p-value threshold. Only motifs with p-values smaller than this threshold will be deemed "significant".
#' @param binom_model binomial model algorithm used. Options = c("one","two). Default = "two".
#'
#' @returns results (list of motif analysis results)
#' @export
#'
#' @examples
#' res = binomial_model(list_of_matrices, bin_thresh = 5, pval_cutoff = 0.05, binom_model = "two")
binomial_model <- function(data, bin_thresh = 0, pval_cutoff = 0.05, binom_model = "two"){

  M <- length(data)
  C <- sapply(1:length(data), function(m) ncol(data[[m]]))
  data <- do.call(cbind, data)

  R <- nrow(data)

  # Calculate Nf, and N_0
  Nf <- ncol(data)

  if (bin_thresh != 0){
    si <- sapply(1:nrow(data),
                 function(r) length(which(data[r,] >= bin_thresh)))
  }else{
    si <- sapply(1:nrow(data),
                 function(r) length(which(data[r,] > bin_thresh)))
  }

  N_0 <- round(pracma::fzero(function(x) Nf - x*(1-prod(1-si/x)), Nf)$x, 0)

  #Divergence Point - calculating p(e)

  if(binom_model == "one"){

    p_e <- pracma::fzero(function(x) Nf - N_0*(1-(1-x)^nrow(data)), 0.5)$x

    p_i <- rep(p_e, R)

    # Compute the expected count for different number of projecting motifs
    expected_number <- sapply(1:nrow(data),
                              function(r){
                                N_0*p_e^(r)*(1-p_e)^(R-r)
                              })

  }else{

    p_i <- si/N_0

  }

  # Find the set of unique clusters
  if (bin_thresh != 0){
    cluster_neuron <- sapply(1:ncol(data),
                             function(c) paste(rownames(data)[data[,c] >= bin_thresh],
                                               collapse = " "))
  }else{
    cluster_neuron <- sapply(1:ncol(data),
                             function(c) paste(rownames(data)[data[,c] > bin_thresh],
                                               collapse = " "))
  }

  # Unique set of cluster
  cluster_label <- unique(cluster_neuron)

  # Region names
  region_names <- rownames(data)

  # Neuron allocation
  allocation <- sapply(1:ncol(data),
                       function(c) which(cluster_label == cluster_neuron[c]))



  `%notin%` <- Negate(`%in%`)
  # Separate allocations
  C_cumsum <- c(0, cumsum(C))

  allocation <- lapply(1:M,
                       function(m) allocation[(C_cumsum[m]+1):(C_cumsum[m+1])])


  cluster_summary <- data.frame(cluster_name = cluster_label,
                                N_projecting_regions = sapply(1:length(cluster_label),
                                                              function(i) length(stringr::str_split(cluster_label[[i]], pattern = ' ')[[1]]))
  )

  if (binom_model == "two"){

    # Expected number of projections
    cluster_summary$size_expect <- sapply(1:nrow(cluster_summary),
                                          function(j){

                                            prod(c(p_i[which(rownames(data) %in%  stringr::str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])],
                                                   1-p_i[which(rownames(data) %notin%  stringr::str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])]))*N_0
                                          }
    )
  }else{

    # Expected number of projections
    cluster_summary$size_expect <- expected_number[cluster_summary$N_projecting_regions]

  }

  # Observed number of projections
  cluster_summary$size_observed <- sapply(1:nrow(cluster_summary),
                                          function(j) length(which(cluster_neuron == cluster_label[j])))

  cluster_summary <- cluster_summary %>%
    dplyr::mutate(cluster.type = ifelse(size_expect < size_observed, 'over-represented', 'under-represented'))

  #Calculate p-values

  if(binom_model == "one"){

    # p-value
    cluster_summary$p_value <- sapply(1:nrow(cluster_summary),
                                      function(j) pbinom(q = cluster_summary$size_observed[j],
                                                         size = ncol(data),
                                                         prob = (p_e)^cluster_summary$N_projecting_regions[j]*(1-p_e)^(R-cluster_summary$N_projecting_regions[j])))
  }else{

    cluster_summary$p_value <- sapply(1:nrow(cluster_summary),
                                      function(j) pbinom(q = cluster_summary$size_observed[j],
                                                         size = N_0,
                                                         prob = prod(c(p_i[which(rownames(data) %in%  stringr::str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])],
                                                                       1-p_i[which(rownames(data) %notin%  stringr::str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])]))
                                      ))

  }

  #update p-values
  updated_pvalues = c()

  for (x in cluster_summary$p_value){
    if (x > 0.5){
      pvalue_updated = 1 - x
      updated_pvalues = c(updated_pvalues, pvalue_updated)
    }else{
      updated_pvalues = c(updated_pvalues, x)
    }
  }

  cluster_summary$p_value = updated_pvalues

  if (binom_model == "one"){

    threshold_p_value <- -log(pval_cutoff/nrow(cluster_summary), base = 10)

  }else{

    threshold_p_value <- -log(pval_cutoff/(M*length(which(cluster_summary$N_projecting_regions > 1))), base = 10)

  }

  cluster_summary$significant <- ifelse(-log(cluster_summary$p_value, base = 10) < threshold_p_value,
                                        'insignificant',
                                        'significant')

  upset = upset_plot(cluster_summary)

  plot_output <- cluster_summary %>%
    dplyr::filter(N_projecting_regions > 1) %>%
    ggplot(mapping = aes(x = log(size_observed/size_expect, base = 2),
                         y = -log(p_value, base = 10)))+
    geom_point()+
    geom_hline(yintercept = threshold_p_value, linetype = 'dashed', color = 'grey')+
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey')+
    theme_bw()


  return(list('allocation' = allocation,
              'estimates' = p_i,
              'cluster_label' = cluster_label,
              'cluster_summary' = cluster_summary,
              'data' = data,
              'plot_output' = plot_output,
              'upset_plot' = upset))
}


#' Generate Motif Analysis Upset Plot
#'
#' @param cluster_summary data.frame describing key features of all identified motifs
#'
#' @returns Motif Analysis Upset Plot
#' @export
#'
#' @examples
#' chart = upset_plot(res$cluster_summary)
upset_plot = function(cluster_summary){

  cluster_summary$motif_num = seq_len(nrow(cluster_summary))

  sig = cluster_summary$significant
  type = cluster_summary$cluster.type

  gr = c()
  for (i in 1:length(sig)){

    sig_val = sig[i]
    type_val = type[i]

    if(sig_val == "significant" & type_val == "over-represented" ){
      gr = c(gr,1)
    }else if (sig_val == "significant" & type_val == "under-represented" ){
      gr = c(gr,2)
    }else if (sig_val == "insignificant"){
      gr = c(gr,3)
    }
  }

  `%notin%` <- Negate(`%in%`)
  #format colour vector
  if(1 %notin% gr){

    gr_cols = c("green", "grey", "brown", "gray")
    gr_cols_comb = c("green", "brown")

  }else if (2 %notin% gr){

    gr_cols = c("red","grey","brown","grey")
    gr_cols_comb = c("red", "brown")

  }else if (3 %notin% gr){

    gr_cols = c("red","grey","green","grey")
    gr_cols_comb = c("red", "green")

  }else{

    gr_cols = c("red", "grey", "green", "gray", "brown", "gray")
    gr_cols_comb = c("red", "green", "brown")

  }

  #reorder based on group
  cluster_summary$group = gr
  cluster_summary <- cluster_summary[order(cluster_summary$group, -cluster_summary$size_observed),]

  motifs = cluster_summary$cluster_name
  mnum = cluster_summary$motif_num
  obs = cluster_summary$size_observed
  exp = cluster_summary$size_expect

  ordered_group = cluster_summary$group

  #create Obs dataframe
  df1 = data.frame(
    motif = motifs,
    condition = "observed",
    count = obs,
    motif_num = mnum,
    group = ordered_group
  )
  df1$motif_num = factor(df1$motif_num, levels = df1$motif_num)

  #create Exp dataframe
  df2 = data.frame(
    motif = motifs,
    condition = "expected",
    count = exp,
    motif_num = mnum,
    group = ordered_group
  )
  df2$motif_num = factor(df2$motif_num, levels = df2$motif_num)

  df3 = rbind(df2, df1)

  #create paired barchart
  barchart = ggplot(df3, aes(fill = interaction(condition, group), y = count, x = motif_num)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.7) +
    ylab("Total Barcode Count") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0,0,0,0.5), 'lines')) +
    scale_fill_manual(values = gr_cols)

  #make combination matrix
  all_motifs = list()

  i = 1
  for (x in motifs){

    x_vector = list(strsplit(x, " "))[[1]]

    all_motifs[as.character(i)] = x_vector

    i = i + 1

  }

  #create combination matrix
  comb_df <- dplyr::bind_rows(lapply(names(all_motifs), function(set_name) {
    data.frame(Set = set_name, Element = all_motifs[[set_name]])
  }), .id = "id")

  comb_df <- dplyr::bind_rows(lapply(seq_along(all_motifs), function(i) {
    data.frame(Set = names(all_motifs)[i],
               Element = all_motifs[[i]],
               Group = ordered_group[i])
  }), .id = "id")

  comb_df$Set = as.integer(comb_df$Set)
  comb_df$Set = as.factor(comb_df$Set)

  # Create a column for grouping lines (each set)
  comb_df <- comb_df %>%
    dplyr::group_by(Set) %>%
    dplyr::mutate(row_number = dplyr::row_number())

  comb_mat = ggplot(comb_df, aes(x = Set, y = Element, group = Set)) +
    geom_point(size = 4, aes(colour=factor(Group))) +  # Dots
    geom_line(aes(group = Set, colour=factor(Group)), linewidth = 1) +  # Connecting lines
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = unit(c(0,0,0,0), 'lines'),
          legend.position = "none") +
    labs(x = "", y = "", title = "") + # no text
    scale_x_discrete(labels= mnum) + #correct x-axis labels
    scale_colour_manual(values = gr_cols_comb)

  figure <- ggpubr::ggarrange(barchart, comb_mat,
                      nrow = 2,
                      align = "v",
                      heights = c(1.5,1))
  return(figure)
}


#' Compare motif analysis results
#'
#' @param results_1 binomial model results object for Brain 1
#' @param results_2 binomial model results object for Brain 2
#'
#' @returns paired barchart comparing the estimated per-region projection probabilities of Brain 1 and Brain 2.
#' @export
#'
#' @examples
#' compare_brains(brain1_results, brain2_results)
compare_brains = function(results_1, results_2){

  estimates_1 = results_1$estimates
  estimates_2 = results_2$estimates

  regions = rownames(results_1$data)

  # Create a data frame in long format for ggplot
  df <- data.frame(
    Region = rep(regions, 2),
    Estimate = c(estimates_1, estimates_2),
    Brain = rep(c("B1", "B2"), each = length(estimates_1))
  )

  # Create the paired bar chart
  comp_figure = ggplot(df, aes(x = Region, y = Estimate, fill = Brain)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    labs(x = "Region", y = "Estimated Projection Probability") +
    theme_minimal() +
    scale_fill_manual(values = c("gray","black"))

  return(comp_figure)

}
