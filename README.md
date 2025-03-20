
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HBMAP

<!-- badges: start -->
<!-- badges: end -->

The goal of HBMAP is to implement a model-based clustering approach
through hierarchical Bayesian mixtures for high-throughput mapping of
axon projections. The model allows for identifying groups of neurons
with unique projection patterns across different brain regions, and
investigates differences across experiments.

## Installation

You can install the development version of HBMAP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("huizizhang949/HBMAP")
```

First, we load the data and necessary functions for MCMC sampling.

``` r
# load required packages
library(ggplot2)
library(ggrepel)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

``` r
library(HBMAP)
```

``` r
data("data_Hans")

M <- length(data_Hans)
C <- sapply(1:M, function(m) ncol(data_Hans[[m]]))
R <- dim(data_Hans[[1]])[1]

# Mouse index
mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]),
                 rep(4, C[4]))
```

Set empirical parameters and initial clustering.

``` r
# ---------------- Empirical parameters --------
#Initialize the clustering and set hyperparameters empirically from the data
data_Hans_cbind <- do.call(cbind, data_Hans)

C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_Hans[[m]]))))

# initialize with k-means with cosine distance
df <- lapply(1:M,function(m){
  normalized_mat <- apply(data_Hans[[m]], 2, function(col) col / max(col))})

df <- t(do.call(cbind, df))

cosine_normalize <- function(mat) {
  # Normalize each row by its magnitude
  return(mat / sqrt(rowSums(mat^2)))
}

df <- cosine_normalize(mat = df)

wss = apply(matrix(seq(6,40,1),ncol=1),1,function(x){
  kmeans_result <- kmeans(df, centers = x, nstart = 25)
  kmeans_result$tot.withinss
})
ggplot() + geom_point(aes(x=seq(6,40,1),y=wss)) + theme_bw()
```

<img src="man/figures/README-empirical estiamtes-1.png" width="60%" style="display: block; margin: auto;" />

``` r

# Based on the plot, 20 seems like a reasonable number of clusters to start with
# initial clustering
kmeans_result <- kmeans(df, centers = 20, iter.max = 100, nstart = 25)$cluster
# Set the truncation to be larger to allow the model to explore more clusters
J = 30

clustinit <- lapply(1:M,
                    function(m) kmeans_result[(C_cumsum[m]+1):C_cumsum[m+1]])

# Empirical choice of alpha parameters 
a_alpha0  <- mean(unlist(lapply(clustinit,function(x){length(unique(x))}))/log(C))
a_alpha  <- length(unique(unlist(clustinit)))/log(sum(unlist(lapply(clustinit,function(x){length(unique(x))}))))

# Set the prior expectation of gamma to scale with the total counts
a_gamma = median(unlist(lapply(data_Hans,colSums)))
b_gamma = 2
```

# Infer the clustering

``` r
# ---- parameters to pass to the main function ------
# mcmc setup
mcmc_list = list(number_iter = 5000, thinning = 1, burn_in = 4000, adaptive_prop = 0.0001,
                 auto_save = FALSE,
                 save_path = NULL,
                 save_frequency = 1000
                 )
# prior parameters, default values will be used if not provided
prior_list = list(a_gamma = a_gamma, b_gamma = b_gamma, lb_gamma = 1, a = 2, tau = 0.4, nu = 1/20,
                  a_alpha = a_alpha, b_alpha = 1, a_alpha0 = a_alpha0, b_alpha0 = 1)


# ------- Run the full model ---------
set.seed(3)
mcmc_all_hans <- HBMAP_mcmc(Y = data_Hans, J = J, mcmc = mcmc_list, prior = prior_list, 
                            Z.init = clustinit, verbose = FALSE)
```

## Traceplots and acceptance probability

``` r
# ------- MCMC check -----------
ind <- seq(1,mcmc_all_hans$output_index,by=1)

Zmat = matrix(unlist(mcmc_all_hans$Z_output), length(mcmc_all_hans$Z_output), sum(C),byrow = TRUE)

## ----- Number of occupied components -------
par(mfrow=c(1,1))
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')
```

<img src="man/figures/README-traceplots1-1.png" width="60%" style="display: block; margin: auto;" />

``` r
## ---- Concentration parameters ------
par(mfrow=c(1,2))
plot(mcmc_all_hans$alpha_zero_output, type = 'l',ylab='',main='alpha0')
plot(mcmc_all_hans$alpha_output, type = 'l',ylab='',main='alpha')
```

<img src="man/figures/README-traceplots2-1.png" width="80%" style="display: block; margin: auto;" />

``` r
## ----- q, gamma and weights may have label switching, better check their convergence from the post-processing step below -------
# j=5
# r=4
# plot(unlist(lapply(mcmc_all_hans$q_star_1_J_output[ind], function(q){q[j,r]})), type = 'l', ylab='', main=paste0('q[',j,',',r,']'))
# j=5
# plot(unlist(lapply(mcmc_all_hans$gamma_star_1_J_output[ind], function(q){q[j]})), type = 'l', ylab='',main=paste0('gamma[',j,']'))
# j=29
# m=4
# plot(unlist(lapply(mcmc_all_hans$omega_J_M_output[ind], function(q){q[j,m]})), type = 'l', ylab='', main=paste0('w[',j,',',m,']'))
# j=1
# plot(unlist(lapply(mcmc_all_hans$omega_output[ind], function(q){q[j]})), type = 'l', ylab='', main=paste0('w[',j,']'))


## ------ Acceptance rate ------
par(mfrow=c(2,3))
plot(unlist(mcmc_all_hans$acceptance_prob$q_star),type = 'l', ylab='', main='acceptance_q')
plot(unlist(mcmc_all_hans$acceptance_prob$gamma_star), type = 'l', ylab='', main='acceptance_gamma')
plot(unlist(mcmc_all_hans$acceptance_prob$alpha), type = 'l', ylab='', main='acceptance_alpha')
plot(unlist(mcmc_all_hans$acceptance_prob$alpha_zero), type = 'l', ylab='', main='acceptance_alpha0')
plot(unlist(mcmc_all_hans$acceptance_prob$omega), type = 'l', ylab='', main='acceptance_w')
```

<img src="man/figures/README-traceplots2-2.png" width="80%" style="display: block; margin: auto;" />

## Optimal clustering

``` r

# Posterior similarity matrix
psm_hans = similarity_matrix(mcmc_run_all_output = mcmc_all_hans)


# Reordered posterior samples of z
hans_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_hans)


# optimal clustering
set.seed(1)
hans_Z <- opt.clustering.comb(z_trace = hans_z_reordered,
                              post_similarity = psm_hans$psm.combined,
                              max.k = max(k))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))

hans_Z <- lapply(1:M,
                function(m) hans_Z[(C_cumsum[m]+1):C_cumsum[m+1]])
```

Posterior similarity matrix can be used to quantify the uncertainty in
clustering.

``` r
# Plot of posterior similarity matrix
psm_hans_plot <- plotpsm(psm.ind = psm_hans$psm.within,
                         psm.tot = psm_hans$psm.combined)

psm_hans_plot$plot.ind
```

<img src="man/figures/README-psm-1.png" width="60%" style="display: block; margin: auto;" />

# Post-processing step with a fixed clustering

``` r

mcmc_list = list(number_iter = 5000, thinning = 1, burn_in = 4000, adaptive_prop = 0.0001,
                 auto_save = FALSE,
                 save_path = NULL,
                 save_frequency = 1000
                 )
prior_list = list(a_gamma = a_gamma, b_gamma = b_gamma, lb_gamma = 1, a = 2, tau = 0.4, nu = 1/20,
                  a_alpha = a_alpha, b_alpha = 1, a_alpha0 = a_alpha0, b_alpha0 = 1)

# for the post-processing step, there is no label switching and we can make inference for each cluster based on q, gamma and weights
post_list = list(run_omega = TRUE, run_q_gamma = TRUE)

set.seed(3)
mcmc_hans_post <- HBMAP_mcmc(Y = data_Hans, mcmc = mcmc_list, prior = prior_list, 
                             Z.fix = hans_Z, post = post_list, verbose = FALSE)
```

We reorder the clusters based on estimated projection strengths.

``` r
## ----- Reorder clusters based on estimated projection strengths ---------
### ----- return reordered samples as well for q, gamma and omegas ------
### ----- return summary statistics for q and gamma, and cluster labels (projected regions) ------
### ----- this is based on the probability of projections strength greater than a threshold (q_tilde) -------------
mcmc_hans_post_reorder <- mcmc_reorder_cluster(post_output = mcmc_hans_post, 
                                               Z = hans_Z, regions.name = rownames(data_Hans[[1]]))

# Neuron allocations after reordering
hans_Z_reordered <- mcmc_hans_post_reorder$Z
```

Below we summarize cluster size, visualize neuron projection strengths
within each cluster and estimated projection strength.

``` r
# you can choose your preferred color palette (change 'col' argument)
## ------ cluster size: number of neurons in each cluster, colored by group (mouse/injection site) ---------
opt.clustering.frequency(clustering = hans_Z_reordered, group.index = mouse.index, group.name = 'mouse', 
                         title = 'Cluster size')
```

<img src="man/figures/README-cluster summary and line plots-1.png" width="70%" style="display: block; margin: auto;" />

``` r

## ---------- Heatmap of empirical projection strength of neurons in each cluster, colored by group (mouse/injection site) ----------
heatmap_ps(Y = data_Hans, Z = hans_Z_reordered, regions.name = rownames(data_Hans[[1]]), 
           group.index = mouse.index, group.name = 'mouse',
           cluster.index = 1:length(unique(unlist(hans_Z_reordered))), title = '')
```

<img src="man/figures/README-cluster summary and line plots-2.png" width="70%" style="display: block; margin: auto;" />

``` r


## --------- Line plot for of empirical projection strengths within each cluster ------------------------
plot_empirical_ps(Y = data_Hans, Z = hans_Z_reordered, 
                  cluster.labels = mcmc_hans_post_reorder$cluster.labels,
                  regions.name = rownames(data_Hans[[1]]),
                  group.index = mouse.index, group.name = 'mouse',
                  cluster.index = 1:length(unique(unlist(hans_Z_reordered))),
                  title = 'Empirical projection strength', facet_ncol = 5)
```

<img src="man/figures/README-cluster summary and line plots-3.png" width="70%" style="display: block; margin: auto;" />

``` r

# ---------- Line plot of estimated projection strength within each cluster ----------
# get summary statistics
ps_summary <- ps_summarize(post_output_reorder = mcmc_hans_post_reorder)
# plot
plot_estimated_ps(ps_summary = ps_summary, cluster.index = 1:length(unique(unlist(hans_Z_reordered))), 
                  title = 'Estimated projection strength', facet_ncol = 5)
```

<img src="man/figures/README-cluster summary and line plots-4.png" width="70%" style="display: block; margin: auto;" />

## Prominent motifs

Below we find prominent motifs based on the posterior samples of the
global weights.

``` r
# ------- Find prominent motifs: Probability of the global weight greater than a threshold ------- 
# prominent motifs are defined as those with a high probability (> 0.95)
# below the function returns the indices of prominent motifs and the posterior probability of the global weight greater than thresh
prominent_motifs <- identify_prominent_motif(post_output_reorder = mcmc_hans_post_reorder, 
                                             thresh = 0.02, prob = 0.95)
#> [1] "With a threshold of 0.02 we identify motifs where we would expect at least 11.06 neurons in that motif across all mice"
```

We could plot the estimated and empirical projection strengths for these
prominent motifs only.

``` r
# redo the lines plots of empirical and estimated projection strengths for prominent motifs only
plot_empirical_ps(Y = data_Hans, Z = hans_Z_reordered, 
                  cluster.labels = mcmc_hans_post_reorder$cluster.labels,
                  regions.name = rownames(data_Hans[[1]]),
                  group.index = mouse.index, group.name = 'mouse', 
                  cluster.index = prominent_motifs$index, 
                  title = 'Empirical projection strength (prominent motif)', 
                  facet_ncol = 4)
```

<img src="man/figures/README-plot prominent motifs-1.png" width="70%" style="display: block; margin: auto;" />

``` r
plot_estimated_ps(ps_summary = ps_summary, cluster.index = prominent_motifs$index, 
                  title = 'Estimated projection strength (prominent motif)', facet_ncol = 5)
```

<img src="man/figures/README-plot prominent motifs-2.png" width="70%" style="display: block; margin: auto;" />

In addition, we can investigate the posterior probabilities of
allocation to these prominent motifs as well as the uncertainty from the
posterior similarity matrix.

``` r
# line plots of empirical projection strengths colored by allocation probabilities
allo.prob = allocation_probability(post_output_reorder = mcmc_hans_post_reorder, 
                                   Y = data_Hans)

projection_probability(Y = data_Hans, Z = hans_Z_reordered, 
                       cluster.labels = mcmc_hans_post_reorder$cluster.labels,
                       regions.name = rownames(data_Hans[[1]]),
                       allocation_prob = allo.prob,
                       cluster.index = 1:length(unique(unlist(hans_Z_reordered))))
```

<img src="man/figures/README-allocation prob-1.png" width="60%" height="60%" style="display: block; margin: auto;" />

``` r
# plot the posterior similarity matrix for prominent motifs only
# to install the superheat package:
# install.packages("devtools")
# devtools::install_github("rlbarter/superheat")
library(superheat)
Z_factored = as.character(unlist(hans_Z_reordered))
ind = sapply(unlist(hans_Z_reordered), function(z){
  any(z == prominent_motifs$index)
})
superheat(psm_hans$psm.combined[ind,ind],
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          heat.pal = c("white", "yellow", "red"),
          heat.pal.values = c(0,.5,1),
          membership.rows = Z_factored[ind],
          membership.cols = Z_factored[ind],
          bottom.label.text.size = 4,
          left.label.text.size = 4)
```

<img src="man/figures/README-superheat heatmap-1.png" width="60%" height="60%" style="display: block; margin: auto;" />

## Variable motifs

We identify variable motifs by comparing the variance of the
brain-specific weights across mice with a null model.

``` r
mcmc_list = list(number_iter = 5000, thinning = 1, burn_in = 4000, adaptive_prop = 0.0001,
                 auto_save = FALSE,
                 save_path = NULL,
                 save_frequency = 1000
)
post_list <- list(run_omega=TRUE, run_q_gamma=FALSE)
# increase N to obtain more reliable results, e.g. N=200
set.seed(2)
local_weights_analysis_vc <- local_weights_analysis(N = 50, Z = mcmc_hans_post_reorder$Z, 
                                                    omega_output = mcmc_hans_post_reorder$omega_output,
                                                    omega_J_M_output = mcmc_hans_post_reorder$omega_J_M_output, 
                                                    prior = prior_list, mcmc = mcmc_list, 
                                                    post = post_list, n_cores=4, verbose = FALSE)
```

The results are summarized in the following figures.

``` r
# local weights variance
w_jm_empirical <- mcmc_hans_post_reorder$omega_J_M_output
w_jm_variance_empirical <- lapply(1:length(w_jm_empirical),
                                  function(t) matrix(apply(w_jm_empirical[[t]], 1, var),
                                                     nrow = 1))
w_jm_variance_empirical <- do.call(rbind, w_jm_variance_empirical)
w_jm_variance_empirical <- colMeans(w_jm_variance_empirical)
local_weights_analysis_vc$variance_empirical <- w_jm_variance_empirical

local_weights_analysis_vc %>%
  ggplot(mapping = aes(x = variance_empirical, y = probability))+
  geom_point()+
  geom_text_repel(aes(label = cluster))+
  theme_bw()+
  xlab('variance of local weights')+
  ylab('probability of observing larger variance')+
  geom_hline(yintercept = 0.95)
```

<img src="man/figures/README-variable motif plot-1.png" width="60%" style="display: block; margin: auto;" />

``` r

# global weights
local_weights_analysis_vc$global_weight = apply(matrix(unlist(mcmc_hans_post_reorder$omega_output), 
                                                       length(unique(unlist(hans_Z_reordered))), 
                                                       length(mcmc_hans_post_reorder$omega_output)),1,mean)
local_weights_analysis_vc$probability_global = prominent_motifs$prob

local_weights_analysis_vc %>%
  ggplot(mapping = aes(x = global_weight, y = probability_global))+
  geom_point()+
  geom_text_repel(aes(label = cluster))+
  theme_bw()+
  xlab('mean of the global weight')+
  ylab('probability of global weight>0.02')+
  geom_hline(yintercept = 0.95) +
  geom_vline(xintercept = 0.02)
```

<img src="man/figures/README-variable motif plot-2.png" width="60%" style="display: block; margin: auto;" />

## Total variation distance between mice

We quantify differences between mice based on the total variation (TV)
distance between brain-specific mixtures.

``` r
## ---- plot a heatmap for the posterior mean (and optionally labelled by credible intervals) of the tv_distance ----
plot_tv_distance(mcmc_run_all_output = mcmc_all_hans, hpd = TRUE, prob = 0.95, text_size = 4)
```

<img src="man/figures/README-TV distance-1.png" width="80%" style="display: block; margin: auto;" />

## Summary

Overall we can summarize the estimated projection strength, global and
local weights as well as posterior expected projection strength.

``` r
# ------ Summary plot ---------
# posterior expected number of counts/N for each region 
eps <- lapply(1:M, function(m) post_epstrength(m,mcmc_all_hans))

# summarize the posterior mean of q, omega, omega_JM (reordered)
params_summ <- list(proj_prob = mcmc_hans_post_reorder$proj_prob_mean, 
                    omega_JM = Reduce('+', mcmc_hans_post_reorder$omega_J_M_output)/length(mcmc_hans_post_reorder$omega_J_M_output),
                    omega = colMeans(do.call(rbind, mcmc_hans_post_reorder$omega_output)))

plot_summary(params = params_summ, eps = eps, prominent_motifs_prob = prominent_motifs$prob, 
             prominent_motifs_thresh = 0.95, global_weight_thresh = 0.01,
             data.source.label = 'Mouse', regions.name = rownames(data_Hans[[1]]), 
             col_bar = c("deepskyblue","darkgrey","darkturquoise","aquamarine"), col_mat = c('white','blue'), 
             legend = TRUE, legend_x = 0.6)
```

<img src="man/figures/README-summary plot-1.png" width="80%" style="display: block; margin: auto;" />

# Posterior predictive checks

For posterior predictive checks, we generate replicated data from the
posterior and compare them with the true observed data.

``` r
## --- a single replicate -----
# compare empirical projection strengths for each mouse
set.seed(3)
ppc_single_result <- ppc_single(mcmc_run_all_output = mcmc_all_hans,
                                Y = data_Hans,
                                regions.name = rownames(data_Hans[[1]]))
m=2
print(ppc_single_result[[m]])
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<img src="man/figures/README-ppc-1.png" width="70%" style="display: block; margin: auto;" />

``` r

## ----- multiple replicate: compare number of zero counts for each region (barplot) ------
### ------ compare distribution of non-zero counts for each region (boxplot) -----------

ppc_multiple_result <- ppc_multiple(mcmc_run_all_output = mcmc_all_hans,
                                   Y = data_Hans,
                                   N = 3,
                                   regions.name = rownames(data_Hans[[1]]))
ppc_multiple_result$zero.plot
```

<img src="man/figures/README-ppc-2.png" width="70%" style="display: block; margin: auto;" />

``` r
ppc_multiple_result$non.zero.plot
```

<img src="man/figures/README-ppc-3.png" width="70%" style="display: block; margin: auto;" />

# Generating data with dissection noise

To generate synthetic data with a given level of dissection noise,

``` r
set.seed(5)
Y_sim <- data_simulation(mcmc_run_all_output = mcmc_all_hans, Y = data_Hans, 
                         regions.name = rownames(data_Hans[[1]]), 
                         M = 4, C = C, noise.levels = list(0.01,0.05))
```
