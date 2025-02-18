# Hans data - also called the VC data
# Make sure data_Hans_5 and data_Hans are loaded before running the code below

setwd("/Users/zhanghuizi/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/data/MAPseq_analysis")
load("data/Han-data/data_Hans.RData")
load("data/Han-data/data_Hans_5.RData")
# data_Hans = lapply(data_Hans,round)

# load required packages
source('huizi/main.R')



M <- length(data_Hans)
C <- sapply(1:M, function(m) ncol(data_Hans[[m]]))
R <- dim(data_Hans[[1]])[1]

# Mouse index
mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]),
                 rep(4, C[4]))




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
ggplot() + geom_point(aes(x=seq(6,40,1),y=wss))

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





# ---- parameters to pass to the main function ------
# mcmc setup
mcmc_list = list(number_iter = 5000, thinning = 1, burn_in = 4000, adaptive_prop = 0.0001,
                 auto_save = TRUE,
                 save_path = paste0('/Users/zhanghuizi/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/data/MAPseq_analysis/huizi/my_result.RData'),
                 save_frequency = 1000
)
# prior parameters, default values will be used if not provided
prior_list = list(a_gamma = a_gamma, b_gamma = b_gamma, lb_gamma = 1, a = 2, tau = 0.4, nu = 1/20,
                  a_alpha = a_alpha, b_alpha = 1, a_alpha0 = a_alpha0, b_alpha0 = 1)



# ------- Run the full model ---------
set.seed(3)
mcmc_all_hans <- HBMAP_mcmc(Y = data_Hans, J = J, mcmc = mcmc_list, prior = prior_list, Z.init = clustinit)

# set.seed(1)
# a <- HBMAP_mcmc(Y = data_Hans, J = J, mcmc = mcmc_list, prior = prior_list, Z.init = clustinit)
# 
# identical(a,b)


# ------- MCMC check -----------
ind <- seq(1,mcmc_all_hans$output_index,by=1)

Zmat = matrix(unlist(mcmc_all_hans$Z_output), length(mcmc_all_hans$Z_output), sum(C),byrow = TRUE)

## ----- Number of occupied components -------
par(mfrow=c(1,1))
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')

## ---- Concentration parameters ------
par(mfrow=c(2,3))
plot(mcmc_all_hans$alpha_zero_output, type = 'l',ylab='',main='alpha0')
plot(mcmc_all_hans$alpha_output, type = 'l',ylab='',main='alpha')
## ----- q, gamma, weights -------
j=5
r=4
# plot(unlist(lapply(mcmc_all_hans$q_star_1_J_output[ind], function(q){q[j,r]})), type = 'l', ylab='', main=paste0('q[',j,',',r,']'))
j=5
# plot(unlist(lapply(mcmc_all_hans$gamma_star_1_J_output[ind], function(q){q[j]})), type = 'l', ylab='',main=paste0('gamma[',j,']'))
j=29
m=4
plot(unlist(lapply(mcmc_all_hans$omega_J_M_output[ind], function(q){q[j,m]})), type = 'l', ylab='', main=paste0('w[',j,',',m,']'))
j=1
plot(unlist(lapply(mcmc_all_hans$omega_output[ind], function(q){q[j]})), type = 'l', ylab='', main=paste0('w[',j,']'))


## ------ Acceptance rate ------
par(mfrow=c(2,3))
plot(unlist(mcmc_all_hans$acceptance_prob$q_star)[ind],type = 'l', ylab='', main='acceptance_q')
plot(unlist(mcmc_all_hans$acceptance_prob$gamma_star)[ind], type = 'l', ylab='', main='acceptance_gamma')
plot(unlist(mcmc_all_hans$acceptance_prob$alpha)[ind], type = 'l', ylab='', main='acceptance_alpha')
plot(unlist(mcmc_all_hans$acceptance_prob$alpha_zero)[ind], type = 'l', ylab='', main='acceptance_alpha0')
plot(unlist(mcmc_all_hans$acceptance_prob$omega)[ind], type = 'l', ylab='', main='acceptance_w')








# ---------- Posterior similarity matrix ---------------
psm_hans = similarity_matrix(mcmc_run_all_output = mcmc_all_hans)

# plot
psm_hans_plot <- plotpsm(psm.ind = psm_hans$psm.within,
                         psm.tot = psm_hans$psm.combined)

psm_hans_plot$plot.ind

# Reordered posterior samples of z
hans_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_hans)


# ---------- Optimal clustering -----------
set.seed(1)
hans_Z <- opt.clustering.comb(z_trace = hans_z_reordered,
                              post_similarity = psm_hans,
                              max.k = max(k))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))

hans_Z <- lapply(1:M,
                 function(m) hans_Z[(C_cumsum[m]+1):C_cumsum[m+1]])

# save(hans_Z,file='huizi/hans_Z.RData')
load(file='huizi/hans_Z.RData')


length(unique(unlist(hans_Z)))




# parameters to pass to the main_mcmc function, all parameters have default values except for save_path (check the function)
mcmc_list = list(number_iter = 5000, thinning = 1, burn_in = 4000, adaptive_prop = 0.0001,
                 auto_save = TRUE,
                 save_path = paste0('/Users/zhanghuizi/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Mac/Project/data/MAPseq_analysis/huizi/my_result_post.RData'),
                 save_frequency = 1000
)
post_list = list(run_omega = TRUE, run_q_gamma = TRUE)

# ------- Post-processing: Sample for q and gamma. Can sample omega and concentration parameters simultaneously or individually later on -------
set.seed(3)
mcmc_hans_post <- HBMAP_mcmc(Y = data_Hans, mcmc = mcmc_list, prior = prior_list, Z.fix = hans_Z, post = post_list)



## ----- Reorder clusters based on estimated projection strengths ---------
### ----- return reordered samples as well for q, gamma and (if sampled) omegas ------
### ----- return summary statistics for q and gamma, and classification of clusters (unicast, bicast and more or broadcast to all) ------
### ----- this based on the probability of projections strength greater than a threshold (q_tilde) -------------
mcmc_hans_post_reorder <- mcmc_reorder_cluster(post_output = mcmc_hans_post, 
                                               Z = hans_Z, regions.name = rownames(data_Hans[[1]]))


# Neuron allocations after reordering
hans_Z_reordered <- mcmc_hans_post_reorder$Z

## ------ cluster size: number of neurons in each cluster, colored by group (mouse/injection site) ---------
# group.index <- rep(c(1,2),c(200,353))
# c('pink','green','yellow','powderblue')

opt.clustering.frequency(clustering = hans_Z_reordered, group.index = mouse.index, group.name = 'mouse', 
                         title = 'Cluster size')
# for a subset of clusters
# opt.clustering.frequency(clustering = hans_Z_reordered, 
#                          group.index = mouse.index, group.name = 'mouse',
#                          cluster.index = 1:3, title = '')


## ---------- Heatmap of empirical projection strength of neurons in each cluster, colored by group (mouse/injection site) ----------

heatmap_ps(Y = data_Hans, Z = hans_Z_reordered, regions.name = rownames(data_Hans[[1]]), 
           group.index = mouse.index, group.name = 'mouse',
           cluster.index = 1:length(unique(unlist(hans_Z_reordered))), title = '')


## --------- Line plot for of empirical projection strengths within each cluster ------------------------
plot_empirical_ps(Y = data_Hans, Z = hans_Z_reordered, 
                  cluster.labels = mcmc_hans_post_reorder$cluster.labels,
                  regions.name = rownames(data_Hans[[1]]),
                  group.index = mouse.index, group.name = 'mouse',
                  cluster.index = 1:length(unique(unlist(hans_Z_reordered))),
                  title = 'Empirical projection strength', facet_ncol = 5)



# ---------- Line plot of estimated projection strength within each cluster ----------
# get summary statistics
ps_summary <- ps_summarize(post_output_reorder = mcmc_hans_post_reorder)
# plot
plot_estimated_ps(ps_summary = ps_summary, cluster.index = 1:length(unique(unlist(hans_Z_reordered))), 
                  title = 'Estimated projection strength', facet_ncol = 5, col = NULL)


# ------- Find prominent motifs: Probability of the global weight greater than a threshold ------- 
# prominent motifs are defined as those with a high probability (> 0.95)
# below the function returns the indices of prominent motifs and the posterior probability of the global weight greater than thresh
## ------ need to get posterior samples of omega from the post-processing step ----------
# if you did not sample it together with q and gamma in the previous step, run the following, remember using the reordered labels
# post_list = list(run_omega = TRUE, run_q_gamma = FALSE)
# mcmc_list$auto_save=FALSE
# set.seed(4)
# mcmc_omega <- HBMAP_mcmc(Y = NULL, mcmc = mcmc_list, prior = prior_list, Z.fix = mcmc_hans_post_reorder$Z, post = post_list)
# prominent_motifs <- identify_prominent_motif(post_output_reorder = mcmc_omega, thresh = 0.02, prob = 0.95)

prominent_motifs <- identify_prominent_motif(post_output_reorder = mcmc_hans_post_reorder, 
                                             thresh = 0.02, prob = 0.95)


# redo the lines plots of empirical and estimated projection strengths for prominent motifs only
plot_empirical_ps(Y = data_Hans, Z = hans_Z_reordered, 
                  cluster.labels = mcmc_hans_post_reorder$cluster.labels,
                  regions.name = rownames(data_Hans[[1]]),
                  group.index = mouse.index, group.name = 'mouse', 
                  cluster.index = prominent_motifs$index, 
                  title = 'Empirical projection strength (prominent motif)', 
                  facet_ncol = 4, col = NULL)
plot_estimated_ps(ps_summary = ps_summary, cluster.index = prominent_motifs$index, 
                  title = 'Estimated projection strength (prominent motif)', facet_ncol = 5, col = NULL)

# line plots of empirical projection strengths colored by allocation probabilities
allo.prob = allocation_probability(post_output_reorder = mcmc_hans_post_reorder, 
                                   Y = data_Hans)
projection_probability(Y = data_Hans,
                       post_output_reorder = mcmc_hans_post_reorder,
                       allocation_prob = allo.prob$allocation_probability,
                       motif_indices = 1:mcmc_hans_post_reorder$J)


# plot the posterior similarity matrix for prominent motifs only
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


# ------- Find variable motifs: Compare the variance of the brain-specific weights across mice with a null model ------- 
# parameters to pass to the main_mcmc function, all parameters have default values except for save_path (check the function)
mcmc_list = list(number_iter = 5000, thinning = 1, burn_in = 4000, adaptive_prop = 0.0001,
                 auto_save = FALSE,
                 save_path = NULL,
                 save_frequency = 1000
)
post_list <- list(run_omega=TRUE, run_q_gamma=FALSE)
set.seed(2)
local_weights_analysis_vc <- local_weights_analysis(N = 50, Z = mcmc_hans_post_reorder$Z, 
                                                    omega_output = mcmc_hans_post_reorder$omega_output,
                                                    omega_J_M_output = mcmc_hans_post_reorder$omega_J_M_output, 
                                                    prior = prior_list, mcmc = mcmc_list, 
                                                    post = post_list, n_cores=4, verbose=FALSE)
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

# ------ Total variation distance between mice -------
## ---- plot a heatmap for the posterior mean (and optionally labelled by credible intervals) of the tv_distance ----
plot_tv_distance(mcmc_run_all_output = mcmc_all_hans, hpd = TRUE, prob = 0.95, text_size = 4)

# ------ Summary plot ---------
# posterior expected number of counts/N for each region 
eps <- lapply(1:M, function(m) post_epstrength(m,mcmc_all_hans))

# summarize the posterior mean of q, omega, omega_JM (reordered)
params_summ <- list(proj_prob = mcmc_hans_post_reorder$proj_prob_mean, 
                    omega_JM = Reduce('+', mcmc_hans_post_reorder$omega_J_M_output)/length(mcmc_hans_post_reorder$omega_J_M_output),
                    omega = colMeans(do.call(rbind, mcmc_hans_post_reorder$omega_output)))

plot_summary(params = params_summ, eps = eps, prominent_motifs_prob = prominent_motifs$prob, global_weight_thresh = 0.01,
             data.source.label = 'Mouse', regions.name = rownames(data_Hans[[1]]), 
             col_bar = c("deepskyblue","darkgrey","darkturquoise","aquamarine"), col_mat = c('white','blue'), 
             legend = TRUE, legend_x = 0.6)


# ------ Posterior predictive checks -------
## --- a single replicate -----
# compare empirical projection strengths for each mouse
set.seed(3)
ppc_single_result <- ppc_single(mcmc_run_all_output = mcmc_all_hans,
                                Y = data_Hans,
                                regions.name = rownames(data_Hans[[1]]))
  
print(ppc_single_result[[2]])

## ----- multiple replicate: compare number of zero counts for each region and mouse (barplot) ------
### ------ compare distribution of non-zero counts for each region (boxplot) -----------

ppc_multiple_result <- ppc_multiple(mcmc_run_all_output = mcmc_all_hans,
                                   Y = data_Hans,
                                   N = 3,
                                   regions.name = rownames(data_Hans[[1]]))
ppc_multiple_result$zero.plot
ppc_multiple_result$non.zero.plot
  
# ---- Generate a synthetic dataset with give noisy levels ----
set.seed(5)
Y_sim <- data_simulation(mcmc_run_all_output = mcmc_all_hans, Y = data_Hans, 
                         regions.name = rownames(data_Hans[[1]]), 
                         M = 4, C = C, noise.levels = list(0.01))








