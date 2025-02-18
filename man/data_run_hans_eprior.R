
# Hans data - also called the VC data
# Make sure data_Hans_5 and data_Hans are loaded before running the code below
load("data/Han-data/data_Hans.RData")
load("data/Han-data/data_Hans_5.RData")
# data_Hans = lapply(data_Hans,round)

# Source required functions
source("R/main.R")
source("R/mcmc_run_all_replusion.R")
source("R/mcmc_steps_repulsion.R")

# for cpp functions
source("huizi/mcmc_run_all_repulsion.R")
source("huizi/mcmc_steps_cpp.R")
Rcpp::sourceCpp("huizi/cpp_funcs.cpp")

# this is the function from sara's mcmc_run_all_repulsion.R, but remove split and merge, and not print anything
# used for time fair time comparisons.
source("huizi/mcmc_run_all.R")

M <- length(data_Hans)
C <- sapply(1:M, function(m) ncol(data_Hans[[m]]))
R <- dim(data_Hans[[1]])[1]

# Gel-plot
gel_plot_data_Hans <- gel_plot(Y = data_Hans)


png(file = './plots/Hans/gelplot_hans.png',
    width = 2000,
    height = 500)

ggarrange(gel_plot_data_Hans[[1]],
          gel_plot_data_Hans[[2]],
          gel_plot_data_Hans[[3]],
          gel_plot_data_Hans[[4]],
          nrow = 1,
          widths = c(1,1,1,1.4))

dev.off()


# Mouse index
mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]),
                 rep(4, C[4]))

#---------------------------------------------------------------------------------------


#Initialize the clustering and set hyperparameters empirically from the data
data_Hans_cbind <- do.call(cbind, data_Hans)

C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_Hans[[m]]))))

# Option 2: initialize with k-means with cosine distance
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

# Run with repulsive prior, remove merge and split, not print info.
library(profvis) # to see time complexity of each step
profvis({

set.seed(3)
# run1 <- mcmc_run_all_cpp(Y = data_Hans,
#                          auto.save = TRUE,
#                          save_frequency = 1000,
#                          partial.save.name = 'huizi/run1.RData',
mcmc_all_hans <- mcmc_run_all_cpp(Y = data_Hans,
                                  J = J,
                              number_iter = 6000,
                              thinning = 1,
                              burn_in = 5000,
                              adaptive_prop = 0.0001,
                              print_Z = FALSE,
                              # iter_update = 100,
                              a_gamma = a_gamma,
                              b_gamma = b_gamma,
                              a_alpha =  a_alpha,
                              b_alpha= 1,
                              a_alpha0 = a_alpha0,
                              b_alpha0= 1,
                              Z.init = clustinit)
})
# run1: 1.455 mins
mcmc_all_hans = run1

#----------------------------------------------------------------------------------------
                              
mcmc_all_hans$output_index
ind <- seq(1,1001,by=1)

Zmat = matrix(unlist(mcmc_all_hans$Z_output), length(mcmc_all_hans$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')

# Trace plots - checking mixing/convergence
plot(mcmc_all_hans$alpha_zero_output, type = 'l')
plot(mcmc_all_hans$alpha_output, type = 'l')
plot(unlist(mcmc_all_hans$acceptance_prob$q_star)[ind], type = 'l')
plot(unlist(mcmc_all_hans$acceptance_prob$gamma_star)[ind], type = 'l')
plot(unlist(mcmc_all_hans$acceptance_prob$alpha)[ind], type = 'l')
plot(unlist(mcmc_all_hans$acceptance_prob$alpha_zero)[ind], type = 'l')
j=5
r=4
plot(unlist(lapply(mcmc_all_hans$q_star_1_J_output[ind], function(q){q[j,r]})), type = 'l')
j=29
m=4
plot(unlist(lapply(mcmc_all_hans$omega_J_M_output[ind], function(q){q[j,m]})), type = 'l')
j=1
plot(unlist(lapply(mcmc_all_hans$omega_output[ind], function(q){q[j]})), type = 'l')
j=5
plot(unlist(lapply(mcmc_all_hans$gamma_star_1_J_output[ind], function(q){q[j]})), type = 'l')

# Posterior similarity matrix
psm_hans = similarity_matrix(mcmc_run_all_output = mcmc_all_hans)


# Reordered posterior samples of z
hans_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_hans)


# optimal clustering
hans_Z <- opt.clustering.comb(z_trace = hans_z_reordered,
                              post_similarity = psm_hans,
                              max.k = max(k))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))

hans_Z <- lapply(1:M,
                function(m) hans_Z[(C_cumsum[m]+1):C_cumsum[m+1]])

#-------------------------------------------------------------------------------------------


# Plot of posterior similarity matrix
psm_hans_plot <- plotpsm(psm.ind = psm_hans$psm.within,
                         psm.tot = psm_hans$psm.combined)


# TODO: Change with superheat to add mouse labels
png(file = './plots/Hans_eprior/heatmap_psm(1).png',
    width = 500,
    height = 400)

psm_hans_plot$plot.ind

dev.off()

#--------------------------------------------  
# MCMC unique
mcmc_unique_hans <- mcmc_run_post(mcmc_run_all_output = mcmc_all_hans,
                                  Z = hans_Z,
                                  thinning = 5,
                                  burn_in = 2000,
                                  number_iter = 12000,
                                  Y = data_Hans,
                                  a_gamma = a_gamma,
                                  b_gamma = b_gamma,
                                  regions.name = rownames(data_Hans[[1]]))

                                 
# Neuron allocations
hans_Z_reordered <- mcmc_unique_hans$Z


# Number of neurons by cluster and mosue
png(file = './plots/Hans_eprior/number_of_neuron_by_m.png',
    width = 600,
    height = 300)

opt.clustering.frequency(clustering = hans_Z_reordered)

dev.off()

# Heatmap of projection strength of neurons in each cluster, colored by mouse
png(file = './plots/Hans_eprior/heatmap_neuron.png',
    width = 400,
    height = 400)

# ------- rename the function !!!!  heatmap.ps -----------
pp.standard.ordering2(Y = data_Hans,
                      Z = hans_Z_reordered,
                      regions.name = rownames(data_Hans[[1]]),
                      mouse.index = mouse.index)

dev.off()

# ---------------------------!!!!! Line plots for of empirical projection strengths within each cluster ------------------------

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

bayesian_motif_pp <- lapply(c(1:mcmc_unique_hans$J),
                            function(j){
                              
                              data.frame(cluster = paste(j, ':', paste(colnames(mcmc_unique_hans$q_tilde_001)[mcmc_unique_hans$q_tilde_001[j,] >= 0.5], collapse = ',')),
                                         pp = as.vector(t(df[which(unlist(mcmc_unique_hans$Z) == j),])),
                                         region.name = rownames(data_Hans[[1]]),
                                         cell.index = rep(which(unlist(mcmc_unique_hans$Z) == j), each = 6),
                                         mouse.index = rep(mouse.index[which(unlist(mcmc_unique_hans$Z) == j)], each = 6))
                            })


bayesian_motif_pp <- do.call(rbind, bayesian_motif_pp)

jpeg(file = './plots/Hans_eprior/line_plots.jpeg',
    width = 650,
    height = 400,
    quality = 100)

bayesian_motif_pp %>%
  mutate(cell.index = as.factor(cell.index),
         mouse.index = as.factor(mouse.index),
         cluster = factor(cluster,
                          levels = unique(bayesian_motif_pp$cluster)),
         region.name = factor(region.name, levels = rownames(data_Hans[[1]]))) %>%
  ggplot(mapping = aes(x = region.name,
                       y = pp,
                       color = mouse.index,
                       group = cell.index))+
  geom_line()+
  facet_wrap(~cluster, ncol = 5)+
  theme_bw()+
  xlab('region')+
  ylab('projection strength')+
  guides(color = guide_legend(title="mouse"))+
  ggtitle('HBMAP')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

#-------------------------------------------- 
# ---------- !!!! Line Plot of estimated projection strength within each cluster ----------

## Add: plot those line graphs color-coded by the number of projecting regions
number_of_projecting_regions <- lapply(sort(unique(unlist(mcmc_unique_hans$Z))),
                                       function(j){
                                         
                                         return(data.frame(cluster = paste('cluster', j),
                                                           number_of_p = length(which(mcmc_unique_hans$q_tilde_001[j,] >= 0.5))))
                                       })

number_of_projecting_regions <- do.call(rbind, number_of_projecting_regions)

plot_df <- mcmc_unique_hans$estimated.projection.df %>%
  left_join(number_of_projecting_regions)

plot_df$number_of_p <- factor(plot_df$number_of_p,
                              levels = 1:R)

# All clusters
regions.name = regions.name = rownames(data_Hans[[1]])

png(file = './plots/Hans_eprior/estimated_pp.png',
    width = 1200,
    height = 800)

ggplot(plot_df,mapping = aes(x = factor(region, levels = regions.name),
                             y = projection.med,
                             group = cluster,
                             color = number_of_p)) +
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin = projection.lower,
                    ymax = projection.upper),
                width = 0.1)+
  theme_bw()+
  ylim(c(0,1))+
  ylab('projection strength')+
  xlab('region')+
  labs(color = 'regions') + 
  facet_wrap(~cluster)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(size=10))

dev.off()

#-------------------------------------------- 

# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_hans,
                                   mcmc_run_post_output = mcmc_unique_hans,
                                   J = mcmc_unique_hans$J,
                                   thinning = 5,
                                   burn_in = 5000,
                                   number_iter = 15000)


# For all clusters, compare motifs weights across mice
png(file = './plots/Hans_eprior/w_jm.png',
    width = 1200,
    height = 800)

omega_JM_mcmc$omega_JM_plot

dev.off()

# ------- !!!!!! Find prominent motifs: Probability of the global weight greater than a threshold ------- 
thresh = 0.02
print(paste('With a threshold of', thresh, 'we identify motifs where we would expect at least', thresh*sum(C), 'neurons in that motif across all mice'))
omega_mat = matrix(unlist(omega_JM_mcmc$omega_output), 
                     nrow = length(omega_JM_mcmc$omega_output), 
                     ncol = mcmc_unique_hans$J,
                     byrow = TRUE)
prob_greater_global = apply(omega_mat>thresh,2, mean)
# ---------------------

# ---------- !!!! find variable motifs: Variance of local weights ---------- 
w_jm_empirical <- omega_JM_mcmc$omega_J_M_output
w_jm_variance_empirical <- lapply(1:length(w_jm_empirical),
                                  function(t) matrix(apply(w_jm_empirical[[t]], 1, var),
                                                     nrow = 1))
w_jm_variance_empirical <- do.call(rbind, w_jm_variance_empirical)
w_jm_variance_empirical <- colMeans(w_jm_variance_empirical)

local_weights_analysis_vc <- local_weights_analysis(N = 200,
                                                    mcmc_run_all_output = mcmc_all_hans,
                                                    mcmc_run_post_output = mcmc_unique_hans,
                                                    mcmc_run_omega_output = omega_JM_mcmc)

local_weights_analysis_vc$variance_empirical <- w_jm_variance_empirical

png(file = './plots/Hans_eprior/local_weight_var.png',
    width = 450,
    height = 400)

local_weights_analysis_vc %>%
  ggplot(mapping = aes(x = variance_empirical, y = probability))+
  geom_point()+
  geom_text_repel(aes(label = cluster))+
  theme_bw()+
  xlab('variance of local weights')+
  ylab('probability of observing larger variance')+
  geom_hline(yintercept = 0.95)

dev.off()

local_weights_analysis_vc$global_weight = apply(matrix(unlist(omega_JM_mcmc$omega_output), mcmc_unique_hans$J, length(omega_JM_mcmc$omega_output)),1,mean)
local_weights_analysis_vc$probability_global = prob_greater_global

png(file = './plots/Hans_eprior/global_weight.png',
    width = 450,
    height = 400)

local_weights_analysis_vc %>%
  ggplot(mapping = aes(x = global_weight, y = probability_global))+
  geom_point()+
  geom_text_repel(aes(label = cluster))+
  theme_bw()+
  xlab('mean of the global weight')+
  ylab('probability of global weight>0.02')+
  geom_hline(yintercept = 0.95) +
  geom_vline(xintercept = 0.02)

dev.off()


# ------- ??? Add to demo: Use superheat to add cluster bars to psm -----------
library(superheat)
png(file = './plots/Hans_eprior/psm_by_z.png',
    width = 500,
    height = 450)

Z_factored = as.character(unlist(hans_Z_reordered))
ind = sapply(unlist(hans_Z_reordered), function(z){
  sum(z == prominent_motifs)>0
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

dev.off()

#-------------------------------------------- 
# -------- !!! Line plot: Repeat plot projection strength, color coded with uncertainty but only for prominent motifs ---------


## Add: plot those line graphs color-coded by the number of projecting regions
plot_df <- mcmc_unique_hans$estimated.projection.df %>%
  left_join(number_of_projecting_regions)

plot_df$number_of_p <- factor(plot_df$number_of_p,
                              levels = 1:R)

# ------ !!!! Find prominent motifs: filter based on prominent motifs ----------- 
prominent_motifs = which(prob_greater_global>0.95)

ind = sapply(plot_df$cluster, function(c){c %in% paste('cluster', prominent_motifs)})
plot_df = plot_df[ind,]

# All clusters
regions.name = regions.name = rownames(data_Hans[[1]])

png(file = './plots/Hans_eprior/estimated_pp_prominent.png',
    width = 600,
    height = 300)

ggplot(plot_df,mapping = aes(x = factor(region, levels = regions.name),
                             y = projection.med,
                             group = cluster,
                             color = number_of_p)) +
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin = projection.lower,
                    ymax = projection.upper),
                width = 0.1)+
  theme_bw()+
  ylim(c(0,1))+
  ylab('projection strength')+
  xlab('region')+
  labs(color = 'regions') + 
  facet_wrap(~cluster)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(size=10))

dev.off()

#---------------------------!!!! Line plots of empirical projection strengths for prominent motifs ------------------------

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

bayesian_motif_pp <- lapply(prominent_motifs,
                            function(j){
                              
                              data.frame(cluster = paste(j, ':', paste(colnames(mcmc_unique_hans$q_tilde_001)[mcmc_unique_hans$q_tilde_001[j,] >= 0.5], collapse = ',')),
                                         pp = as.vector(t(df[which(unlist(mcmc_unique_hans$Z) == j),])),
                                         region.name = rownames(data_Hans[[1]]),
                                         cell.index = rep(which(unlist(mcmc_unique_hans$Z) == j), each = 6),
                                         mouse.index = rep(mouse.index[which(unlist(mcmc_unique_hans$Z) == j)], each = 6))
                            })


bayesian_motif_pp <- do.call(rbind, bayesian_motif_pp)

jpeg(file = './plots/Hans_eprior/line_plots_prominent.jpeg',
     width = 650,
     height = 400,
     quality = 100)

bayesian_motif_pp %>%
  mutate(cell.index = as.factor(cell.index),
         mouse.index = as.factor(mouse.index),
         cluster = factor(cluster,
                          levels = unique(bayesian_motif_pp$cluster)),
         region.name = factor(region.name, levels = rownames(data_Hans[[1]]))) %>%
  ggplot(mapping = aes(x = region.name,
                       y = pp,
                       color = mouse.index,
                       group = cell.index))+
  geom_line()+
  facet_wrap(~cluster, ncol = 3)+
  theme_bw()+
  xlab('region')+
  ylab('projection strength')+
  guides(color = guide_legend(title="mouse"))+
  ggtitle('HBMAP')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()


#-------------------------------------------- 
# Investigate if differences in total projection counts across clusters
# Distribution of N_{i,m} for neurons within the same cluster

data_Hans_N <- lapply(1:length(data_Hans),
                      function(m) colSums(data_Hans[[m]]))

df <- data.frame(N = unlist(data_Hans_N),
                 motif = unlist(mcmc_unique_hans$Z))

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_hans$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')



#-------------------------------------------- 
# ------------- !!! ppc  Model checking: Posterior predictive check with multiple replicates ------------
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_hans,
                      Y = data_Hans,
                      N = 3,
                      regions.name = rownames(data_Hans[[1]]))


png(file = './plots/Hans_eprior/ppc_zero.png',
    width = 600,
    height = 300)

ppc_multiple$zero.plot

dev.off()


png(file = './plots/Hans_eprior/ppc_nonzero.png',
    width = 600,
    height = 300)

ppc_multiple$non.zero.plot

dev.off()



# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_hans,
                           Y = data_Hans,
                           regions.name = rownames(data_Hans[[1]]))


for(m in 1:M){
  
  png(file = paste0('./plots/Hans_eprior/ppc_single_mouse_', m, '.png'),
      width = 800,
      height = 400)
  
  print(ppc_single[[m]])
  
  dev.off()
}

#-------------------------------------------- 
# Posterior expected empirical projection strength

# Given a fixed number of total counts across all regions
total_counts = 100
# Maybe we want to use something that reflects the typical number of total counts across neurons
# total_counts =median(unlist(lapply(data_Hans,colSums)))

#Mouse 1
ps_m1 = post_epstrength(1,total_counts,mcmc_all_hans)

#Mouse 2
ps_m2 = post_epstrength(2,total_counts,mcmc_all_hans)

#Mouse 3
ps_m3 = post_epstrength(3,total_counts,mcmc_all_hans)

#Mouse 4
ps_m4 = post_epstrength(4,total_counts,mcmc_all_hans)

# Projection strength across mice

group.colors <- c('1'= "deepskyblue",'2'= "darkgrey",'3'= "darkturquoise", '4' = "aquamarine")
df = data.frame(p = c(ps_m1$eps,ps_m2$eps,ps_m3$eps,ps_m4$eps))
df$Region = factor(rep(rownames(data_Hans[[1]]), M), levels = rownames(data_Hans[[1]]))
df$Mouse = c(rep('1', R),rep('2', R),rep('3', R),rep('4',R))
df$lower = c(ps_m1$eps_ci[1,],ps_m2$eps_ci[1,],ps_m3$eps_ci[1,],ps_m4$eps_ci[1,])
df$upper = c(ps_m1$eps_ci[2,],ps_m2$eps_ci[2,],ps_m3$eps_ci[2,],ps_m4$eps_ci[2,])

png(file = './plots/Hans_eprior/ps_mouse.png',
    width = 600,
    height = 500)

ggplot(df, mapping = aes(x = Region, y= p, color = Mouse, group = Mouse)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_color_manual(values=group.colors) +
  ylab('E[y/N]') +
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

################################
## Heat plot to summarize HBMAP results
################################
## TODO: create a function to create this summary plot

library(ComplexHeatmap)
library(circlize)
library(colorRamp2)

## Replace with results of the data!!!
mcmc_unique = mcmc_unique_hans #MCMC results for motif parameters given optimal clustering
mcmc_omega = omega_JM_mcmc #MCMC results for weights given optimal clustering

# Matrix of expected projection strength across clusters (regions x clusters)
df = data.frame(matrix(mcmc_unique$estimated.projection.df$projection.mean, R, mcmc_unique$J),row.names = mcmc_unique$estimated.projection.df$region[1:R])
names(df) = seq(1:mcmc_unique$J)
# Matrix of expected weight across clusters for the global and local weights for each mouse  (clusters x M+1)
weightmat = matrix(cbind(mcmc_omega$mean_omega_J, mcmc_omega$mean_omega_J_M), nrow = mcmc_unique$J, ncol = M +1)
# Matrix of expected projection strength for each mouse (regions x M)
psmat = matrix(cbind(ps_m1$eps,ps_m2$eps,ps_m3$eps,ps_m4$eps), nrow = R, ncol = M)
# vector of probability that the global weight is greater than a threshold (length = num of clusters)
pg = prob_greater_global

# Colorscale for heatmap
col_fun = colorRamp2(c(0, 1), c( "white", "blue"))

# Top annotation: expected global and local weights for each mouse
# Note: if more mice, add another color. Red indicates global weights
colnames(weightmat) = c("global", paste("mouse", c(1:M)))
column_ha  = HeatmapAnnotation(w = anno_barplot(weightmat, 
                                                beside = TRUE, attach = TRUE, 
                                                barwidth = 0.99,
                                                gp = gpar(fill = c("red","deepskyblue","darkgrey","darkturquoise","aquamarine")),
                                                height = unit(6, "cm"))
)
# If you only want the global weights in the top annotation
#column_ha = HeatmapAnnotation(global = anno_barplot(mcmc_omega$mean_omega_J))

# Left annotation: expected projection strength for each mouse 
row_ha  = rowAnnotation(p = anno_barplot(psmat, 
                                         beside = TRUE, attach = TRUE, 
                                         barwidth = 0.99,
                                         gp = gpar(fill = c("deepskyblue","darkgrey","darkturquoise","aquamarine")),
                                         border = FALSE,
                                         width = unit(2, "cm"))
)


# Save as size 1500 by 750

# Heat map of projection strength for each cluster with annotation
Heatmap(df,
        cluster_rows = FALSE, 
        row_names_side = "left",
        cluster_columns = FALSE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title = ""),
        col = col_fun,
        column_order = order(pg, decreasing = TRUE),
        column_split = as.factor(2*(pg<0.95&weightmat[,1]<0.02) + (pg<0.95&weightmat[,1]>0.02)),
        column_title_gp = gpar(fontsize=0),
        border = TRUE,
        top_annotation = column_ha,
        left_annotation = row_ha
)

lgd = Legend(labels = c("Global","Mouse 1","Mouse 2","Mouse 3","Mouse 4"), title = "", 
             legend_gp = gpar(fill = c("red","deepskyblue","darkgrey","darkturquoise","aquamarine")),
             border = "black")

draw(lgd,x = unit(0.9, "npc"), y = unit(0.95, "npc"),just = c("right", "top"))



#---------------------------??? To do: Uncertainty in the neuron allocation ------------------------


# Probability in allocation
allo.prob = allocation_probability(mcmc_unique_hans,
                                  mcmc_omega,
                                  data_Hans,
                                  mcmc_unique_hans$Z)
projection_probability(data_Hans,
                       mcmc_unique_hans$Z,
                       rownames(data_Hans[[1]]),
                       allo.prob$allocation_probability,
                       1:mcmc_unique_hans$J)

# ---------------------------??? To do: TV distance between mice ------------------------


# Compute the total variation distance between the mixing measures for each mouse pair

mytv_dist = function(x,ind){
  xdim = dim(mcmc_all_hans$omega_J_M_output[[1]])
  y = matrix(x[,ind],xdim[1],xdim[2])
  return(0.5*colSums(abs(x-y)))
}

m =  1 # index of mouse
tv_dist = lapply(mcmc_all_hans$omega_J_M_output,mytv_dist, ind = m )
tv_dist <- data.frame(matrix(unlist(tv_dist), nrow=length(tv_dist), byrow=TRUE))
names(tv_dist) = c("1", "2", "3", "4")
tv_dist = tv_dist[-m]

tv_dist = tv_dist %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Mouse", 
    values_to = "TV"
  )


group.colors <- c('1'= "deepskyblue",'2'= "darkgrey",'3'= "darkturquoise", '4' = "aquamarine")


png(file = './plots/Hans_eprior/mouse1_tv.png',
    width = 350,
    height = 300)

ggplot(tv_dist) +
  geom_density(aes(x=TV,col=Mouse, fill=Mouse), alpha=0.25) +
  ggtitle(paste('Mouse',m,'comparisons')) + 
  theme_bw()+
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors)

dev.off()



m =  2 # index of mouse
tv_dist = lapply(mcmc_all_hans$omega_J_M_output,mytv_dist, ind = m )
tv_dist <- data.frame(matrix(unlist(tv_dist), nrow=length(tv_dist), byrow=TRUE))
names(tv_dist) = c("1", "2", "3", '4')
tv_dist = tv_dist[-m]

tv_dist = tv_dist %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Mouse", 
    values_to = "TV"
  )




png(file = './plots/Hans_eprior/mouse2_tv.png',
    width = 350,
    height = 300)

ggplot(tv_dist) +
  geom_density(aes(x=TV,col=Mouse, fill=Mouse), alpha=0.25) +
  ggtitle(paste('Mouse',m,'comparisons')) + 
  theme_bw()+
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors)

dev.off()



m =  3 # index of mouse
tv_dist = lapply(mcmc_all_hans$omega_J_M_output,mytv_dist, ind = m )
tv_dist <- data.frame(matrix(unlist(tv_dist), nrow=length(tv_dist), byrow=TRUE))
names(tv_dist) = c("1", "2", "3", '4')
tv_dist = tv_dist[-m]

tv_dist = tv_dist %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Mouse", 
    values_to = "TV"
  )



png(file = './plots/Hans_eprior/mouse3_tv.png',
    width = 350,
    height = 300)

ggplot(tv_dist) +
  geom_density(aes(x=TV,col=Mouse, fill=Mouse), alpha=0.25) +
  ggtitle(paste('Mouse',m,'comparisons')) + 
  theme_bw()+
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors)

dev.off()

m =  4 # index of mouse
tv_dist = lapply(mcmc_all_hans$omega_J_M_output,mytv_dist, ind = m )
tv_dist <- data.frame(matrix(unlist(tv_dist), nrow=length(tv_dist), byrow=TRUE))
names(tv_dist) = c("1", "2", "3", '4')
tv_dist = tv_dist[-m]

tv_dist = tv_dist %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Mouse", 
    values_to = "TV"
  )



png(file = './plots/Hans_eprior/mouse4_tv.png',
    width = 350,
    height = 300)

ggplot(tv_dist) +
  geom_density(aes(x=TV,col=Mouse, fill=Mouse), alpha=0.25) +
  ggtitle(paste('Mouse',m,'comparisons')) + 
  theme_bw()+
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors)

dev.off()


# Compute the posterior mean
tv_mean = matrix(0,4,4)
for (m in c(1:4)){
  tv_dist_m = lapply(mcmc_all_hans$omega_J_M_output,mytv_dist, ind = m )
  tv_dist_m = data.frame(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
  tv_mean[m,] = colMeans(tv_dist_m)
}
tv_mean = data.frame(tv_mean, row.names = c("1", "2", "3", '4'))
names(tv_mean) = c("1", "2", "3", '4')
print(tv_mean)

tv_mean = data.frame(tv_mean, "Mouse 1" = c("1", "2", "3", '4'))
names(tv_mean)[1:4] = c("1", "2", "3", '4')
tv_mean <-  pivot_longer(tv_mean,
                         cols = !Mouse.1,
                         names_to = "Mouse.2", 
                         values_to = "TV"
)

ggplot(tv_mean, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
  geom_tile() +
  labs(x = "Mouse",
       y = "Mouse") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "orange", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(tv_mean$TV,3)), color = "black", size = 4) +
  coord_fixed()

##### Compute HPD interval
# Compute HPD intervals
M <- 4
tv_lower = matrix(0, mcmc_all_hans$M,mcmc_all_hans$M)
tv_upper = matrix(0, mcmc_all_hans$M,mcmc_all_hans$M)
for (m in c(1:mcmc_all_hans$M)){
  tv_dist_m = lapply(mcmc_all_hans$omega_J_M_output,mytv_dist, ind = m )
  tv_dist_m = as.mcmc(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
  tv_hpd =  HPDinterval((tv_dist_m))
  tv_lower[m,] = tv_hpd[,1]
  tv_upper[m,] = tv_hpd[,2]
}

tv_lower = data.frame(tv_lower, row.names = as.factor(c(1:mcmc_all_hans$M)))
names(tv_lower) = as.factor(c(1:mcmc_all_hans$M))
print(tv_lower)

tv_lower = data.frame(tv_lower, "Mouse 1" = as.factor(c(1:mcmc_all_hans$M)))
names(tv_lower)[1:mcmc_all_hans$M] = as.factor(c(1:mcmc_all_hans$M))
tv_lower <-  pivot_longer(tv_lower,
                          cols = !Mouse.1,
                          names_to = "Mouse.2", 
                          values_to = "TV"
)

tv_upper = data.frame(tv_upper, row.names = as.factor(c(1:mcmc_all_hans$M)))
names(tv_upper) = as.factor(c(1:mcmc_all_hans$M))
print(tv_upper)

tv_upper = data.frame(tv_upper, "Mouse 1" = as.factor(c(1:mcmc_all_hans$M)))
names(tv_upper)[1:mcmc_all_hans$M] = as.factor(c(1:mcmc_all_hans$M))
tv_upper <-  pivot_longer(tv_upper,
                          cols = !Mouse.1,
                          names_to = "Mouse.2", 
                          values_to = "TV"
)

labs = paste0(round(tv_mean$TV,3),' [',round(tv_lower$TV,2),',', round(tv_upper$TV,2),']')


png(file = './plots/Hans_eprior/total_variation.png',
    width = 750,
    height = 500)


ggplot(tv_mean, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
  geom_tile() +
  labs(x = "Mouse",
       y = "Mouse") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "orange", low="white", midpoint = 0.5) +
  geom_text(aes(label = labs), color = "black", size = 4) +
  coord_fixed()

dev.off()


