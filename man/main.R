##---------- Load all required packages --------------------

suppressPackageStartupMessages({
  
  library(mvtnorm)
  library(extraDistr)
  library(SciViews)
  library(doParallel)
  library(foreach)
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
  library(fields)
  library(Matrix)
  library(Rtsne)
  library(ggpubr)
  library(ggrepel)
  library(parallel)
  library(progress)
  library(pheatmap)
  library(RColorBrewer)
  library(readxl)
  library(truncnorm)
  library(clevr)
  library(pbapply)
  library(grid)
  library(gridExtra)
  library(coda)
  library(mcclust)
  library(ClusterR)
  library(latex2exp)
  library(pracma)
  library(cascsim)
  
  # summary plot
  library(ComplexHeatmap)
  library(colorRamp2)
})

`%notin%` <- Negate(`%in%`)

# for cpp functions and combined mcmc functions
source("huizi/HPMAP_mcmc.R")
source("huizi/mcmc_steps_cpp.R")
Rcpp::sourceCpp("huizi/cpp_funcs.cpp")


# Posterior similarity matrix and plot
source('R/similarity_matrix.R')
source('R/posterior_similarity_plot.R')


# Optimal clustering
source('R/optimal_clustering.R')
source('R/clustering_estimate.R')
source('huizi/cluster_size.R')

# Post-processing: reorder cluster labels based on q and rearrange samples q, gamma, omega
source("huizi/cluster_reorder.R")


# Empirical projection strength by cluster
source('huizi/heatmap_ps.R')


# Line plot of empirical projection strengths within each cluster
# Line plot of estimated projection strengths for each cluster
source('huizi/line_plot.R')

# Identify prominent motifs
source('huizi/prominent_motif.R')
# Line plot of empirical projection strengths within each cluster colored by allocation probability
source('huizi/allocation_prob.R')

# Identify variable motifs
source('huizi/weight_analysis.R')

# Total variation distance between mice
source('huizi/tv_distance.R')

# Summary plot of q, omega and expected number of counts/N
source('huizi/summary_plot.R')

# Generate synthetic data with noise
source('huizi/synthetic_data_generation.R')

# Posterior predictive checks
source('huizi/ppc.R')
