% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allocation_prob.R
\name{allocation_probability}
\alias{allocation_probability}
\title{Compute average allocation probabilities across all MCMC samples}
\usage{
allocation_probability(Y, post_output_reorder)
}
\arguments{
\item{Y}{a list of matrices. Each is a region-by-neuron matrix of projection counts for individual mouse.}

\item{post_output_reorder}{output from \code{mcmc_reorder_cluster}.}
}
\value{
a list of numeric vectors. Each vector represents the probabilities of neurons in a mouse belonging to their optimal cluster.
}
\description{
Compute average allocation probabilities across all MCMC samples
}
