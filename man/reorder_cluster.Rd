% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_reorder.R
\name{reorder_cluster}
\alias{reorder_cluster}
\title{Reorder the clusters and relabel the allocations based on an estimate of projection strengths}
\usage{
reorder_cluster(proj_strength, Z)
}
\arguments{
\item{proj_strength}{a \eqn{J \times R} matrix of estimated projection strengths.}

\item{Z}{a list of allocations (integers). Each is a vector of allocations for individual mouse.}
}
\value{
The function returns the reordered original labels \code{label.reorder} and the allocations after relabelling \code{Z.relabel}.
}
\description{
Reorder the clusters and relabel the allocations based on an estimate of projection strengths
}
