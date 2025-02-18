% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cl_against_random.R
\name{cl_against_random}
\alias{cl_against_random}
\title{Compare transitivity against simulated networks}
\usage{
cl_against_random(x, n.sim, n.cor)
}
\arguments{
\item{x}{igraph newtork.}

\item{n.sim}{integer. Number of simulated networks to create.}

\item{n.cor}{integer. Number of processing cores.}
}
\value{
List with three objects: 1. Data frame with two columns: Transitivity values and Source  (observed or simulated), where each row represents an independent simulation. 2. Numeric vector of lenght 1 with the z-score value. 3. Numeric vector of lenght 1 with the observed clustering coefficient of the input graph.
}
\description{
Calculates transitivity of a given network and then perform \code{n.sims} degree-preserving simulations of the input network and calculates its corresponsing transitivity. This is useful to check the validity of a network, as if it is different from the random expectation, then there must be a biologically plausible reason behind generating the observed network transitivity.
}
\examples{
library(igraph)
my_network <- barabasi.game(50)
output <- cl_against_random(my_network, 100, 2)
}
