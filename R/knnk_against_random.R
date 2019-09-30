#' Compare average degree of the neighbors-KNN against simulated networks
#'
#' Calculates average degree of the neighbors-KNN of a given network and then perform \code{n.sims} degree-preserving simulations of the input network and calculates its corresponsing KNN. This is useful to check the validity of a network, as if it is different from the random expectation, then there must be a biologically plausible reason behind generating the observed network KNN
#'
#' @param x igraph newtork.
#' @param n.sim integer. Number of simulated networks to create.
#' @param n.cor integer. Number of processing cores.
#'
#' @return Data frame whose first column is the degree range (\eqn{k_{min}, \dots, k_{max}}) of x, the second column is the observed knnk and the rest are knnks of random simulations.
#' @export
#' @import igraph parallel
#'
knnk_against_random <-
function(x, n.sim, n.cor) {
  set.seed(0) # To ensure reproducible results from random samples

  if (is.igraph(x) != T) {
    stop("Network is not igraph class")
  }

  # Simplifying
  x <- simplify(x)

  # Extracting the number of vertices from the original network
  verte.x <- vcount(x)

  #The null model is the Erdos-Renyi random graph
  KNNK.simulated <- mclapply(1:n.sim,
    function(i) {
      knnk_random <- knn(rewire(x,
        with = keeping_degseq(niter = verte.x * 10)))$knnk
        knnk_random_df <- data.frame(k=1:length(knnk_random), Random = knnk_random)
        return(knnk_random_df)
      },
      mc.cores = n.cor
    )

  #Initialise the df with the observed degree and observed knnk
  df <- list(data.frame(k=1:max(degree(x)), Observed = knn(x)$knnk))

  KNNK.simulated <- c(df, KNNK.simulated)

  KNNK.simulated <- plyr::join_all(KNNK.simulated, by='k', type='left')

  names(KNNK.simulated) <- c("k", "Observed", 1:(length(KNNK.simulated)-2))

  return(KNNK.simulated)
  }
