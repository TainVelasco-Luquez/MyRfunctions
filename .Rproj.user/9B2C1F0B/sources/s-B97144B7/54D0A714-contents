#' Compare transitivity against simulated networks
#'
#' Calculates transitivity of a given network and then perform \code{n.sims} degree-preserving simulations of the input network and calculates its corresponsing transitivity. This is useful to check the validity of a network, as if it is different from the random expectation, then there must be a biologically plausible reason behind generating the observed network transitivity.
#'
#' @param x igraph newtork.
#' @param n.sim integer. Number of simulated networks to create.
#' @param n.cor integer. Number of processing cores.
#'
#' @return List with three objects: 1. Data frame with two columns: Transitivity values and Source  (observed or simulated), where each row represents an independent simulation. 2. Numeric vector of lenght 1 with the z-score value. 3. Numeric vector of lenght 1 with the observed clustering coefficient of the input graph.
#' @export
#' @import igraph parallel
#'
#' @examples
#' library(igraph)
#' my_network <- barabasi.game(50)
#' output <- cl_against_random(my_network, 100, 2)
cl_against_random <-
  function(x, n.sim, n.cor) {
    set.seed(0) # To ensure reproducible results from random samples

    if (is.igraph(x) != T) {
      stop("Network is not igraph class")
    }

    # Extracting the number of vertices from the original network
    verte.x <- vcount(x)

    # The null model is the Erdos-Renyi random graph
    trans.simulated <- mclapply(1:n.sim,
      function(i) {
        transitivity(rewire(x,
          with = keeping_degseq(niter = verte.x * 10)
        ))
      },
      mc.cores = n.cor
    ) # Alternatively one can set niter to vcount(x) * 10 or verte.x * 10 as adviced by the author. Espinoza http://biogrid.engr.uconn.edu/REU/Reports_12/Final_Reports/Max.pdf recommend 100*e.count as the minimum

    # Coercing the list class to numeric class, and this in turn into a data frame for ggplot2 compatibility
    trans.simulated <- as.data.frame(as.numeric(unlist(trans.simulated)))
    trans.simulated$Source <- "simulated"
    names(trans.simulated) <- c("Transitivity", "Source")

    # Adding the observed transitivity to trans.simulated
    trans.x <- transitivity(x) # Observed global clustering coefficient

    # Calculating the Z-score = ((X - mean.random) / satndardeviation.random), to see how many standar deviation a given data set is from its mean, between the observed transitivity and the average transitivity from random samples. Where mean(trans.simulated) is the simulated average global clustering coefficient and sd(trans.simulated) its standard deviation
    z.score <- (trans.x - mean(trans.simulated$Transitivity)) / sd(trans.simulated$Transitivity)

    # P values for the z.zscore
    p.value.oneside <- pnorm(-abs(z.score))
    p.value.twoside <- 2*pnorm(-abs(z.score))

    trans.simulated <- rbind(trans.simulated, data.frame(Transitivity = trans.x, Source = "Observed", stringsAsFactors = F))

    # Result message
    suppressWarnings(cat("\nThe observed Clustering Coefficient (C) is:", trans.x, "\nThe simulated average Clustering Coefficient (C_s) is:", mean(trans.simulated$Transitivity), "\nThe Z-score is:", z.score, "\nWith one-sided pvalue = ", p.value.oneside, "\nand two-sided pvalue =", p.value.twoside, "\n"))

    list_to_return <- list(simulated_df = trans.simulated, zscore = z.score, observed_cl = trans.x, p_oneside = p.value.oneside, p_twoside = p.value.twoside)

    return(list_to_return)
  }
