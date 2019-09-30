#' Targeted network robustness
#'
#'
#' @param network
#' @param n.cor
#'
#' @return
#' @export
#' @import igraph parallel
#'
#' @examples
my_targeted_robustness <-
  function(network, n.cor) {
    # Numeric sorted vector with the degree for each vertex
    net_deg <- sort(degree(network), decreasing = T)

    # Calculate cl after eqach deletion
    cl_per_deleted <- mclapply(1:length(net_deg),
      function(i) {
        network <- delete.vertices(network, names(net_deg[1:i]))
        transitivity(network)
      },
      mc.cores = n.cor
    )
    # Numeric conversion
    cl_per_deleted <- as.numeric(unlist(cl_per_deleted))

    # Calculate LL after eqach deletion
    LCC_per_deleted <- mclapply(1:length(net_deg),
      function(i) {
        network <- delete.vertices(network, names(net_deg[1:i]))
        max(components(network)$csize)
      },
      mc.cores = n.cor
    )
    # Numeric conversion
    LCC_per_deleted <- as.numeric(unlist(LCC_per_deleted))

    # Calculate percent of removed nodes in each iteration
    net_n <- vcount(network)
    percent_removed <- mclapply(1:length(net_deg),
      function(i) {
        network <- delete.vertices(network, names(net_deg[1:i]))
        (vcount(network) * 100) / net_n
      },
      mc.cores = n.cor
    )
    # Numeric conversion
    percent_removed <- as.numeric(unlist(percent_removed))

    # Joinning together and exporting the dataframe
    x <- cbind.data.frame(net_deg, cl_per_deleted, LCC_per_deleted, percent_removed)
    row.names(x) <- names(net_deg)
    names(x) <- c("Degree", "cl", "LCC", "PercentRemoved")
    return(x)
  }
