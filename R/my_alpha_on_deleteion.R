#' Calculate alpha after targeted deletion of nodes
#'
#' @param network
#' @param n.cor
#'
#' @return
#' @export
#' @import poweRlaw
#'
#' @examples
my_alpha_on_deletion <-
  function(network, n.cor) {
    require(igraph)
    # Numeric sorted vector with the degree for each vertex
    net_deg <- sort(degree(network), decreasing = T)

    # Calculate cl after eqach deletion
    alpha_per_deleted <- mclapply(1:length(net_deg),
      function(i) {
        network <- delete.vertices(network, names(net_deg[1:i]))
        deg <- degree(network)
        power_all <- fit_power_law(deg)
        power_all$alpha
      },
      mc.cores = n.cor
    )
    # Numeric conversion
    alpha_per_deleted <- as.numeric(unlist(alpha_per_deleted))

    # Calculate percent of retained nodes in each iteration
    net_n <- vcount(network)
    percent_retained <- mclapply(1:length(net_deg),
      function(i) {
        network <- delete.vertices(network, names(net_deg[1:i]))
        (vcount(network) * 100) / net_n
      },
      mc.cores = n.cor
    )
    # Numeric conversion
    percent_retained <- as.numeric(unlist(percent_retained))

    # Joinning together and exporting the dataframe
    x <- cbind.data.frame(net_deg, alpha_per_deleted, percent_retained)
    row.names(x) <- names(net_deg)
    names(x) <- c("Degree", "alpha", "PercentRetained")
    return(x)
  }
