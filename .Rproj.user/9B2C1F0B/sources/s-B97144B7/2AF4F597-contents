#' Custom per gene plotCounts
#'
#' Creates a customised ggplot2 plot of counts per gene from a \code{DESeqDataSet} object. It also add the pvalue significance to a user specified list of contrasts with \link[ggpubr]{stat_compare_means}.
#'
#' @param ensembl Character. One or many gene identifiers in the same format (i.e ensembl ID, entrez ID, symbol, etc) as in the \code{rownames(assay(dds.object))}.
#' @param dds.object DESeqDataSet. Constructed using \code{DESeq()}.
#' @param plot.title Character. Vector of the same length and order as \code{ensembl}, specifying the tile name for the plot (e.g gene name). If a vector is supplied and it does not match with \code{ensembl}, then a \code{NA} will be assigned to the plot tile. Default prints no plot title.
#' @param pvalue.contrasts Character. Each contrast must be of the form \code{c("tretament_1", "treatment_2")} and many contrasts can be especified by creating a list of contrasts. See \link[ggpubr]{stat_compare_means}.
#' @param ... other arguments to pass to \link[ggpubr]{stat_compare_means}.
#'
#' @return A list of ggplot2 plots.
#' @export
#'
#' @seealso \code{\link[DESeq2]{plotCounts}}, \code{\link[ggpubr]{stat_compare_means}}, \code{\link[ggsignif]{geom_signif}}
#'
#' @importFrom DESeq2 plotCounts
#' @importFrom gridExtra grid.arrange
#'
#' @examples
#' \dontrun{
#' genes <- c("ENSG00000026025", "ENSG00000198836")
#' gene.names <- c("vim", "opa1")
#' my.contrasts <- list(c("tibolona", "vehiculo"), c("tibolona", "tib_pal"), c("tibolona", "palmitato"), c("tibolona", "DMEM"))
#' aesir <- my_plotCounts(genes, dds, gene.names, my.contrasts)
#' aesir
#' }
#'
#' \dontrun{
#' # Supliying additional arguments to stat_compare_means()
#' genes <- c("ENSG00000026025", "ENSG00000198836")
#' gene.names <- c("vim", "opa1")
#' my.contrasts <- list(c("tibolona", "vehiculo"), c("tibolona", "tib_pal"), c("tibolona", "palmitato"), c("tibolona", "DMEM"))
#' aesir <- my_plotCounts(genes, dds, gene.names, my.contrasts, method = "t.test")
#' aesir
#' }
#'
my_plotCounts <- function(ensembl, dds.object, plot.title = NULL, pvalue.contrasts = NULL, pdf.name = NULL, ...) {
  library(ggplot2)
  library(gridExtra)

  # Checking input
  stopifnot(is.character(ensembl), class(dds.object) == "DESeqDataSet")

  # Create a list of plots
  plots.list <- lapply(1:length(ensembl), function(genei) {
    # Ensure the identifier supplied is in the expression object
    gene.name <- grep(as.character(ensembl[genei]), rownames(dds.object), value = T)

    # Extract the counts and sample_name for each gene i
    d <- DESeq2::plotCounts(dds.object, gene = gene.name, intgroup = "condition", returnData = TRUE)

    if (is.null(pvalue.contrasts)) {
      d.plot <- ggplot(d, aes(x = condition, y = count, colour = condition)) +
        geom_boxplot() + theme_classic() + theme(legend.position = "none", axis.title.x = element_blank()) + geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(1), aes(fill = condition)) + ggtitle(plot.title[genei])
    } else { # Create the plot with ggplot2 if a list of contrast is specified and add the correspondings pvalues
      d.plot <- ggplot(d, aes(x = condition, y = count, colour = condition)) +
        geom_boxplot() + theme_classic() + theme(legend.position = "none", axis.title.x = element_blank()) + geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(1), aes(fill = condition)) + ggtitle(plot.title[genei]) + ggpubr::stat_compare_means(comparisons = pvalue.contrasts, ...)
    }
    return(d.plot)
  })

  if (is.null(pdf.name)){
    return(plots.list)
  } else {
    # Chech there is a character path
    stopifnot(is.character(pdf.name))
    # For each entry (plot) in the list create a page in a pdf
    pdf(pdf.name, onefile = TRUE)
    print(plots.list)
    dev.off()
    return(message("\nPlease check if your pdf was created.\n"))
  }


  return(message("All done!"))
}
