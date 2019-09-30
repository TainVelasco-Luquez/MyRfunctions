#' Extract number of up and down regulated genes per contrast
#'
#' Creates a useful summary dataframe contrast wise to compare thenumber of upregulated and downregulated genes in a bulk RNAseq experiment analysed with DESEq2. This data frame can then be easily ploted to get insight about the contrast with more changes.
#'
#' @param dds DESeqDataSet object. Can be constructed with the \code{DESeq} function.
#' @param contrasts list. Each element must be a valid contrast name or contrast formula as expected by \code{DESeq2::results()}.
#'
#' @return Data frame with 5 columns: 1. "up" the number of upregulated genes in the specified contrast, 2. "down" the number of downregulated genes in the specified contrast, 3. the percentage of upregulated genes out of the total DE genes, 4. the percentage of downregulated genes out of the total DE genes and 5. the contrast name.
#' @export
#' @importFrom DESeq2 results
#' @importFrom dplyr filter count bind_rows
#'
#' @examples
#' data(dds)
#' my_contrasts <- resultsNames(dds)
#' up_down_contrasts <- MyRfunctions::my_up_down(dds, my_contrasts)
#' up_down_contrasts <- up_down_contrasts %>% tidyr::gather("Type_genes", "No_genes", up, down) %>% mutate(Percentage = ifelse(Type_genes == "up", up_percent, down_percent)) %>% select(-c(up_percent, down_percent))
#' ggplot(up_down_contrasts, aes(x = contrast, y = No_genes, fill = Type_genes)) + geom_bar(stat = "identity", position=position_dodge()) + theme_classic() + theme(axis.text.x = element_text( angle=45, hjust=1 )) + geom_text(aes(label=paste0(round(Percentage, 0), "%")), vjust=-0.2, size=3, position = position_dodge(width = 1)) + labs(x = element_blank(), y = "Number of genes", fill = "Expression")
#'
my_up_down <- function(dds, contrasts){

  if (is.list(contrasts) != T) {
    stop("Contrasts must be an object of class list")
  }

  up_down_list <- lapply(1:length(contrasts), function(i){

    # Contrast specified by name
    if (length(contrasts[[i]]) == 1) {

      x <- results(dds, name = contrasts[[i]], tidy = T)

      y_up <- x %>% filter(log2FoldChange > 0 & padj <0.1) %>% dplyr::count()
      y_down <- x %>% filter(log2FoldChange < 0 & padj <0.1) %>% dplyr::count()

      y_up_percent <- y_up * 100/sum(x$baseMean > 0)
      y_down_percent <- y_down * 100/sum(x$baseMean > 0)

      y_contrast_name <- contrasts[[i]]

      y <- data.frame(y_up, y_down, y_up_percent, y_down_percent, y_contrast_name, stringsAsFactors = F)

      return(y)

    } else {
      # Contrast specified by "contrast"

      x <- results(dds, contrast = c(contrasts[[i]][1], contrasts[[i]][2], contrasts[[i]][3]), tidy = T)

      y_up <- x %>% filter(log2FoldChange > 0 & padj <0.1) %>% dplyr::count()
      y_down <- x %>% filter(log2FoldChange < 0 & padj <0.1) %>% dplyr::count()

      y_up_percent <- y_up * 100/sum(x$baseMean > 0)
      y_down_percent <- y_down * 100/sum(x$baseMean > 0)

      y_contrast_name <- paste(contrasts[[i]][1], contrasts[[i]][3], "vs", contrasts[[i]][2], sep = "_")

      y <- data.frame(y_up, y_down, y_up_percent, y_down_percent, y_contrast_name, stringsAsFactors = F)

      return(y)

    }

  })

  up_down_df <- dplyr::bind_rows(up_down_list)
  names(up_down_df) <- c("up", "down", "up_percent", "down_percent", "contrast")

  return(up_down_df)

}
