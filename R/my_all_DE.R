#' Title
#'
#' @param dds_path
#' @param contrast_name
#' @param my_FDR
#' @param my_alpha
#'
#' @return
#' @export
#'
#' @examples
my_all_DE <- function(dds_path, contrast_name, my_FDR, my_alpha) {
  my_dds <- readRDS(dds_path)

  DE_genes_per_contrast <-
    lapply(1:length(contrast_name), function(i) {
      # Extract the results table

      res <- results(my_dds, contrast = contrast_name[[i]], lfcThreshold = my_FDR, alpha = my_alpha)
      #res <- lfcShrink(my_dds, contrast = contrast_name[[i]], type="apeglm", lfcThreshold = my_FDR, alpha = my_alpha)

      # To ease its handling, lets convert the rownames to a column
      res <-
        res %>% as.data.frame()

      # Create the unified dataframe
      all_DE_genes <-
        res %>% filter(log2FoldChange >= my_FDR,
                       log2FoldChange <= -my_FDR,
                       padj <= my_alpha) %>% select(-c(lfcSE, stat, pvalue))

      Contrast_name <- paste0(contrast_name[[i]][2], "_", contrast_name[[i]][3])
      colnames(all_DE_genes) <- paste(colnames(all_DE_genes), Contrast_name,  sep = "_")
      all_DE_genes <- all_DE_genes %>% rownames_to_column("Ensembl")

      return(all_DE_genes)
    })
  DE_genes_per_contrast <- plyr::join_all(DE_genes_per_contrast, by='Ensembl', type='left')
  return(DE_genes_per_contrast)

}

#odin <- my_all_DE(dds_path = "Data/Vehiculo2/dds.rds", contrast_name = contrast_name, my_FDR= 1, my_alpha = 0.05)
