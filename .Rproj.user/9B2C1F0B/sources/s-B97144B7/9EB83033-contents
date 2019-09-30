#' Vulcano plot for each contrast in DESeq2
#'
#' After doing the differential expression step, the next step is to prioritise genes to validate or focus on. This function enables the quick ploting of those genes inusing vulcano plots. Additionally, the arguments \code{alpha} and \code{lfcThreshold}for the function \link[DESeq2]{results} can be tunned to be more liberal or restrictive when choosing interesting DE genes.
#'
#' @param dds_path DESeqDataSet. Obtained from running \link[DESeq2]{DESeq}. The object itself can be suplied or a quoted string defining the path to the \code{dds.rds} object (e.g. \code{my/folder/dds.rds}).
#' @param contrast_name list of \code{contrast} as described in \code{\link[DESeq2:results]{results}}. Level names can be seen by calling \code{resultsNames(dds)}.
#' @param my_lfcThreshold integer. Decision rule cutoff (e.g. DESeq2 is 0.10).The smaller the more restrictive the threshold will be and fewer genes will pop up.
#' @param my_alpha integer. Cutoff for the Log2 Fold Change (e.g. DESeq2 default is 0, however 1.5 is also very common). The bigger the more restrictive the threshold will be and fewer genes will pop up.
#' @param ... Additional arguments to the ploting function \link[EnhancedVolcano]{EnhancedVolcano}. A coomon one is to define a hard lim for the x axis by \code{xlim = c(-5, 5)}.
#'
#' @return Either a list of ggplot \code{EnhancedVolcano} plots or a message indicating the path to the plots. In the former, there will be a \code{ggplot} object per contrast.
#' @export
#'
#' @importFrom DESeq2 results
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom dplyr filter
#'
#' @examples
#'\dontrun{
#'dds <- makeExampleDESeqDataSet()
#'dds <- DESeq(dds)
#'
#'df <- as.data.frame(colData(dds))
#'}
#' data(dds)

#' contrast_name <- list(c("condition","palmitato", "vehiculo"), c("condition","tib_pal", "vehiculo"), c("condition", "tibolona", "tib_pal"), c("condition","tibolona", "vehiculo"), c("condition", "palmitato", "tibolona"))
#'
#'  vulcano_per_contrast <- my_vulcano(
#'  dds_path = dds,
#'  contrast_name = contrast_name,
#'  my_alpha = 0.05,
#'  my_lfcThreshold = 1
#' )
#'
#'
my_vulcano <-
  function(dds_path,
             contrast_name,
             my_lfcThreshold = 0,
             my_alpha = 0.1,
             output_pdf_path,
             ...) {

    # Reading DESeqDataSet.
    if (exists("dds_path")) {
      my_dds <- dds_path
    } else {
      my_dds <- readRDS(dds_path)
    }

    # https://support.bioconductor.org/p/86788/ on why to set alpha and lfcThreshold in the results table
    vulcano_plot_per_contrast <-
      lapply(1:length(contrast_name), function(i) {
        # Extract the results table

        res <-
          DESeq2::results(
            my_dds,
            contrast = contrast_name[[i]],
            lfcThreshold = my_lfcThreshold,
            alpha = my_alpha
          )


        all_DE_genes <-
          res %>%
          as.data.frame() %>%
          dplyr::filter(
            log2FoldChange >= my_lfcThreshold |
              log2FoldChange <= -my_lfcThreshold,
            padj <= my_alpha
          )
        n_DE <- nrow(all_DE_genes)

        Contrast_name <-
          paste0(contrast_name[[i]][2], "_", contrast_name[[i]][3])


        if (n_DE == 0) {
          print(
            paste0(
              "Your LFC threshold and alpha are too restrictive and no DE genes were found for the contrast ",
              Contrast_name,
              ". Try setting a more liberal combination of them."
            )
          )
        } else {
          vulcano_plot <-
            EnhancedVolcano::EnhancedVolcano(
              res,
              lab = rownames(res),
              x = "log2FoldChange",
              y = "padj",
              FCcutoff = my_lfcThreshold,
              pCutoff = my_alpha,
              title = Contrast_name,
              subtitle = paste0(
                "Number of DE genes: ",
                n_DE,
                ". Alpha: ",
                my_alpha,
                ". LFC treshold: ",
                my_lfcThreshold
              ),
              ...
            )
        }
        return(vulcano_plot)
      })

    if (!is.null(output_pdf_path)) {
      file_name <-
        paste0(
          output_pdf_path,
          "vulcano_plots_alpha",
          my_alpha,
          "_lfcThreshold",
          my_lfcThreshold,
          ".pdf"
        )
      pdf(file_name)
      invisible(lapply(vulcano_plot_per_contrast, print))
      dev.off()


      return(print(paste0(
        "The file has been saved to: ",
        file_name
      )))
    } else {
      return(vulcano_plot_per_contrast)
    }
  }
