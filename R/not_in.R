#' Not in
#'
#' Opposite of base \code{%in%}, converting all \code{TRUE}s in \code{FALSE}s.
#'
#' @param x First R object of any kind.
#' @param y Second R object of any kind.
#'
#' @return Logical.
#' @export
#'
#' @examples
#' x <- 1:10
#' y <- 11:20
#' x %!in% y
#'
#' @references \url{https://stackoverflow.com/a/5831829}
`%!in%` <-
  function(x, y) {
    !("%in%"(x, y))
  }
