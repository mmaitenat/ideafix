#' Check required arguments exist
#'
#' Checks if required arguments for a function are present
#'
#' It stops the program if the argument is not present.
#'
#' @param argument Name of the argument to be checked
#'
#' @return None
#' @export
#'
#' @examples
check_arguments <- function(defined, passed) {
  if (any(!defined %in% passed)) {
    stop(paste("Missing values for", paste(setdiff(defined, passed), collapse = ", ")))
  }
}
