#' oneHotEncoding
#'
#' @param variant_descriptors Matrix containing categorical variables
#'
#' @return Matrix with values one-hot encoded
#' @import dplyr
#' @import tidyr
#' @examples
one_hot_encoding <- function(variant_descriptors) {
  categorical_vars <- variant_descriptors %>%
    select_if(function(x) is.factor(x) | is.character(x)) %>%
    colnames()
  categorical_vars <- setdiff(categorical_vars, c("id", "complete_id"))
  if (length(categorical_vars)) {
    dedup_X_by_var <- lapply(categorical_vars, function(var) {
      variant_descriptors %>%
        mutate(yesno = 1) %>%
        spread(var, yesno, fill = 0, sep = "_")
    })
    dedup_X_by_var %>%
      purrr::reduce(.f = full_join) -> dedup_X
    # This final step discards id/complete_id column
    dedup_X <- dedup_X %>%
      select_if(is.numeric)
    return(dedup_X)
  } else {
    return(dplyr::select(variant_descriptors, -ends_with("id")))
  }
}
