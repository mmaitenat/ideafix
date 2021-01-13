#' filter_variants_XGBoost
#'
#' @param variant_descriptors Matrix containing the descriptors of the variants to be analyzed and subjected to filtering
#'
#' @details Note that the function only works if \strong{variant_descriptors} includes all the descriptors that were used to train the model
#' @return Matrix containing three columns: variant id, probability of deamination and predicted label: deamination or non-deamination
#' @export
#' @import dplyr
#' @examples
filter_variants_XGBoost <- function(variant_descriptors) {
  XGBoost_model_filename <- system.file("extdata", "XGBoost_final_model.RDS", package = "ideafix")
  XGBoost_model <- readRDS(XGBoost_model_filename)
  variant_descriptors_XGBoost <- one_hot_encoding(variant_descriptors)
  XGBoost_pred <- stats::predict(XGBoost_model,
                                 as.matrix(variant_descriptors_XGBoost),
                                 type = "prob")
  thr <- 0.5766285
  XGBoost_pred_class <- ifelse(XGBoost_pred$X1 > thr, "deamination", "non-deamination") %>%
    factor(levels = c("non-deamination", "deamination"))
  predictions <- tibble(CHROM = gsub(":.*", "", variant_descriptors$id),
                        POS = gsub(".*:", "", variant_descriptors$id),
                        REF = variant_descriptors$ref.allele,
                        ALT = variant_descriptors$alt.allele,
                        DEAM_SCORE = XGBoost_pred$X1,
                        DEAMINATION = XGBoost_pred_class)
  return(predictions)
}

#' filter_variants_RF
#'
#' @param variant_descriptors Matrix containing the descriptors of the variants to be analyzed and subjected to filtering
#'
#' @return Matrix containing three columns: variant id, probability of deamination and predicted label: deamination or non-deamination
#'
#' @export
#' @details As opposed to filter_variants_XGBoost, filter_variants_RF can run even if not all the descriptors in the model are included in \strong{variant_descriptors}. However, good performance is not guaranteed.
#' @import h2o
#' @import dplyr
#' @examples

filter_variants_RF <- function(variant_descriptors) {
  RF_model_parts <- list.files(system.file("extdata", package = "ideafix"), "RF_final_model.tar.gz.parta*", full.names = TRUE)
  join_uncompress_targz(RF_model_parts)
  RF_model_filename <- system.file("extdata", "RF_final_model", package = "ideafix")
  Sys.unsetenv("http_proxy")
  h2o.init(
    port = 54325,
    nthreads = 2,
    max_mem_size = "8G")
  RF_model <- h2o.loadModel(RF_model_filename)
  variant_descriptors_RF <- dplyr::select(variant_descriptors, -ends_with("id"))
  variant_descriptors_h2o <- as.h2o(variant_descriptors_RF)
  RF_pred <- h2o.predict(object = RF_model,
                         newdata = variant_descriptors_h2o) %>%
    as_tibble()
  predictions <- tibble(CHROM = gsub(":.*", "", variant_descriptors$id),
                        POS = gsub(".*:", "", variant_descriptors$id),
                        REF = variant_descriptors$ref.allele,
                        ALT = variant_descriptors$alt.allele,
                        DEAM_SCORE = RF_pred$X1,
                        DEAMINATION = RF_pred$predict) %>%
    mutate(DEAMINATION = factor(ifelse (DEAMINATION == "X1", "deamination", "non-deamination")))
  return(predictions)
}

#' check_descriptors
#'
#' @param variant_descriptors Matrix containing the descriptors of the variants to be analyzed and subjected to filtering
#' @param algorithm One of the algorithms to use to filter the variants. Can be RF or XGBoost.
#'
#' @return None
#' @export
#'
#' @examples
check_descriptors <- function(variant_descriptors, algorithm) {
  model_descriptors <- c("allele.freq", "alt.bases", "norm.alt.bases", "ref.bases", "norm.ref.bases", "ref.allele", "alt.allele", "base.qual", "base.qual.frac", "frag.length", "pos.from.end", "map.qual", "FdeamC", "SOB", "SBGuo", "SBGATK", "norm.pos.from.end", "before.2", "before.1", "after.1", "after.2", "before", "after", "is.repeat.region")
  data_descriptors <- colnames(variant_descriptors)
  missing_data_descriptors <- data_descriptors[!grepl("*.id$", data_descriptors)] %>%
    setdiff(model_descriptors, .)
  if (length(missing_data_descriptors)) {
    if (algorithm == "RF") {
      warning(sprintf("variant_descriptors has less descriptors than expected. Missing descriptors are: %s. Filtering process will be run but performance is not guaranteed.", paste(missing_data_descriptors, collapse = ", ")))
    } else if (algorithm == "XGBoost") {
      stop(sprintf("variant_descriptors has less descriptors than expected. Missing descriptors are: %s. Filtering process can not be run.", paste(missing_data_descriptors, collapse = ", ")))
    } else {
      stop("Algorithm not available. Please choose one of the following: RF, XGBoost")
    }
  } else {
    message("Descriptor check went successfully. Proceeding with filtering.")
  }
}

#' check_algorithm
#'
#' @param algorithm One of the algorithms to use to filter the variants. Can be RF or XGBoost.
#'
#' @return None
#' @export
#'
#' @examples
check_algorithm <- function(algorithm) {
  if (algorithm %in% c("RF", "XGBoost")) {
    message(sprintf("Running variant filtering with %s", algorithm))
  } else {
    stop("Algorithm not available. Please choose one of the following: RF, XGBoost")
  }
}


#' filter_variants
#'
#' @param variant_descriptors Matrix containing the descriptors of the variants to be analyzed and subjected to filtering
#' @param algorithm One of the algorithms to use to filter the variants. Can be RF or XGBoost. Defaults to RF.
#'
#' @return Matrix containing three columns: variant id, probability of deamination and predicted label: deamination or non-deamination
#' @export
#'
#' @examples
filter_variants <- function(variant_descriptors, algorithm = "RF") {
  defined <- ls()
  passed <- c(names(as.list(match.call())[-1]), "algorithm") # manually add algorithm because it does alwalys have a default value, which is RF
  check_arguments(defined = defined, passed = passed)
  check_algorithm(algorithm = algorithm)
  check_descriptors(variant_descriptors = variant_descriptors, algorithm = algorithm)
  predictions <- if(algorithm == "RF") {
    filter_variants_RF(variant_descriptors)
  } else {
    filter_variants_XGBoost(variant_descriptors)
  }
  message("Filtering went successfully.")
  return(predictions)
}

#' Join split tar files and uncompress them
#'
#' @param split_targz_files Vector containing the full paths of the parts in which a tar file has been split
#'
#' @return None
#' @export
#'
#' @examples
join_uncompress_targz <- function(split_targz_files) {
  message("Uncompressing RF model...")
  outdir <- system.file("extdata", package = "ideafix")
  cmd <- sprintf("cat %s | tar -xzvf - -C %s", paste(split_targz_files, collapse = " "), outdir)
  system(cmd, intern = FALSE)
}
