#' Classify variants as deaminations or non-deaminations with XGBoost
#'
#' \code{classify_variants_XGBoost} classifies a set of C:G > T:A variants into deaminations or non-deaminations based on a series of relevant variant descriptors, using an XGBoost-based model.
#'
#' @param variant_descriptors tibble containing the variants to be classified together with their values for a series of descriptors obtained by \code{\link{get_descriptors}}.
#'
#' @details \code{classify_variants_XGBoost} is able to run only when \code{variant_descriptors} contains variant values for the whole set of descriptors included in the XGBoost model, namely VAF, number of alternate bases, normalized number of alternate bases, number of reference bases, normalized number of reference bases, reference allele, alternate allele, base quality, base quality fraction, fragment length, median position from read end, mapping quality, FDeamC, SOB, SB-GUO, SB-GATK, normalized median position from read end, base two positions before, base one position before, base two positions after, base one position after, dinucleotide before and dinucleotide after. If any of these is not present in the \code{variant_descriptors} object, an error message is thrown and the variant classification is stopped. Notice that all these descriptor values are automatically retrieved from a Mutect2 vcf file using \code{\link{get_descriptors}}.
#'
#'
#' @return Tibble with six columns: CHROM, POS, REF, ALT, DEAM_SCORE, DEAMINATION. CHROM and POS identify the variant position, REF and ALT describe the reference and alternate alleles. DEAM_SCORE equals to the deamination score yielded by the selected classification algorithm (RF or XGBoost). Note that these values should not be interpreted as ordinary probabilities. DEAMINATION contains the label ideafix has assigned to the variant based on an optimized classification threshold.
#' @export
#' @import dplyr
#' @import xgboost
#' @import tibble
#' @examples
classify_variants_XGBoost <- function(variant_descriptors) {
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

#' Classify variants as deaminations or non-deaminations with Random Forest
#'
#' \code{classify_variants_RF} classifies a set of C:G > T:A variants into deaminations or non-deaminations based on a series of relevant variant descriptors, using a Random Forest-based model.
#'
#' @param variant_descriptors tibble containing the variants to be classified together with their values for a series of descriptors obtained by \code{\link{get_descriptors}}.
#'
#' @return Tibble with six columns: CHROM, POS, REF, ALT, DEAM_SCORE, DEAMINATION. CHROM and POS identify the variant position, REF and ALT describe the reference and alternate alleles. DEAM_SCORE equals to the deamination score yielded by the selected classification algorithm (RF or XGBoost). Note that these values should not be interpreted as ordinary probabilities. DEAMINATION contains the label ideafix has assigned to the variant based on an optimized classification threshold.
#'
#' @export
#' @details As opposed to classify_variants_XGBoost, \code{classify_variants_RF} can run with a \code{variant_descriptors} object with fewer descriptors than those included in the model. The whole set of descriptors includes VAF, number of alternate bases, normalized number of alternate bases, number of reference bases, normalized number of reference bases, reference allele, alternate allele, base quality, base quality fraction, fragment length, median position from read end, mapping quality, FDeamC, SOB, SB-GUO, SB-GATK, normalized median position from read end, base two positions before, base one position before, base two positions after, base one position after, dinucleotide before and dinucleotide after. However, if that is the case, good performance is not guaranteed.
#' @import h2o
#' @import dplyr
#' @import tibble
#' @examples

classify_variants_RF <- function(variant_descriptors) {
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

#' Check available descriptors
#'
#' \code{check_descriptors} checks which variant descriptors are present in the input object \code{variant_descriptors} and decides if these are sufficient to be used with the model provided in \code{algorithm}.
#'
#' @param variant_descriptors tibble containing the variants to be classified together with their values for a series of descriptors obtained by \code{\link{get_descriptors}}.
#' @param algorithm character string naming the algorithm to use to classify the variants. Can be "RF" or "XGBoost".
#'
#' @return None
#' @export
#' @details \code{check_descriptors} checks if the descriptors present in \code{variant_descriptors} are sufficient to run the model in \code{algorithm} on them. If \code{variant_descriptors} object contains fewer descriptors than those included in the models and \code{algorithm} is "XGBoost", it throws an error message warning about the impossibility of executing the classification using that algorithm, and if the classification has been called, it stops. If \code{algorithm} equals to "RF" instead, a warning message alerts about this and the absence of performance guarantees, but the classification is executed if called.
#'
#' @examples
check_descriptors <- function(variant_descriptors, algorithm) {
  model_descriptors <- c("allele.freq", "alt.bases", "norm.alt.bases", "ref.bases", "norm.ref.bases", "ref.allele", "alt.allele", "base.qual", "base.qual.frac", "frag.length", "pos.from.end", "map.qual", "FdeamC", "SOB", "SBGuo", "SBGATK", "norm.pos.from.end", "before.2", "before.1", "after.1", "after.2", "before", "after", "is.repeat.region")
  data_descriptors <- colnames(variant_descriptors)
  missing_data_descriptors <- data_descriptors[!grepl("*.id$", data_descriptors)] %>%
    setdiff(model_descriptors, .)
  if (length(missing_data_descriptors)) {
    if (tolower(algorithm) == "rf") {
      warning(sprintf("variant_descriptors has less descriptors than expected. Missing descriptors are: %s. Classification process will be run but performance is not guaranteed.", paste(missing_data_descriptors, collapse = ", ")))
    } else if (tolower(algorithm) == "xgboost") {
      stop(sprintf("variant_descriptors has less descriptors than expected. Missing descriptors are: %s. Classification process can not be run with XGBoost. Consider using RF model instead.", paste(missing_data_descriptors, collapse = ", ")))
    } else {
      stop("Algorithm not available. Please choose one of the following: RF, XGBoost")
    }
  } else {
    message("Descriptor check went successfully. Proceeding with classification.")
  }
}

#' Check classification algorithm name
#'
#' \code{check_algorithm} checks if a valid \code{algorithm} value has been provided.
#'
#' @param algorithm Character string naming the algorithm to use to classify the variants.
#'
#' @return None
#' @export
#' @details Valid \code{algorithm} values are "RF" and "XGBoost". Valid value check is case-insensitive. For instance, "RF", "rf", "Rf" and "rF" are all valid \code{algorithm} argument values.
#'
#' If an invalid \code{algorithm} value is provided, an error-message is displayed and the process is stopped.
#'
#' @examples
check_algorithm <- function(algorithm) {
  if (tolower(algorithm) %in% tolower(c("RF", "XGBoost"))) {
    message(sprintf("Running variant classification with %s...", algorithm))
  } else {
    stop("Algorithm not available. Please choose one of the following: RF, XGBoost")
  }
}

#' Join and uncompress tar.gz files
#'
#' \code{join_uncompress_targz} joins multiple split tar.gz files to produce a single uncompressed file
#'
#' @param split_targz_files character vector containing the full path values of the tar.gz files to be joint together.
#'
#' @return None
#' @export
#' @details \code{join_uncompress_targz} makes an internal system call to the program \code{tar} (\code{bsdtar}/GNU \code{tar}).
#'
#' @examples
join_uncompress_targz <- function(split_targz_files) {
  message("Uncompressing RF model...")
  outdir <- system.file("extdata", package = "ideafix")
  cmd <- sprintf("cat %s | tar -xzvf - -C %s", paste(split_targz_files, collapse = " "), outdir)
  system(cmd, intern = FALSE)
}

#' Classify variants as deaminations or non-deaminations
#'
#' \code{classify_variants} classifies a set of C:G > T:A variants into deaminations or non-deaminations based on a series of relevant variant descriptors.
#'
#' @param variant_descriptors tibble containing the variants to be classified together with their values for a series of descriptors obtained by \code{\link{get_descriptors}}.
#' @param algorithm character string naming the algorithm to use to classify the variants. Can be "RF" or "XGBoost". Defaults to "RF".
#'
#' @return Tibble with six columns: CHROM, POS, REF, ALT, DEAM_SCORE, DEAMINATION. CHROM and POS identify the variant position, REF and ALT describe the reference and alternate alleles. DEAM_SCORE equals to the deamination score yielded by the selected classification algorithm (RF or XGBoost). Note that these values should not be interpreted as ordinary probabilities. DEAMINATION contains the label ideafix has assigned to the variant based on an optimized classification threshold.
#' @export
#' @details \code{classify_variants} takes as an input \code{variant_descriptors}, which is a tibble created after calling the function \code{\link{get_descriptors}}. This tibble contains, for a collection of C:G > T:A variants, the values for a series of descriptors that have shown to be relevant for the classification of these type of variants into deaminations or non-deaminations. See the documentation of \code{\link{get_descriptors}} for more details on these descriptors.
#'
#' \code{classify_variants} also takes the name of the algorithm to be used to classify the variants with the \code{algorithm} argument. Valid values are "RF" and "XGBoost", but case is irrelevant (value is case-insensitive). If an invalid \code{algorithm} value is provided, an error-message is displayed and the process is stopped.
#'
#' @examples
classify_variants <- function(variant_descriptors, algorithm = "RF") {
  defined <- ls()
  passed <- c(names(as.list(match.call())[-1]), "algorithm") # manually add algorithm because it does alwalys have a default value, which is RF
  check_arguments(defined = defined, passed = passed)
  check_algorithm(algorithm = algorithm)
  check_descriptors(variant_descriptors = variant_descriptors, algorithm = algorithm)
  predictions <- if(tolower(algorithm) == "rf") {
    classify_variants_RF(variant_descriptors)
  } else {
    classify_variants_XGBoost(variant_descriptors)
  }
  message("Classification went successfully.")
  return(predictions)
}
