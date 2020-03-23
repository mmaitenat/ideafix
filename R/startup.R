# Cosas que no se donde van

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ideafix")
  packageStartupMessage(paste0("version ", utils::packageVersion("ideafix")))
}

# You need the suggested package for this function: h2o and xgboost
my_fun <- function(a, b) {
  if (!requireNamespace("pkg", quietly = TRUE)) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
}
