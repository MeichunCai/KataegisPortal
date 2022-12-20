
require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package ", pkg, " must be installed to use this function.",
      call. = FALSE
    )
  }
}
