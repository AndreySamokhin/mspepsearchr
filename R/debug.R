#==============================================================================#
#' Call MSPepSearch from R
#'
#' @description
#'   Call MSPepSearch from R for testing and debugging purposes.
#'
#' @param ...
#'   Strings or character vectors. Arguments passed to the MSPepSearch tool.
#' @param os_arch
#'   A string. The architecture of the operating system (i.e., \code{"x32"} or
#'   \code{"x64"}) or \code{NULL} (to detect it automatically).
#'
#' @return
#'   Return \code{NULL}.
#'
#' @noRd
#'
#==============================================================================#
.CallMsPepSearch <- function(..., os_arch = NULL) {

  temp <- list(...)
  if (!all(vapply(temp, is.character, logical(1L)))) {
    stop("'...' must be a set of character vectors.")
  }
  args <- unlist(temp)

  exe_path <- .GetExePath(os_arch)
  if (!file.exists(exe_path)) {
    stop("'", basename(exe_path), "' is missing.")
  }

  # err_code <- system2(command = exe_path, args = args, stderr = NULL)
  err_code <- system2(command = exe_path, args = args)
  cat(paste0("error code: ", err_code))

  return(invisible(NULL))
}


