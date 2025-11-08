#==============================================================================#
#' Call MSPepSearch from R
#'
#' @description
#'   Call MSPepSearch from R for testing and debugging purposes.
#'
#' @param cli_args
#'   A character vector. Arguments passed to the MSPepSearch tool.
#' @param os_arch
#'   A string. The architecture of the operating system (i.e., \code{"x32"} or
#'   \code{"x64"}) or \code{NULL} (to detect it automatically).
#'
#' @noRd
#==============================================================================#
.CallMsPepSearch <- function(cli_args, os_arch = NULL) {

  #--[ Check input arguments ]--------------------------------------------------

  if (!is.character(cli_args) || length(cli_args) == 0L) {
    stop("'cli_args' must be a character vector.")
  }


  #--[ Call MSPepSearch ]-------------------------------------------------------

  exe_path <- .GetExePath(os_arch)
  if (!file.exists(exe_path)) {
    stop("'", basename(exe_path), "' is missing.")
  }

  if (.Platform$OS.type == "windows") {
    err_code <- system2(exe_path, args = cli_args, stderr = NULL)
  } else if (.Platform$OS.type == "unix") {
    err_code <- system2("wine", args = c(exe_path, cli_args), stderr = NULL)
  } else {
    stop("Only Windows or Unix-like OS (with Wine installed) are supported.")
  }

  if (err_code != 0L) {
    stop("'MSPepSearch' returned non-zero exit status: ", err_code, ".")
  }

  return(invisible(NULL))
}


