#==============================================================================#
#' Get operating system architecture
#'
#' @return
#'   A string. \code{"x32"}, \code{"x64"}, or \code{NA_character_}
#'
#' @noRd
#'
#==============================================================================#
.GetOsArch <- function() {
  os_arch <- NA_character_
  if (.Platform$OS.type == "windows") {
    # 'Sys.getenv("PROCESSOR_ARCHITEW6432")' is ignored, so if R is running in
    # 32-bit mode, the 32-bit version of 'MSPepSearch' will be used.
    os_arch <- Sys.getenv("PROCESSOR_ARCHITECTURE")
  } else if (.Platform$OS.type == "unix") {
    if (system("which wine", ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
      stop("Wine is not installed.")
    }
    os_arch <- trimws(tryCatch(
      # 'wine cmd /c echo %PROCESSOR_ARCHITECTURE%'
      system2("wine", c("cmd", "/c", "echo", "%PROCESSOR_ARCHITECTURE%"),
              stdout = TRUE, stderr = NULL),
      error = function(e) { return(NA_character_) }
    ))
  }

  if (os_arch == "AMD64") {
    return("x64")
  } else if (os_arch == "x86") {
    return("x32")
  } else {
    stop("Unable to detect OS architecture.")
  }
}



#==============================================================================#
#' Get path to the MSPepSearch tool
#'
#' @param os_arch
#'   A string. The architecture type: \code{"x32"} or \code{"x64"}.
#'
#' @return
#'   A string. The full path to the MSPepSearch tool.
#'
#' @noRd
#'
#==============================================================================#
.GetExePath <- function(os_arch = NULL) {
  if (is.null(os_arch)) {
    os_arch <- .GetOsArch()
  }
  if (os_arch == "x64") {
    path <- system.file("bin", "2024_03_15_MSPepSearch_x64_0.9.7.1",
                       "MSPepSearch64.exe", package = "mspepsearchr")
  } else if (os_arch == "x32") {
    path <- system.file("bin", "2024_03_15_MSPepSearch_x32_0.9.7.1",
                       "MSPepSearch.exe", package = "mspepsearchr")
  } else {
    stop("'os_arch' is invalid.")
  }
  return(path)
}


