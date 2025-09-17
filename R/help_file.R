#==============================================================================#
#' Open MSPepSearch help file
#'
#' @description
#'   Open the \emph{MSPepSearch64.exe.hlp.txt} file included with the
#'   MSPepSearch (NIST) tool in the default text editor.
#'
#' @return
#'   Return \code{NULL}.
#'
#' @export
#'
#==============================================================================#
OpenHelpFile <- function() {
  help_file <-
    system.file("bin", "2024_03_15_MSPepSearch_x64_0.9.7.1",
                "MSPepSearch64.exe.hlp.txt", package = "mspepsearchr")
  if (!file.exists(help_file)) {
    stop("Help file does not exist.")
  }
  return(.OpenTxtFile(help_file))
}



#==============================================================================#
#' Open txt file
#'
#' @description
#'   Open txt file
#'
#' @return
#'   Return \code{NULL}.
#'
#' @noRd
#'
#==============================================================================#
.OpenTxtFile <- function(txt_file) {
  sys_info <- Sys.info()
  tryCatch({
    switch(sys_info[["sysname"]],
           Windows = shell.exec(txt_file),
           Linux   = system2("xdg-open", txt_file),
           Darwin  = system2("open", txt_file),
           stop("Unsupported OS: ", txt_file[["sysname"]]))
  }, error = function(e) {
    message("Text file failed to open automatically.")
  })
  return(invisible(NULL))
}


