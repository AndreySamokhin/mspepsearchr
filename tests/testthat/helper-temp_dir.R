CreateTempDir <- function(prefix = "mspepsearchr") {
  temp_dir <- file.path(tempdir(), paste0(prefix, "_", Sys.getpid()))
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir)
  }
  return(normalizePath(temp_dir, mustWork = TRUE))
}

