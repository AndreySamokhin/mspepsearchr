#==============================================================================#
#' Split spectra into multiple files for parallel computation
#'
#' @description
#'   A helper function that writes mass spectra to disk and split the output
#'   into several files when \code{n_treads > 1L}.
#'
#' @return A character vector. Names of the temporary output files.
#'
#' @importFrom mssearchr WriteMsp
#'
#' @noRd
#==============================================================================#
.SplitSpectraForParallel <- function(spectra,
                                     temp_dir,
                                     n_threads) {

  #--[ Check input arguments ]--------------------------------------------------

  # 'spectra'
  if (is.character(spectra)) {
    if (length(spectra) != 1L) {
      stop("'spectra' must be a single file path.")
    }
    # Other sanity checks are performed within 'normalizePath()'.
  } else {
    if (!is.list(spectra) || !is.list(spectra[[1]]) ||
        is.null(spectra[[1]]$mz) || is.null(spectra[[1]]$intst) ||
        is.null(spectra[[1]]$name)) {
      stop("'spectra' is invalid.")
    }
  }

  # 'temp_dir'
  if (!is.character(temp_dir) || length(temp_dir) != 1) {
    stop("'temp_dir' must be a string.")
  }
  if (!dir.exists(temp_dir)) {
    stop("The directory ", "'", temp_dir, "'", " does not exist.")
  }

  # 'n_threads'
  .CheckNumThreads(n_threads)


  #--[ Helper functions ]-------------------------------------------------------

  .Split <- function(len, n) {
    # Split sequence into approximately equal parts and return start indices for
    # each part. Code was adopted from 'parallel::splitIndices()'.
    if (len <= n) {
      return(seq_len(len))
    }
    fuzz <- 0.4
    start_idxs <- ceiling(seq(1 - fuzz, by = (len - 1 + 2 * fuzz) / n,
                              length.out = n))
    return(start_idxs)
  }


  #--[ Split input spectra ]----------------------------------------------------

  if (is.character(spectra)) { # i.e., 'spectra' is a path to an MSP/MGF file
    if (n_threads == 1L) {
      return(normalizePath(spectra))
    } else {
      all_lines <- readLines(spectra)
      ms_start_idxs <- grep("^name:", all_lines, ignore.case = TRUE)
      if (length(ms_start_idxs) > 0L) { # i.e., MSP
        file_extension <- ".msp"
      } else { # i.e., MGF
        begin_ions_idxs <- grep("^begin ions", all_lines, ignore.case = TRUE)
        if (length(begin_ions_idxs) == 0L) {
          stop("No spectra found in file: ", spectra)
        }
        stop("Multiprocessing is not yet supported for MGF files.")
      }
      if (n_threads > length(ms_start_idxs)) {
        n_threads <- length(ms_start_idxs)
      }
      ms_split_idxs <- .Split(length(ms_start_idxs), n_threads)
      chunk_start_idxs <- ms_start_idxs[ms_split_idxs]
      chunk_end_idxs <- c(chunk_start_idxs[-1L] - 1L, length(all_lines))
      input_files <- vapply(seq_len(n_threads), function(file_no) {
        line_idxs <- seq(chunk_start_idxs[[file_no]], chunk_end_idxs[[file_no]])
        temp_file_name <- tempfile("mspepsearchr_", tmpdir = temp_dir,
                                   fileext = file_extension)
        writeLines(all_lines[line_idxs], temp_file_name)
        return(temp_file_name)
      }, character(1L))
      return(input_files)
    }
  } else { # i.e., 'spectra' is a list of nested lists
    if (n_threads == 1L) {
      input_file <- tempfile("mspepsearchr_", tmpdir = temp_dir,
                             fileext = ".msp")
      mssearchr::WriteMsp(spectra, input_file)
      return(input_file)
    } else {
      start_idxs <- .Split(length(spectra), n_threads)
      end_idxs <- c(start_idxs[-1L] - 1L, length(spectra))
      input_files <- vapply(seq_len(n_threads), function(file_no) {
        spectra_idxs <- seq(start_idxs[[file_no]], end_idxs[[file_no]])
        temp_file_name <- tempfile("mspepsearchr_", tmpdir = temp_dir,
                                   fileext = ".msp")
        mssearchr::WriteMsp(spectra[spectra_idxs], temp_file_name)
        return(temp_file_name)
      }, character(1L))
      return(input_files)
    }
  }
}



#==============================================================================#
#' Prepare CLI jobs
#'
#' @description
#'   Internal helper that assembles CLI arguments and output file paths.
#'
#' @return
#'   A list of job descriptors, each containing:
#'   \itemize{
#'     \item \code{cli_args} - A character vector of CLI arguments.
#'     \item \code{output_file} - The path to the output file returned by
#'     MSPepSearch.
#'   }
#'
#' @noRd
#==============================================================================#
.PrepareJobs <- function(
    spectra,
    libraries,
    algorithm,
    presearch,
    addl_columns,
    valid_presearch_vals,
    precursor_ion_tol = NULL,
    product_ions_tol = NULL,
    ignore_precursor_ion_tol = NULL,
    nominal_mw = NULL,
    n_hits = 100L,
    search_method = c("standard", "reverse", "pss"),
    best_hits_only = FALSE,
    penalize_rare_compounds = FALSE,
    min_abundance = 1L,
    mz_min = NULL,
    mz_max = NULL,
    ri_column_type = NULL,
    load_in_memory = FALSE,
    temp_dir = NULL,
    n_threads = 1L,
    addl_cli_args = NULL
) {

  #--[ Local helper functions ]-------------------------------------------------

  .IsToleranceValid <- function(tol) {
    if (is.null(tol)) {
      return(TRUE)
    }
    if (!is.list(tol) || length(tol) != 2L ||
        !is.numeric(tol[[1L]]) || length(tol[[1L]]) != 1L ||
        !is.character(tol[[2L]]) || length(tol[[2L]]) != 1L ||
        !(tol[[2L]] %in% c("mz", "ppm"))) {
      return(FALSE)
    }
    return(TRUE)
  }


  #--[ Reorder 'presearch$type' ]-----------------------------------------------

  if (is.list(presearch) && length(presearch) == 2L &&
      match("type", names(presearch), nomatch = 1L) != 1L) {
    presearch <- presearch[c(2L, 1L)]
  }


  #--[ Check input arguments ]--------------------------------------------------

  # 'spectra'
  # It is checked elsewhere, see '.SplitSpectraForParallel()'.


  # 'libraries'
  # Sanity checks are performed within 'normalizePath()'.
  if (length(libraries) > 16L) {
    stop("The maximum number of mass spectral libraries is 16.")
  }
  if (length(libraries) != length(unique(libraries))) {
    stop("'libraries' contains duplicate libraries.")
  }
  if (sum(tolower(basename(libraries)) == "mainlib") > 1L) {
    stop("'libraries' contains two 'mainlib' entries.")
  }
  if (sum(tolower(basename(libraries)) == "replib") > 1L) {
    stop("'libraries' contains two 'replib' entries.")
  }

  # 'algorithm'
  # It is checked elsewhere.

  # 'valid_presearch_vals'
  all_presearch_vals <-
    c("default", "default_msms", "fast", "off", "mw", "inchikey")
  if (!is.character(valid_presearch_vals) ||
      length(valid_presearch_vals) == 0L) {
    stop("'valid_presearch_vals' must be a character vector.")
  }
  if (!all(valid_presearch_vals %in% all_presearch_vals)) {
    stop("'valid_presearch_vals' is invalid.")
  }

  # 'presearch'
  if (!((is.character(presearch) && length(presearch) == 1L) ||
        (is.list(presearch) && length(presearch) <= 2L))) {
    stop("'presearch' must be a string or a list of one or two elements.")
  }
  if (!is.character(presearch[[1L]]) || length(presearch[[1L]]) != 1L ||
      !(presearch[[1L]] %in% all_presearch_vals)) {
    stop("'presearch' is invalid.")
  }
  if (!(presearch[[1L]] %in% valid_presearch_vals)) {
    stop("'presearch' value is not supported by the selected algorithm.")
  }
  if (presearch[[1L]] == "mw") {
    if (!is.null(names(presearch)) && is.na(match("mw", names(presearch)))) {
      stop("'presearch' does not contain 'mw' element.")
    }
    if (!is.numeric(presearch[[2L]]) || length(presearch[[2L]]) != 1L ||
        abs(presearch[[2L]] - as.integer(presearch[[2L]])) > 1e-10) {
      stop("'presearch$mw' must be an integer value.")
    }
  }
  if (presearch[[1L]] == "inchikey") {
    if (!is.null(names(presearch)) &&
        is.na(match("n_segments", names(presearch)))) {
      stop("'presearch' does not contain 'n_segments' element.")
    }
    if (!is.numeric(presearch[[2L]]) || length(presearch[[2L]]) != 1L ||
        !(as.integer(presearch[[2L]]) %in% c(1L, 2L, 3L))) {
      stop("'presearch$n_segments' must be 1, 2, or 3.")
    }
  }

  # 'addl_columns'
  if (!is.character(addl_columns) || length(addl_columns) == 0L) {
    stop("'addl_columns' must be a character vector.")
  }

  # 'precursor_ion_tol'
  if (!.IsToleranceValid(precursor_ion_tol)) {
    stop("'precursor_ion_tol' is invalid.")
  }

  # 'ignore_precursor_ion_tol'
  if (!.IsToleranceValid(ignore_precursor_ion_tol)) {
    stop("'ignore_precursor_ion_tol' is invalid.")
  }

  # 'product_ions_tol'
  if (!.IsToleranceValid(product_ions_tol)) {
    if (algorithm == "identity_hires") {
      stop("'mz_tolerance' is invalid.")
    } else {
      stop("'product_ions_tol' is invalid.")
    }
  }

  # 'nominal_mw'
  if (!is.null(nominal_mw)) {
    if (!is.numeric(nominal_mw) || length(nominal_mw) != 1L ||
        abs(nominal_mw - as.integer(nominal_mw)) > 1e-10) {
      stop("'nominal_mw' must be an integer value.")
    }
  }

  # 'n_hits'
  if (!is.numeric(n_hits) || length(n_hits) != 1L ||
      abs(n_hits - as.integer(n_hits)) > 1e-10 || as.integer(n_hits) < 1L) {
    stop("'n_hits' must be a non-negative integer value.")
  }
  if (n_hits > 100L) {
    stop("'n_hits' must be less than or equal to 100.")
  }

  # 'search_method'
  # The 'match.arg()' function is used.

  # 'best_hits_only'
  if (!is.logical(best_hits_only) || length(best_hits_only) != 1L) {
    stop("'best_hits_only' must be a logical value.")
  }

  # 'penalize_rare_compounds'
  if (!is.logical(penalize_rare_compounds) ||
      length(penalize_rare_compounds) != 1L) {
    stop("'penalize_rare_compounds' must be a logical value.")
  }

  # 'min_abundance'
  if (!is.numeric(min_abundance) || length(min_abundance) != 1L ||
      abs(min_abundance - as.integer(min_abundance)) > 1e-10) {
    stop("'min_abundance' must be an integer value.")
  }
  if (min_abundance < 1 || min_abundance > 999) {
    stop("'min_abundance' must be in the range [1; 999].")
  }

  # 'mz_min'
  if (!is.null(mz_min)) {
    if (!is.numeric(mz_min) || length(mz_min) != 1L ||
        abs(mz_min - as.integer(mz_min)) > 1e-10 || as.integer(mz_min) < 0L) {
      stop("'mz_min' must be a non-negative integer value.")
    }
  }

  # 'mz_max'
  if (!is.null(mz_max)) {
    if (!is.numeric(mz_max) || length(mz_max) != 1L ||
        abs(mz_max - as.integer(mz_max)) > 1e-10 || as.integer(mz_max) < 0L) {
      stop("'mz_max' must be a non-negative integer value.")
    }
    if (!is.null(mz_min)) {
      if (as.integer(mz_min) >= as.integer(mz_max)) {
        stop("'mz_max' must be greater than 'mz_min'.")
      }
    }
  }

  # 'ri_column_type'
  # It is checked elsewhere.

  # 'load_in_memory'
  if (!is.logical(load_in_memory) || length(load_in_memory) != 1L) {
    stop("'load_in_memory' must be a logical value.")
  }

  # 'temp_dir'
  # It is checked elsewhere, see '.SplitSpectraForParallel()'.

  # 'n_threads'
  # It is checked elsewhere, see '.SplitSpectraForParallel()'.

  # 'addl_cli_args'
  if (!is.null(addl_cli_args)) {
    if (!is.character(addl_cli_args) || length(addl_cli_args) == 0L) {
      stop("'addl_cli_args' must be a character vector.")
    }
  }


  #--[ CLI arguments ]----------------------------------------------------------

  # [{sdfmk[n]}][aijnopqr[2]vx][{uyz}][{leh[n]}[{IQSHLMPGD}]
  # [/Z[PPM] z] [/M[PPM] m] [/ZI[PPM] zi] [/W[PPM] w] ...
  # [/PATH path]
  # [/LIB lib]
  # [/WRK wrk]
  # /INP InputFile
  # [/OUT[TAB] OutputFile] [/HITS max_hits] [/MinMF minmf] ...

  cli_args <- list(encoded = "",
                   params  = character(0L),
                   libs    = character(0L),
                   output  = character(0L),
                   addl    = NULL)

  temp <- unlist(lapply(addl_cli_args, strsplit, split = " ", fixed = TRUE))
  if ("x" %in% temp) {
    cli_args$encoded <- paste0(cli_args$encoded, "x")
    cli_args$addl <- temp[!(temp %in% c("", "x"))]
  } else {
    cli_args$addl <- temp[temp != ""]
  }

  forbidden_options <-
    c("/PATH", "/MAIN", "/REPL", "/LIB", "/WRK", "/INP", "/OUT", "/OUTTAB")
  mask <- (forbidden_options %in% toupper(cli_args$addl))
  if (any(mask)) {
    stop("The following CLI options are not intended to be set manually: ",
         paste(forbidden_options[mask], collapse = ", "))
  }


  #--[ Input files ]------------------------------------------------------------

  if (is.null(temp_dir)) {
    temp_dir <- tempdir()
  }
  input_file_names <- .SplitSpectraForParallel(spectra = spectra,
                                               temp_dir = temp_dir,
                                               n_threads = n_threads)


  #--[ CLI arguments (search options) ]-----------------------------------------

  #.. Algorithm ................................................................
  # LoRes search types:
  #   I - identity
  #   Q - quick identity
  #   S - simple similarity
  #   H - hybrid similarity
  #   L - neutral loss similarity
  #   M - ms/ms in EI library
  # HiRes search type:
  #   P - peptide (use peak annotations and peak weighting)- Default
  #   G - generic (no peak annotations or peptide peak weighting used)
  #   D - the score is a pure dot product
  # HiRes options:
  #   u - don't match precursor m/z (this is HiRes No Precursor spectrum search)
  #   z - ms/ms spectrum search (match precursor m/z within precursor ion m/z
  #       uncertainty /Z[PPM])
  #   y - hybrid similarity search for ms/ms spectra
  #   a - alternative peak matching in computing match factors (recommended)
  #   l/e/h - low, medium or high ms/ms search threshold
  encoded_algorithm <- switch(algorithm,
                              identity_normal = "I",
                              identity_quick = "Q",
                              similarity_simple = "S",
                              similarity_hybrid = "H",
                              similarity_msms_hybrid = "aylG",
                              similarity_neutral_loss = "L",
                              similarity_msms_in_ei = "M",
                              identity_hires = "aulG",
                              identity_msms = "azlG",
                              stop("'algorithm' is invalid."))
  cli_args$encoded <- paste0(cli_args$encoded, encoded_algorithm)

  #.. Presearch ................................................................
  encoded_presearch <- switch(presearch[[1L]],
                              default = "d",
                              default_msms = "m",
                              fast = "f",
                              off = "s",
                              mw = paste("/MwPresearch", presearch[[2L]]),
                              inchikey = paste0("k", presearch[[2L]]),
                              stop("'presearch' is invalid."))
  if (startsWith(encoded_presearch, "/")) {
    cli_args$params <- c(cli_args$params, encoded_presearch)
  } else {
    cli_args$encoded <- paste0(cli_args$encoded, encoded_presearch)
  }

  #.. Search method ............................................................
  search_method <- match.arg(search_method)
  if (search_method == "reverse") {
    cli_args$encoded <- paste0(cli_args$encoded, "r")
  } else if (search_method == "pss") {
    cli_args$encoded <- paste0(cli_args$encoded, "r2")
  }

  #.. Penalize rare compounds ..................................................
  if (penalize_rare_compounds) {
    # This option is deprecated in MS Search 3.0 (NIST23).
    cli_args$encoded <- paste0(cli_args$encoded, "p")
  }

  #.. Precursor ion tolerance, /Z or /ZPPM .....................................
  if (!is.null(precursor_ion_tol)) {
    if (precursor_ion_tol[[2L]] == "mz") {
      cli_args$params <- c(cli_args$params, "/Z", precursor_ion_tol[[1L]])
    } else { # i.e., 'precursor_ion_tol[[2L]] == "ppm"'
      cli_args$params <- c(cli_args$params, "/ZPPM", precursor_ion_tol[[1L]])
    }
  }

  #..Ignoring peaks around precursor, /ZI or /ZIPPM ............................
  if (!is.null(ignore_precursor_ion_tol)) {
    cli_args$encoded <- paste0(cli_args$encoded, "i")
    if (ignore_precursor_ion_tol[[2L]] == "mz") {
      cli_args$params <-
        c(cli_args$params, "/ZI", ignore_precursor_ion_tol[[1L]])
    } else { # i.e., 'ignore_precursor_ion_tol[[2L]] == "ppm"'
      cli_args$params <-
        c(cli_args$params, "/ZIPPM", ignore_precursor_ion_tol[[1L]])
    }
  }

  #.. Product ion tolerance, /M or /MPPM .......................................
  if (!is.null(product_ions_tol)) {
    if (product_ions_tol[[2L]] == "mz") {
      cli_args$params <- c(cli_args$params, "/M", product_ions_tol[[1L]])
    } else { # i.e., 'product_ions_tol[[2L]] == "ppm"'
      cli_args$params <- c(cli_args$params, "/MPPM", product_ions_tol[[1L]])
    }
  }

  #.. Nominal MW or precursor m/z, /MwForLoss ..................................
  if (!is.null(nominal_mw)) {
    cli_args$params <- c(cli_args$params, "/MwForLoss", nominal_mw)
  }

  #.. Minimum abundance, /MinInt ...............................................
  if (min_abundance > 1L) {
    cli_args$params <- c(cli_args$params, "/MinInt", min_abundance)
  }

  #.. Minimum and maximum m/z, /MzLimits .......................................
  if (!is.null(mz_min) || !is.null(mz_max)) {
    if (is.null(mz_min)) {
      mz_min <- -1L
    }
    if (is.null(mz_max)) {
      mz_max <- -1L
    }
    cli_args$params <- c(cli_args$params, "/MzLimits", mz_min, mz_max)
  }


  #--[ CLI arguments (mass spectral libraries) ]--------------------------------

  libraries <- normalizePath(libraries)

  #.. MainLib, /MAIN ...........................................................
  mainlib_idx <- which(tolower(basename(libraries)) == "mainlib")
  if (length(mainlib_idx) == 1L) {
    cli_args$libs <- c(cli_args$libs, "/MAIN", libraries[[mainlib_idx]])
  } else if (length(mainlib_idx) > 1L) {
    stop("Several 'mainlib' libraries are specified.")
  }

  #.. RepLib, /REPL ............................................................
  replib_idx <- which(tolower(basename(libraries)) == "replib")
  if (length(replib_idx) == 1L) {
    cli_args$libs <- c(cli_args$libs, "/REPL", libraries[[replib_idx]])
  } else if (length(replib_idx) > 1L) {
    stop("Several 'replib' libraries are specified.")
  }

  #.. Other libraries, /LIB ....................................................
  lib_idxs <- which(!(tolower(basename(libraries)) %in% c("mainlib", "replib")))
  for (i in lib_idxs) {
    cli_args$libs <- c(cli_args$libs, "/LIB", libraries[[i]])
  }


  #--[ CLI arguments (Output and Other Options) ]-------------------------------

  #.. Best hits only ...........................................................
  if (best_hits_only) {
    cli_args$output <- c(cli_args$output, "/OutBestHitsOnly")
  }

  #.. The number of hits, /HITS ................................................
  cli_args$output <- c(cli_args$output, "/HITS", n_hits)

  #.. Output all hit lists, /All ...............................................
  cli_args$output <- c(cli_args$output, "/All")

  #.. Retention indices, /RI ...................................................
  if ("/RI" %in% cli_args$addl) {
    warning("'ri_column_type' is ignored, ",
            "because 'addl_cli_args' contains '/RI'.")
  } else if (!is.null(ri_column_type)) {
    encoded_column <- switch(ri_column_type,
                             n = "n",
                             s = "s",
                             p = "p",
                             stdnp = "n",
                             semistdnp = "s",
                             stdpolar = "p",
                             stop("'ri_column_type' is invalid."))
    cli_args$output <- c(cli_args$output, "/RI", paste0(encoded_column, "ux"))
  }

  #.. Additional columns .......................................................
  # "/OutRevMF"    - Reverse Match Factor (the same as 'v')
  # "/OutRRevMF"   - PSS Match Factor
  # "/OutMW"       - Nominal MW
  # "/OutCAS"      - CAS registry number
  # "/OutChemForm" - Chemical formula
  # "/OutIK"       - InChIKey
  # "x"            - Do not output hit probabilities
  cli_args$output <- c(cli_args$output, addl_columns)


  #.. Load libraries in memory .................................................
  if (load_in_memory) {
    cli_args$output <- c(cli_args$output, "/LibInMem")
  }


  #--[ Check 'cli_args' ]-------------------------------------------------------

  all_cli_args <-
    c(cli_args$params, cli_args$libs, cli_args$output, cli_args$addl)
  mask <- startsWith(all_cli_args, "/")
  all_cli_options <- all_cli_args[mask]
  duplicated_cli_options <- unique(all_cli_options[duplicated(all_cli_options)])
  if (length(duplicated_cli_options) != 0L) {
    stop("The following CLI flags are duplicated: ",
         paste(duplicated_cli_options, collapse = ", "))
  }


  #--[ Compile CLI call ]-------------------------------------------------------

  out <- lapply(input_file_names, function(input_file) {
    output_file <- normalizePath(
      tempfile("mspepsearchr_", tmpdir = temp_dir, fileext = ".txt"),
      mustWork = FALSE
    )
    # TODO: rename
    cli_args <- c(cli_args$encoded,
                  cli_args$params,
                  cli_args$libs,
                  "/INP", input_file,
                  "/OUTTAB", output_file,
                  cli_args$output,
                  cli_args$addl)
    return(list(cli_args = cli_args, output_file = output_file))
  })
  return(out)
}



#==============================================================================#
#' Execute prepared CLI jobs and parse hit lists
#'
#' @description
#'   Internal helper that runs a list of prepared jobs. For each job, it calls
#'   the external CLI tool, parses the output files, and combines the results
#'   into a list of data frames.
#'
#' @param jobs A list of jobs returned by \code{.PrepareJobs()}.
#' @param n_threads An integer value. Number of threads to use.
#'
#' @return
#'   Return a list of data frames. Each data frame represents a hit list.
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#'
#' @noRd
#==============================================================================#
.ExecuteJobs <- function(jobs,
                         n_threads) {

  #--[ Check input arguments ]--------------------------------------------------

  # 'jobs'
  # It is checked elsewhere, see '.CallAndParse()'.

  # 'n_threads'
  .CheckNumThreads(n_threads)


  #--[ Execute jobs ]-----------------------------------------------------------

  if (n_threads == 1L) {
    res <- lapply(jobs, .CallAndParse)
  } else if (length(jobs) == n_threads) {
    cl <- parallel::makeCluster(n_threads)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c(".CallAndParse", ".ParseOutputFile"),
                            envir = asNamespace("mspepsearchr"))
    res <- parallel::parLapply(cl, jobs, .CallAndParse)
  } else {
    stop("The number of jobs is not equal to the number of threads.")
  }


  #--[ Output ]-----------------------------------------------------------------

  out <- unlist(res, recursive = FALSE)
  attr(out, "cli") <- lapply(res, attr, which = "cli")
  return(out)
}



#==============================================================================#
#' Execute a single CLI job and parse the output
#'
#' @description
#'   Internal helper function that runs the external CLI tool for a single job,
#'   reads the output file, and parses it into a list of data frames (i.e.,
#'   hit lists).
#'
#' @param job
#'   A single job descriptor, typically returned by \code{.PrepareJobs()}.
#'
#' @return
#'   A data frame containing the parsed results for the job.
#'
#' @noRd
#==============================================================================#
.CallAndParse <- function(job, os_arch = NULL) {

  #--[ Check input arguments ]--------------------------------------------------

  if (is.null(job$cli_args) | !is.character(job$cli_args) |
      is.null(job$output_file) | !is.character(job$output_file) |
      length(job$output_file) != 1L) {
    stop("'job' is invalid.")
  }


  #--[ Call MSPepSearch ]-------------------------------------------------------

  exe_path <- .GetExePath(os_arch)
  if (!file.exists(exe_path)) {
    stop("'", basename(exe_path), "' is missing.")
  }

  if (.Platform$OS.type == "windows") {
    err_code <- system2(exe_path, args = job$cli_args, stderr = NULL)
  } else if (.Platform$OS.type == "unix") {
    err_code <- system2("wine", args = c(exe_path, job$cli_args), stderr = NULL)
  } else {
    stop("Only Windows or Unix-like OS (with Wine installed) are supported.")
  }

  if (err_code != 0L) {
    stop("'MSPepSearch' returned non-zero exit status: ", err_code, ".")
  }


  #--[ Output ]-----------------------------------------------------------------

  hitlists <- .ParseOutputFile(job$output_file)
  return(invisible(hitlists))
}


