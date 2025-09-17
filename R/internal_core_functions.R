#==============================================================================#
#' Internal interface to MSPepSearch
#'
#' @description
#'   Manage communication with MSPepSearch (NIST), an external CLI application
#'   that performs the actual library search computations. This function is
#'   called by intermediate wrappers and not intended for direct use by package
#'   users.
#'
#' @importFrom mssearchr WriteMsp
#'
#' @noRd
#'
#==============================================================================#
.LibrarySearch <- function(
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
  if (!is.null(temp_dir)) {
    if (!is.character(temp_dir) || length(temp_dir) != 1) {
      stop("'temp_dir' must be a string.")
    }
    if (!dir.exists(temp_dir)) {
      stop("The directory ", "'", temp_dir, "'", " does not exist.")
    }
  }

  # 'addl_cli_args'
  if (!is.null(addl_cli_args)) {
    if (!is.character(addl_cli_args) || length(addl_cli_args) == 0L) {
      stop("'addl_cli_args' must be a character vector.")
    }
  }


  #--[ Input and output files ]-------------------------------------------------

  if (is.null(temp_dir)) {
    temp_dir <- tempdir()
  }

  if (is.character(spectra)) {
    input_file <- normalizePath(spectra)
  } else { # i.e., 'is.list(spectra) && ...'
    input_file <- tempfile("mspepsearchr_", tmpdir = temp_dir, fileext = ".msp")
    mssearchr::WriteMsp(spectra, input_file)
  }

  libraries <- normalizePath(libraries)

  output_file <- tempfile("mspepsearchr_", tmpdir = temp_dir, fileext = ".txt")


  #--[ CLI arguments ]----------------------------------------------------------

  cli_args <- list()

  #.. Algorithm ................................................................
  encoded_algorithm <- switch(algorithm,
                              identity_normal = "I",
                              identity_quick = "Q",
                              similarity_simple = "S",
                              similarity_hybrid = "H",
                              similarity_msms_hybrid = c("G", "y", "a", "l"),
                              similarity_neutral_loss = "L",
                              similarity_msms_in_ei = "M",
                              identity_hires = c("G", "u", "a", "l"),
                              identity_msms = c("G", "z", "a", "l"),
                              stop("'algorithm' is invalid."))
  cli_args <- c(cli_args, encoded_algorithm)
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

  #.. Presearch ................................................................
  encoded_presearch <- switch(presearch[[1L]],
                              default = "d",
                              default_msms = "m",
                              fast = "f",
                              off = "s",
                              mw = paste("/MwPresearch", presearch[[2L]]),
                              inchikey = paste0("k", presearch[[2L]]),
                              stop("'presearch' is invalid."))
  cli_args <- c(cli_args, encoded_presearch)

  #.. Input file, /INP .........................................................
  cli_args <- c(cli_args, paste("/INP", input_file))

  #.. MainLib, /MAIN ...........................................................
  idx_mainlib <- which(tolower(basename(libraries)) == "mainlib")
  if (length(idx_mainlib) > 0L) {
    cli_args <- c(cli_args, paste("/MAIN", libraries[[idx_mainlib]]))
  }

  #.. RepLib, /REPL ............................................................
  idx_replib <- which(tolower(basename(libraries)) == "replib")
  if (length(idx_replib) > 0L) {
    cli_args <- c(cli_args, paste("/REPL", libraries[[idx_replib]]))
  }

  #.. Other libraries, /LIB ....................................................
  mask_other_libs <- !(tolower(basename(libraries)) %in% c("mainlib", "replib"))
  if (sum(mask_other_libs) > 0L) {
    cli_args <- c(cli_args, paste("/LIB", libraries[mask_other_libs]))
  }

  #.. Precursor ion tolerance, /Z or /ZPPM .....................................
  if (!is.null(precursor_ion_tol)) {
    if (precursor_ion_tol[[2L]] == "mz") {
      cli_args <- c(cli_args, paste("/Z", precursor_ion_tol[[1L]]))
    } else { # i.e., 'precursor_ion_tol[[2L]] == "ppm"'
      cli_args <- c(cli_args, paste("/ZPPM", precursor_ion_tol[[1L]]))
    }
  }

  #..Ignoring peaks around precursor, /ZI[PPM]
  if (!is.null(ignore_precursor_ion_tol)) {
    cli_args <- c(cli_args, "i")
    if (ignore_precursor_ion_tol[[2L]] == "mz") {
      cli_args <- c(cli_args, paste("/ZI", ignore_precursor_ion_tol[[1L]]))
    } else { # i.e., 'ignore_precursor_ion_tol[[2L]] == "ppm"'
      cli_args <- c(cli_args, paste("/ZIPPM", ignore_precursor_ion_tol[[1L]]))
    }
  }

  #.. Product ion tolerance, /M or /MPPM .......................................
  if (!is.null(product_ions_tol)) {
    if (product_ions_tol[[2L]] == "mz") {
      cli_args <- c(cli_args, paste("/M", product_ions_tol[[1L]]))
    } else { # i.e., 'product_ions_tol[[2L]] == "ppm"'
      cli_args <- c(cli_args, paste("/MPPM", product_ions_tol[[1L]]))
    }
  }

  #.. Nominal MW or precursor m/z, /MwForLoss ..................................
  if (!is.null(nominal_mw)) {
    cli_args <- c(cli_args, paste("/MwForLoss", nominal_mw))
  }

  #.. The number of hits, /HITS ................................................
  cli_args <- c(cli_args, paste("/HITS", n_hits))

  #.. Search method ............................................................
  search_method <- match.arg(search_method)
  if (search_method == "reverse") {
    cli_args <- c(cli_args, "r")
  } else if (search_method == "pss") {
    cli_args <- c(cli_args, "r2")
  }

  #.. Best hits only ...........................................................
  if (best_hits_only) {
    cli_args <- c(cli_args, "/OutBestHitsOnly")
  }

  #.. Penalize rare compounds ..................................................
  if (penalize_rare_compounds) {
    # This option is deprecated in MS Search 3.0 (NIST23).
    cli_args <- c(cli_args, "p")
  }

  #.. Minimum abundance, /MinInt ...............................................
  if (min_abundance > 1L) {
    cli_args <- c(cli_args, paste("/MinInt", min_abundance))
  }

  #.. Minimum and maximum m/z, /MzLimits .......................................
  if (!is.null(mz_min) || !is.null(mz_max)) {
    if (is.null(mz_min)) {
      mz_min <- -1L
    }
    if (is.null(mz_max)) {
      mz_max <- -1L
    }
    cli_args <- c(cli_args, paste("/MzLimits", mz_min, mz_max))
  }

  #.. Output all hit lists, /All ...............................................
  cli_args <- c(cli_args, "/All")

  #.. Output file, /OUTTAB .....................................................
  cli_args <-
    c(cli_args, paste("/OUTTAB", normalizePath(output_file, mustWork = FALSE)))

  #.. Additional columns .......................................................
  # "/OutRevMF"    - Reverse Match Factor (the same as 'v')
  # "/OutRRevMF"   - PSS Match Factor
  # "/OutMW"       - Nominal MW
  # "/OutCAS"      - CAS registry number
  # "/OutChemForm" - Chemical formula
  # "/OutIK"       - InChIKey
  # "/RI x"        - Retention indices (I, S, L, Q, H searches)
  # "x"            - Do not output hit probabilities
  cli_args <- c(cli_args, addl_columns)
  temp <- unlist(strsplit(paste(addl_cli_args, collapse = " "), " ",
                          fixed = TRUE))
  if (all(trimws(temp) != "/RI")) {
    if (!is.null(ri_column_type)) {
      encoded_ri_column_type <- switch(ri_column_type,
                                       n = "n",
                                       s = "s",
                                       p = "p",
                                       stdnp = "n",
                                       semistdnp = "s",
                                       stdpolar = "p",
                                       stop("'ri_column_type' is invalid."))
      cli_args <- c(cli_args, paste0("/RI ", encoded_ri_column_type, "ux"))
    }
  } else {
    warning("'ri_column_type' is ignored, ",
            "because 'addl_cli_args' contains '/RI'.")
  }

  #.. Load libraries in memory .................................................
  if (load_in_memory) {
    cli_args <- c(cli_args, "/LibInMem")
  }

  #.. Additional arguments .....................................................
  cli_args <- c(cli_args, addl_cli_args)


  #--[ Check 'cli_args' ]-------------------------------------------------------

  temp <- vapply(strsplit(unlist(cli_args), " ", fixed = TRUE), function(a1) {
    a1[[1L]]
  }, character(1L))
  duplicated_args <- unique(temp[duplicated(temp)])
  mask <- !(duplicated_args == "/LIB")
  if (sum(mask) > 0L) {
    stop("Some CLI arguments are duplicated: ",
         paste(duplicated_args[mask], collapse = ", "))
  }


  #--[ Call MSPepSearch ]-------------------------------------------------------

  exe_path <- .GetExePath()
  if (!file.exists(exe_path)) {
    stop("'", basename(exe_path), "' is missing.")
  }

  if (.Platform$OS.type == "windows") {
    err_code <- system2(exe_path, args = unlist(cli_args), stderr = NULL)
  } else if (.Platform$OS.type == "unix") {
    err_code <- system2("wine", args = c(exe_path, unlist(cli_args)),
                        stderr = NULL)
  } else {
    stop("Only Windows or Unix-like OS (with Wine installed) are supported.")
  }

  if (err_code != 0L) {
    stop("'MSPepSearch' returned non-zero exit status: ", err_code, ".")
  }


  #--[ Output ]-----------------------------------------------------------------

  hitlists <- ParseOutputFile(output_file)
  return(invisible(hitlists))
}


