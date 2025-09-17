#==============================================================================#
#' Parse MSPepSearch output file
#'
#' @description
#'   Parse an output file containing the results of a library search performed
#'   with the MSPepSearch tool.
#'
#' @param input_file
#'   A string. The full path to a file containing the library search results
#'   (\code{/OUTTAB}).
#' @param rm_empty_cols
#'   A logical value. If \code{TRUE}, some empty columns are removed.
#' @param comments
#'   Any R object. Additional information (e.g., library search options, the
#'   list of libraries used, etc.). It is saved as the 'comments' attribute of
#'   the returned list.
#'
#' @return
#'   Return a list of data frames. Each data frame represents a hitlist.
#'
#' @importFrom utils read.table
#'
#' @noRd
#'
#==============================================================================#
ParseOutputFile <- function(input_file,
                            rm_empty_cols = TRUE,
                            comments = NULL) {

  #--[ Check input arguments ]--------------------------------------------------

  # 'input_file'
  if (!is.character(input_file) || length(input_file) != 1) {
    stop("'input_file' must be a string.")
  }
  if (!file.exists(input_file)) {
    stop("The file ", "'", input_file, "'", " does not exist.")
  }

  # 'rm_empty_cols'
  if (!is.logical(rm_empty_cols) || length(rm_empty_cols) != 1L) {
    stop("'rm_empty_cols' must be a logical value.")
  }

  # 'comments'
  # Any R object can be passed.


  #--[ Read lines ]-------------------------------------------------------------

  all_lines <- readLines(input_file)


  #--[ Parse comments ]---------------------------------------------------------

  meta_lines <- grep("^>", all_lines, value = TRUE)
  if (length(meta_lines) != 5L) {
    stop("Five comment lines were expected in the file.")
  }
  start_message <- substring(meta_lines[[1L]], 3L)
  call_cmd <- substring(meta_lines[[2L]], 3L)
  end_message <- substring(meta_lines[[5L]], 3L)


  #--[ Get column names ]-------------------------------------------------------

  header_idx <- min(which(!grepl("^>", all_lines)))
  raw_col_names <- unlist(strsplit(all_lines[[header_idx]], "\t", fixed = TRUE))
  clean_col_names <- gsub("\\W", "", gsub("[ .]", "_", tolower(raw_col_names)))
  if (all(c("name", "peptide") %in% clean_col_names)) {
    stop("Both 'name' and 'peptide' columns are present.")
  }
  name_map <-
    c(library = "lib", lib_mw = "mw", r_match = "rmf", pss_match = "pss_mf",
      dot_product = "dot", revdot = "rdot", pssdot = "pss_dot",
      mass = "exact_mass", deltamz = "delta_mz", peptide = "name",
      precursor_mz = "u_prec_mz", lib_precursor_mz = "prec_mz")
  final_col_names <- ifelse(clean_col_names %in% names(name_map),
                            name_map[clean_col_names], clean_col_names)

  # : Original          : Auto             : Manual     :
  # :-------------------:------------------:------------:
  # : Unknown           : unknown          :            :
  # : Rank              : rank             :            :
  # : Library           : library          : lib        :
  # : Id                : id               :            :
  # : Mass              : mass             : exact_mass :
  # : Lib MW            : lib_mw           : mw         :
  # : InLib             : inlib            :            :
  # : MF                : mf               : *          :
  # : R.Match           : r_match          : rmf        :
  # : PSS.Match         : pss_match        : pss_mf     :
  # : Score             : score            :            :
  # : DotProd           : dot_product      : dot        :
  # : RevDot            : revdot           : rdot       :
  # : PSS-Dot           : pssdot           : pss_dot    :
  # : Prob(%)           : prob             : *          :
  # : Name              : name             : *          :
  # : Charge            : charge           :            :
  # : Formula           : formula          : *          :
  # : FlankRes          : flankres         :            :
  # : Mods              : mods             :            :
  # : Pep               : pep              :            :
  # : CAS               : cas              : *          :
  # : Nreps             : nreps            :            :
  # : Tfratio           : tfratio          :            :
  # : Protein           : protein          :            :
  # : Peptide           : peptide          : name       :
  # : uInChIKey         : uinchikey        :            :
  # : InChIKey          : inchikey         : *          :
  # : RI                : ri               :            :
  # : Precursor m/z     : precursor_mz     : u_prec_mz  :
  # : Delta(m/z)        : deltamz          : delta_mz   :
  # : Lib Precursor m/z : lib_precursor_mz : prec_mz    :
  # : DeltaMass         : deltamass        :            :
  # : o.Match           : o_match          :            :
  # : o.R.Match         : o_r_match        :            :
  # : o.PSS.Match       : o_pss_match      :            :
  # : o.Score           : o_score          :            :
  # : o.DotProd         : o_dotprod        :            :


  #--[ Read table ]-------------------------------------------------------------

  integer_cols <- c("rank", "id", "mw", "inlib", "mf", "rmf", "pss_mf", "score",
                    "dot", "rdot", "pss_dot", "charge", "o_match", "o_r_match",
                    "o_pss_match", "o_score", "o_dotprod")
  numeric_cols <- c("exact_mass", "prob", "delta_mz", "prec_mz", "u_prec_mz",
                    "deltamass")
  col_classes <- rep("character", length(final_col_names))
  col_classes[final_col_names %in% c(integer_cols, numeric_cols)]  <- "numeric"

  # The 'fill = TRUE' option is used to correctly handle empty hitlists,
  # which appear as lines with a value only in the first column.
  df <- utils::read.table(input_file, sep = "\t", header = TRUE, fill = TRUE,
                          quote = "", colClasses = col_classes,
                          comment.char = ">")
  colnames(df) <- final_col_names

  # A stricter conversion to integers is implemented here because 'o.R.Match'
  # was sometimes exported as doubles (e.g., '7.00') instead of plain integers,
  # which caused 'utils::read.table()' to fail. This issue was observed in the
  # '@examples' section of the 'SimilaritySearchEiHybrid()' function.
  for (i in which(final_col_names %in% integer_cols)) {
    if (any(df[[i]] != floor(df[[i]]), na.rm = TRUE)) {
      stop("Column '", final_col_names[[i]],
           "' contains non-integer numeric values.")
    }
    df[[i]] <- as.integer(df[[i]])
  }


  #--[ Remove empty columns ]---------------------------------------------------

  if (rm_empty_cols) {
    if (all(is.na(df$charge))) {
      df$charge <- NULL
    }
    mask <- vapply(df, function(a1) {
      is.character(a1) && all(a1 == "")
    }, logical(1L))
    df <- df[!mask]
  }


  #--[ Parse RI column ]--------------------------------------------------------

  if (!is.null(df$ri)) {
    # It is supposed that RI is "<VALUE>-<TYPE>", "-", or "".
    temp <- vapply(strsplit(df$ri, "-", fixed = TRUE), function(a1) {
      if ((length(a1) == 1L && a1 == "") || length(a1) == 0L) {
        return(c(NA_character_, ""))
      } else {
        return(tolower(a1))
      }
    }, character(2L))
    df$ri <- as.integer(temp[1L, ])
    df$ri_type <- temp[2L, ]
  }


  #--[ Customize columns ]------------------------------------------------------

  # Reorder columns
  preferred_order <- c("name", "mf", "rmf", "pss_mf", "score", "dot", "rdot",
                       "pss_dot", "prob", "deltamass", "formula", "mw",
                       "exact_mass", "charge", "prec_type", "delta_mz",
                       "prec_mz", "inchikey", "cas", "lib", "id", "ri",
                       "ri_type")
  present_columns <- preferred_order[preferred_order %in% colnames(df)]
  remaining_columns <- setdiff(colnames(df), present_columns)
  df <- df[c(present_columns, remaining_columns)]


  #--[ Parse table ]------------------------------------------------------------

  cols_to_remove <- c("rank", "unknown", "inlib", "uinchikey", "u_prec_mz")
  hl_mask <- !(colnames(df) %in% cols_to_remove)

  # It is needed to empty hitlists and replace 'NA' values.
  df$rank[is.na(df$rank)] <- -1L

  # 'which(diff(df$rank) != 1L)' does not work for '/OutBestHitsOnly'
  temp <- which(diff(df$rank) < 1L)
  first_idxs <- c(1L, temp + 1L)
  last_idxs <- c(temp, nrow(df))

  hitlists <- lapply(seq_along(first_idxs), function(a1) {
    row_idxs <- seq(first_idxs[[a1]], last_idxs[[a1]])
    if (first_idxs[[a1]] == last_idxs[[a1]] && df$rank[[row_idxs]] == -1L) {
      hl <- df[0L, ]
    } else {
      hl <- df[row_idxs, hl_mask]
      # if (!is.null(hl$exact_mass)) {
      #   hl$exact_mass[hl$exact_mass == 0.0] <- NA_real_
      # }
      # if (!is.null(hl$cas)) {
      #   hl$cas[hl$cas == "0"] <- NA_character_
      # }
      # if (!is.null(hl$formula)) {
      #   hl$formula[hl$formula == ""] <- NA_character_
      # }
    }
    attr(hl, "unknown_name") <- df$unknown[[first_idxs[[a1]]]]
    if (!is.null(df$uinchikey)) {
      attr(hl, "unknown_inchikey") <- df$uinchikey[[first_idxs[[a1]]]]
    } else {
      attr(hl, "unknown_inchikey") <- NA_integer_
    }
    if (!is.null(df$inlib)) {
      attr(hl, "inlib") <- df$inlib[[first_idxs[[a1]]]]
    }
    if (!is.null(df$u_prec_mz)) {
      attr(hl, "prec_mz") <- df$u_prec_mz[[first_idxs[[a1]]]]
    }
    rownames(hl) <- seq_len(nrow(hl))
    return(hl)
  })

  attr(hitlists, "cli") <- list(start_message = start_message,
                                end_message = end_message,
                                call_cmd = call_cmd,
                                output_file = normalizePath(input_file))
  attr(hitlists, "comments") <- comments
  return(invisible(hitlists))
}


