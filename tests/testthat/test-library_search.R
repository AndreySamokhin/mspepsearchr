#==[ Functions ]================================================================

CompareHitLists <- function(hl_mspepsearchr,
                            hl_mssearchr,
                            col_names = c("name", "mf", "rmf", "lib")) {

  all(vapply(seq_along(hl_mssearchr), function(a1) {

    # In the Similarity EI Hybrid algorithm, match factors of candidates from
    # hitlists returned by the MSPepSearch tool are sometimes zero. As a result,
    # these hitlists may differ from those obtained using MS Search.
    if (!is.null(hl_mspepsearchr[[a1]]$mf)) {
      mask <- (hl_mspepsearchr[[a1]]$mf > 0L)
    } else {
      mask <- (hl_mspepsearchr[[a1]]$score > 0L)
    }

    # When 'reverse_search = TRUE', the order of candidates may vary when the
    # match factors of adjacent candidates are identical. So, ordering
    # is required to align vectors.
    new_order <- match(hl_mspepsearchr[[a1]]$id[mask], hl_mssearchr[[a1]]$id)

    if (is.null(attr(hl_mspepsearchr[[a1]], "inlib"))) {
      attr(hl_mspepsearchr[[a1]], "inlib") <- NA_integer_
    }

    all(
      vapply(col_names, function(a2) {
        if (is.double(hl_mspepsearchr[[a1]][[a2]])) {
          if (a2 == "prob") {
            if (any(hl_mssearchr[[a1]][[a2]] < 0.0)) {
              return(TRUE)
            } else {
              # The probability (i.e., 'prob') is calculated differently in
              # MS Search 2.4 (2020) and MSPepSearch 0.9.7.1 (2023).
              return(all(abs(hl_mspepsearchr[[a1]][[a2]][mask] -
                               hl_mssearchr[[a1]][[a2]][new_order]) < 0.15))
            }
          } else {
            stop("'", a2, "' is double.")
          }
        } else {
          return(identical(hl_mspepsearchr[[a1]][[a2]][mask],
                           hl_mssearchr[[a1]][[a2]][new_order]))
        }
      }, logical(1L)),
      identical(attr(hl_mssearchr[[a1]], "unknown_name"),
                attr(hl_mspepsearchr[[a1]], "unknown_name")),
      identical(attr(hl_mssearchr[[a1]], "inlib"),
                attr(hl_mspepsearchr[[a1]], "inlib"))
    )
  }, logical(1L)))
}


UpdatePaths <- function(fun_args) {
  dir_names <- unlist(strsplit(fun_args$spectra, "/"))
  if (dir_names[[1L]] == "extdata") {
    fun_args$spectra <- system.file(fun_args$spectra,
                                    package = "mspepsearchr")
    fun_args$libraries <- system.file(fun_args$libraries,
                                      package = "mspepsearchr")
  } else if (dir_names[[1L]] == "data") {
    fun_args$spectra <- test_path(fun_args$spectra)
    fun_args$libraries <- test_path(fun_args$libraries)
  }
  return(fun_args)
}


PrintCallInfo <- function(fun_name,
                          fun_args) {
  call_expr <- as.call(c(as.name(fun_name), fun_args))
  call_str <- deparse1(call_expr, "\n", width.cutoff = 20)
  return(paste0("Check failed for:\n", call_str))
}



#==[ Temporary directory ]======================================================

temp_dir <-
  normalizePath(file.path(tempdir(), paste0("mspepsearchr_", Sys.getpid())),
                mustWork = FALSE)
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}
on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
temp_dir_arg <- list(temp_dir = temp_dir)



#==[ Snapshot-based tests ]=====================================================

#.. 'synthetic_spectra_results.rda' ...........................................#

load(test_path("data/rda/synthetic_spectra_results.rda"))
synthetic_spectra_results <- lapply(synthetic_spectra_results, function(res) {
  res$fun_args <- UpdatePaths(res$fun_args)
  return(res)
})

test_that(".LibrarySearch(), synthethic mass spectra", {
  for (res in synthetic_spectra_results) {
    hl <- do.call(res$fun_name, c(res$fun_args, temp_dir_arg))
    expect(CompareHitLists(hl, res$hl),
           PrintCallInfo(res$fun_name, res$fun_args))
  }
})


#.. 'alkanes_ei_lr_results.rda' ...............................................#

load(test_path("data/rda/alkanes_ei_lr_results.rda"))
alkanes_ei_lr_results <-
  lapply(alkanes_ei_lr_results, function(res) {
    res$fun_args <- UpdatePaths(res$fun_args)
    return(res)
  })

test_that(".LibrarySearch(), alkanes_ei_lr.msp", {
  for (res in alkanes_ei_lr_results) {
    hl <- do.call(res$fun_name, c(res$fun_args, temp_dir_arg))
    expect(CompareHitLists(hl, res$hl),
           PrintCallInfo(res$fun_name, res$fun_args))
  }
})


#.. 'alkyl_fragments_msms_hr_results.rda' .....................................#

load(test_path("data/rda/alkyl_fragments_msms_hr_results.rda"))
alkyl_fragments_msms_hr_results <-
  lapply(alkyl_fragments_msms_hr_results, function(res) {
    res$fun_args <- UpdatePaths(res$fun_args)
    return(res)
  })

test_that(".LibrarySearch(), alkanes_ei_lr.msp", {
  for (res in alkyl_fragments_msms_hr_results) {
    hl <- do.call(res$fun_name, c(res$fun_args, temp_dir_arg))
    expect(CompareHitLists(hl, res$hl),
           PrintCallInfo(res$fun_name, res$fun_args))
  }
})


#.. 'chlorine_compounds_ei_hr_results.rda' ....................................#

load(test_path("data/rda/chlorine_compounds_ei_hr_results.rda"))
chlorine_compounds_ei_hr_results <-
  lapply(chlorine_compounds_ei_hr_results, function(res) {
    res$fun_args <- UpdatePaths(res$fun_args)
    return(res)
  })

test_that(".LibrarySearch(), alkanes_ei_lr.msp", {
  for (res in chlorine_compounds_ei_hr_results) {
    hl <- do.call(res$fun_name, c(res$fun_args, temp_dir_arg))
    expect(CompareHitLists(hl, res$hl),
           PrintCallInfo(res$fun_name, res$fun_args))
  }
})


#.. 'chlorine_compounds_ei_lr_results.rda' ....................................#

load(test_path("data/rda/chlorine_compounds_ei_lr_results.rda"))
chlorine_compounds_ei_lr_results <-
  lapply(chlorine_compounds_ei_lr_results, function(res) {
    res$fun_args <- UpdatePaths(res$fun_args)
    return(res)
  })

test_that(".LibrarySearch(), alkanes_ei_lr.msp", {
  for (res in chlorine_compounds_ei_lr_results) {
    hl <- do.call(res$fun_name, c(res$fun_args, temp_dir_arg))
    expect(CompareHitLists(hl, res$hl),
           PrintCallInfo(res$fun_name, res$fun_args))
  }
})


#.. 'mw288_compounds_msms_hr_results.rda' .....................................#

load(test_path("data/rda/mw288_compounds_msms_hr_results.rda"))
mw288_compounds_msms_hr_results <-
  lapply(mw288_compounds_msms_hr_results, function(res) {
    res$fun_args <- UpdatePaths(res$fun_args)
    return(res)
  })

test_that(".LibrarySearch(), alkanes_ei_lr.msp", {
  for (res in mw288_compounds_msms_hr_results) {
    hl <- do.call(res$fun_name, c(res$fun_args, temp_dir_arg))
    expect(CompareHitLists(hl, res$hl),
           PrintCallInfo(res$fun_name, res$fun_args))
  }
})



#==[ Inline tests (01) ]========================================================

dummy_spectrum <- list(name = "dummy alkane",
                       mz = c(57, 71, 85),
                       intst = c(999, 500, 300))

lib_path <- system.file("extdata/libraries/massbank_subset_ei_lr/",
                        package = "mspepsearchr")


#.. 'presearch = list("mw", ...)' .............................................#

test_that(".LibrarySearch(..., presearch = list(\"mw\", ...))", {
  hl <- SimilaritySearchEiSimple(list(dummy_spectrum), lib_path,
                                 presearch = list("mw", 142),
                                 temp_dir = temp_dir)
  expect_true(nrow(hl[[1L]]) == 0L)

  hl <- SimilaritySearchEiSimple(list(dummy_spectrum), lib_path,
                                 presearch = list("mw", 156),
                                 temp_dir = temp_dir)
  expect_true(all(nrow(hl[[1L]]) > 0L, hl[[1L]]$mw == 156L))
})


#.. 'presearch = list("inchikey", ...)' .......................................#

dummy_spectra <-
  list(c(dummy_spectrum, inchikey = "RSJKGSCJYJTIGS-UHFFFAOYSA-N"),
       c(dummy_spectrum, inchikey = "RSJKGSCJYJTIGS-UHFFFAOYSA-X"),
       c(dummy_spectrum, inchikey = "RSJKGSCJYJTIGS-XXXXXXXXXX-X"),
       c(dummy_spectrum, inchikey = "XXXXXXXXXXXXXX-XXXXXXXXXX-X"))

test_that(".LibrarySearch(..., presearch = list(\"inchikey\", ...))", {
  hl <- SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                 presearch = list("inchikey", 1),
                                 temp_dir = temp_dir)
  expect_true(all(
    nrow(hl[[1L]]) > 0L, hl[[1L]]$inchikey == "RSJKGSCJYJTIGS-UHFFFAOYSA-N",
    nrow(hl[[2L]]) > 0L, hl[[2L]]$inchikey == "RSJKGSCJYJTIGS-UHFFFAOYSA-N",
    nrow(hl[[3L]]) > 0L, hl[[3L]]$inchikey == "RSJKGSCJYJTIGS-UHFFFAOYSA-N",
    nrow(hl[[4L]]) == 0L
  ))

  hl <- SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                 presearch = list("inchikey", 2),
                                 temp_dir = temp_dir)
  expect_true(all(
    nrow(hl[[1L]]) > 0L, hl[[1L]]$inchikey == "RSJKGSCJYJTIGS-UHFFFAOYSA-N",
    nrow(hl[[2L]]) > 0L, hl[[2L]]$inchikey == "RSJKGSCJYJTIGS-UHFFFAOYSA-N",
    nrow(hl[[3L]]) == 0L, nrow(hl[[4L]]) == 0L
  ))

  hl <- SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                 presearch = list("inchikey", 3),
                                 temp_dir = temp_dir)
  expect_true(all(
    nrow(hl[[1L]]) > 0L, hl[[1L]]$inchikey == "RSJKGSCJYJTIGS-UHFFFAOYSA-N",
    nrow(hl[[2L]]) == 0L, nrow(hl[[3L]]) == 0L, nrow(hl[[4L]]) == 0L
  ))
})


#.. 'n_hits' ..................................................................#

test_that(".LibrarySearch(..., n_hits = ...)", {
  hl <- SimilaritySearchEiSimple(list(dummy_spectrum), lib_path,
                                 n_hits = 3L, temp_dir = temp_dir)
  expect_true(nrow(hl[[1L]]) == 3L)
})



#==[ Inline tests (02) ]========================================================

dummy_spectra <- list(list(name = "dummy compound",
                           mz = c(57, 69, 99, 147, 149),
                           intst = c(999, 800, 700, 600, 500)))
lib_path <- system.file("extdata/libraries/massbank_subset_ei_lr_with_ri/",
                        package = "mspepsearchr")


#.. Retention indices .........................................................#

test_that(".LibrarySearch(..., ri_column_type = ..., ...))", {

  temp <- SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                   ri_column_type = "n")
  hl <- temp[[1L]]
  new_order <- order(hl$id)
  expect_identical(hl$ri[new_order],
                   c(1100L, NA_integer_, 1433L, NA_integer_, 1300L, 1300L,
                     1200L, 2296L, 1100L, 1200L))
  expect_identical(hl$ri_type[new_order],
                   c("a", "", "n", "", "a", "a", "a", "n", "a", "a"))

  temp <- SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                   ri_column_type = "s")
  hl <- temp[[1L]]
  new_order <- order(hl$id)
  expect_identical(hl$ri[new_order],
                   c(1100L, NA_integer_, 1434L, 1938L, 1300L, 1300L, 1200L,
                     2334L, 1100L, 1200L))
  expect_identical(hl$ri_type[new_order],
                   c("a", "", "s", "s", "a", "a", "a", "s", "a", "a"))

  temp <- SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                   ri_column_type = "p")
  hl <- temp[[1L]]
  new_order <- order(hl$id)
  expect_identical(hl$ri[new_order],
                   c(1100L, NA_integer_, NA_integer_, 2688L, 1300L, 1300L,
                     1200L, NA_integer_, 1100L, 1200L))
  expect_identical(hl$ri_type[new_order],
                   c("a", "", "", "p", "a", "a", "a", "", "a", "a"))

  temp <- SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                   ri_column_type = NULL)
  hl <- temp[[1L]]
  expect_null(hl$ri)
  expect_null(hl$ri_type)

  temp <- suppressWarnings(SimilaritySearchEiSimple(dummy_spectra, lib_path,
                                                    addl_cli_args = "/RI anux"))
  hl <- temp[[1L]]
  new_order <- order(hl$id)
  expect_identical(hl$ri[new_order],
                   c(1100L, NA_integer_, 1433L, 1938L, 1300L, 1300L, 1200L,
                     2296L, 1100L, 1200L))
  expect_identical(hl$ri_type[new_order],
                   c("a", "", "n", "s", "a", "a", "a", "n", "a", "a"))

})




