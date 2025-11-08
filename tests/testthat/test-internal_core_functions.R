#==[ Temporary directory ]======================================================

temp_dir <- CreateTempDir()
on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)



#==[ .SplitSpectraForParallel() ]===============================================

msp_input_file <- test_path("data/dummy_spectra/dummy01.msp")
mgf_input_file <- test_path("data/dummy_spectra/dummy01.mgf")
jdx_input_file <- test_path("data/dummy_spectra/dummy01.jdx")
msp_objs <- list(
  list(name = "001", mz = 101.0, intst = 999.0),
  list(name = "002", mz = 102.0, intst = 999.0),
  list(name = "003", mz = 103.0, intst = 999.0),
  list(name = "004", mz = 104.0, intst = 999.0),
  list(name = "005", mz = 105.0, intst = 999.0)
)


test_that(".SplitSpectraForParallel(), single thread", {
  # msp-file, 1 thread
  expect_identical(
    basename(.SplitSpectraForParallel(msp_input_file, temp_dir, 1L)),
    basename(msp_input_file)
  )

  # R object, 1 thread
  temp_file <- .SplitSpectraForParallel(msp_objs, temp_dir, 1L)
  expect_identical(readLines(temp_file), readLines(msp_input_file))
})


test_that(".SplitSpectraForParallel(), parallelization", {
  n_threads <- parallel::detectCores()
  skip_if(n_threads < 3L, "Skipping: fewer than 3 threads available.")

  # msp-file, 3 threads
  temp_files <- .SplitSpectraForParallel(msp_input_file, temp_dir, 3L)
  all_lines <- readLines(msp_input_file)
  expect_identical(length(temp_files), 3L)
  expect_identical(lapply(temp_files, readLines),
                   list(all_lines[1:8], all_lines[9:12], all_lines[13:20]))

  # mgf-file, 3 threads
  expect_error(.SplitSpectraForParallel(mgf_input_file, temp_dir, 3L),
               "Multiprocessing is not yet supported")

  # jdx-file, 3 threads
  expect_error(.SplitSpectraForParallel(jdx_input_file, temp_dir, 3L),
               "No spectra found")

  # R object, 3 threads
  temp_files <- .SplitSpectraForParallel(msp_objs, temp_dir, 3L)
  all_lines <- readLines(msp_input_file)
  expect_identical(length(temp_files), 3L)
  expect_identical(lapply(temp_files, readLines),
                   list(all_lines[1:8], all_lines[9:12], all_lines[13:20]))
})

