test_that("OpenHelpFile()", {
  local_mocked_bindings(.OpenTxtFile = function(...) { return(NULL) })
  expect_no_condition(OpenHelpFile())
  # local_mocked_bindings(file.exists = function(...) { return(FALSE) },
  #                       .package = "base")
  # expect_error(OpenHelpFile(), "Help file does not exist.")
})


