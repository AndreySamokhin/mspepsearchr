test_that("ParseOutputFile(), apostrophe as prime symbol", {
  hitlists <- ParseOutputFile(test_path("data/hitlists/prime_symbol.txt"))
  expect_identical(hitlists[[2L]]$name,
                   c("2,2'-bis...", "5,5'-dimethyl...", "N,N'-..."))
  expect_identical(hitlists[[2L]]$mf, c(800L, 700L, 600L))
})


test_that("ParseOutputFile(), quotation mark as double prime symbol", {
  hitlists <-
    ParseOutputFile(test_path("data/hitlists/double_prime_symbol.txt"))
  expect_identical(hitlists[[2L]]$name,
                   c("2,2'-bis...", "O,N',N\"-...", "N,N'-..."))
  expect_identical(hitlists[[2L]]$mf, c(801L, 701L, 601L))
})


test_that("ParseOutputFile(), empty hitlist", {
  hitlists <- ParseOutputFile(test_path("data/hitlists/empty_hitlist.txt"))
  expect_identical(nrow(hitlists[[1L]]), 0L)
})

