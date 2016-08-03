library(TCGAome)
context("Tests config.R functions")

test_that("configure_logging() in a non existing folder", {
  expect_equal(TCGAome::configure_logging("afolderthatdoesnotexist"),
               NULL)
  expect_warning(expect_error(futile.logger::flog.info("test logs")))
})

test_that("configure_logging() in an existing folder", {
  existing_folder = "."
  log_file = file.path(existing_folder, "TCGAome.log")
  expect_equal(TCGAome::configure_logging(existing_folder),
               NULL)
  expect_output(futile.logger::flog.info("test logs"))
  expect_true(file.exists(log_file))
  file.remove(log_file)
})

test_that("get_package_folder() in a non existing folder", {
  expect_is(TCGAome::get_package_folder("afolderthatdoesnotexist"), "character")
})

test_that("get_package_folder() in NULL", {
  expect_is(TCGAome::get_package_folder(NULL), "character")
})

test_that("get_results_folder()", {
  results_folder = TCGAome::get_results_folder()
  expect_is(results_folder, "character")
  expect_true(file.exists(results_folder))
  file.remove(results_folder)
})
