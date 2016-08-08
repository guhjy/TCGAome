context("Tests config.R functions")

test_that("configure_logging() in a non existing folder", {
    testthat::expect_equal(TCGAome::configure_logging("afolderthatdoesnotexist"),
               NULL)
    testthat::expect_warning(expect_error(futile.logger::flog.info("test logs")))
})

test_that("configure_logging() in an existing folder", {
    existing_folder = "."
    log_file = file.path(existing_folder, "TCGAome.log")
    testthat::expect_equal(TCGAome::configure_logging(existing_folder), NULL)
    testthat::expect_output(futile.logger::flog.info("test logs"))
    testthat::expect_true(file.exists(log_file))
    file.remove(log_file)
})

test_that("get_package_folder() in a non existing folder", {
    subfolder <- "afolderthatdoesnotexist"
    folder <- TCGAome::get_package_folder(subfolder)
    testthat::expect_is(folder, "character")
    testthat::expect_true(file.exists(folder))
    file.remove(folder)
    testthat::expect_false(file.exists(folder))
})

test_that("get_package_folder() in NULL", {
    folder <- TCGAome::get_package_folder(NULL)
    testthat::expect_is(folder, "character")
    testthat::expect_true(file.exists(folder))
})

test_that("get_results_folder()", {
    results_folder = TCGAome::get_results_folder()
    testthat::expect_is(results_folder, "character")
    testthat::expect_true(file.exists(results_folder))
    file.remove(results_folder)
})
