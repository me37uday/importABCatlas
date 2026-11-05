# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(importABCatlas)

# Keep the normal progress report while suppressing testthat's randomized
# celebratory praise, which is unrelated to test results and may contain emoji.
test_check(
  "importABCatlas",
  reporter = ProgressReporter$new(show_praise = FALSE)
)
