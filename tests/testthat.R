Sys.setenv("R_TESTS" = "")

library("testthat")
library("VERSO")

test_check("VERSO")
